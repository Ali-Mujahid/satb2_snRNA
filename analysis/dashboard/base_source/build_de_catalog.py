#!/usr/bin/env python3
"""Package completed DESeq2 artifacts for a portable, lazy offline dashboard.

Reads produced tables/receipts only. No fits, raw-nucleus aggregation or altered
memberships. PRIMARY and requested sensitivity endpoints retain exact population
atoms and coordinate membership; arbitrary unions cannot acquire invented fits.
"""
from __future__ import annotations
import argparse
from collections import OrderedDict
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time

import numpy as np
import pandas as pd

BRANCHES=["original_full","paired_original","paired_cellbender"]
TOOLS=["scanpy","seurat"]
POPS=["allcells","Excitatory","Inhibitory","CA1_ProS"]
SEVEN=["gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
MIME="text/tab-separated-values;charset=utf-8"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for x in iter(lambda:f.read(8*1024*1024),b""):h.update(x)
    return h.hexdigest()

def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)):return [clean(v) for v in value]
    if isinstance(value,np.generic):return clean(value.item())
    if value is None or value is pd.NA:return None
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,Path):return str(value)
    return value

def dumps(value):return json.dumps(clean(value),ensure_ascii=True,separators=(",",":"),allow_nan=False)

def member_sha(ids):return hashlib.sha256("".join(str(i)+"\n" for i in sorted(ids)).encode()).hexdigest()

def records(frame):return clean(frame.to_dict("records"))

class Shards:
    def __init__(self,root):self.root=Path(root);self.index={};self.used={}
    def add(self,kind,payload,source_sha256=None,fixed_id=None):
        import base64
        raw=dumps(payload).encode();digest=hashlib.sha256(raw).hexdigest()
        ident=fixed_id or f"{kind}__{digest}"
        packed=gzip.compress(raw,compresslevel=6,mtime=0)
        encoded=base64.b64encode(packed).decode("ascii")
        entry={"encoding":"gzip-base64","data":encoded}
        script=";(function(w){var id="+json.dumps(ident)+",entry="+dumps(entry)+";(w.SATB2_SHARDS=w.SATB2_SHARDS||{})[id]=entry;if(typeof w.registerShard==='function'){w.registerShard(id,entry.data);}})(window);\n"
        path=self.root/"data/shards"/(ident+".js");path.parent.mkdir(parents=True,exist_ok=True)
        content=script.encode()
        if path.exists():
            if path.read_bytes()!=content:raise ValueError("Refusing different shard content at existing ID: "+ident)
        else:path.write_bytes(content)
        self.index[ident]={"path":str(path.relative_to(self.root)),"sha256":sha(path),"encoding":"gzip-base64","uncompressed_bytes":len(raw),"gzip_bytes":len(packed),"payload_sha256":digest,"source_sha256":source_sha256}
        return ident
    def download(self,e,path,role):
        path=e.bind(path);source_digest=e.hashes[str(path)]
        raw=path.read_bytes();text=raw.decode("utf-8")
        if hashlib.sha256(text.encode("utf-8")).hexdigest()!=source_digest:raise ValueError("Source download text did not preserve UTF-8 bytes")
        return self.add("download",{"text":text,"filename":f"{role}_{source_digest[:16]}.tsv","mime":MIME,"sha256":source_digest},source_digest)

def readiness(cfg):
    r=Path(cfg["results_root"])
    return {"complete_branches":[b for b in BRANCHES if (r/"deseq2"/b/"COMPLETE.json").exists()],"comparison_complete":(r/"comparison/COMPLETE.json").exists(),"count_validation_complete":(r/"count_validation/COMPLETE.json").exists()}

def streaming_members(e,path,index):
    path=e.bind(path);groups={p:[] for p in index.population_key}
    with gzip.open(path,"rt",encoding="utf-8",newline="") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["population_key"] not in groups:raise ValueError("Unregistered population membership")
            groups[row["population_key"]].append(row["nucleus_id"])
    for row in index.itertuples(index=False):
        ids=groups[row.population_key]
        if len(ids)!=len(set(ids)) or len(ids)!=int(row.n_nuclei) or member_sha(ids)!=row.membership_sha256:raise ValueError("Native population IDs/hash do not match")
    return groups

def coordinate_registry(e,path):
    x=e.table(path,dtype={k:str for k in ["dataset_id","tool","branch","population","nucleus_id"]})
    required={"dataset_id","tool","branch","population","nucleus_id","row_index"}
    if not required.issubset(x):raise ValueError("Coordinate index lacks required fields")
    out={}
    for ident,z in x.groupby("dataset_id",sort=False):
        if z[["tool","branch","population"]].drop_duplicates().shape[0]!=1:raise ValueError("Mixed coordinate dataset identity")
        a=z.row_index.to_numpy(float)
        if not np.isfinite(a).all() or not np.array_equal(a,np.floor(a)) or (a<0).any():raise ValueError("Coordinate row indexes are not nonnegative integers")
        z=z.assign(row_index=a.astype(np.int64)).sort_values("row_index")
        ids=z.nucleus_id.tolist()
        if ids!=sorted(ids) or len(ids)!=len(set(ids)) or not np.array_equal(z.row_index,np.arange(len(z))):raise ValueError("Coordinate row axis must be sorted unique IDs, zero-based, gap-free")
        tool,branch,pop=z[["tool","branch","population"]].iloc[0]
        if ident!=f"{tool}__{branch}__{pop}":raise ValueError("Unexpected canonical dataset ID")
        out[ident]={"id":ident,"tool":tool,"branch":branch,"population":pop,"ids":ids,"lookup":dict(zip(ids,range(len(ids))))}
    expected={f"{t}__{b}__{p}" for t in TOOLS for b in BRANCHES for p in POPS}
    if set(out)!=expected:raise ValueError("Coordinate index does not cover all24 exact datasets")
    return out

def settings(e,path):
    if not path.exists():return {}
    x=e.table(path,dtype=str)
    if x.setting.duplicated().any():raise ValueError("Duplicate model setting")
    return x.set_index("setting").value.to_dict()

def review_metadata(e,cfg,root,shards,coordinates,accept_final,bundle_root):
    review={"comparison_completed":accept_final,"statistical_implementation_independence":False,"same_shared_core_counts_scope":"COUNT_AND_SAME_DESEQ2_IMPLEMENTATION_REPRODUCIBILITY;NOT_INDEPENDENT_BIOLOGICAL_REPLICATION",
      "core_reliability_definition":"BH<0.05 with the same nonzero sign in all3 PRIMARY count branches; both original partition effects finite/tested/same sign; all12 original leave-one-mouse-out effects finite/tested/same sign. Sensitivity refits need not each have BH<0.05.",
      "fine_scope":"PRIMARY_ONLY;FINE_LABEL_CONFIDENCE_SEPARATE_FROM_MOUSE_MODEL_ELIGIBILITY;BH_WITHIN_EACH_ENDPOINT",
      "physical_layer_validated":False,"program_definition_scope":"WT_ONLY_ASSOCIATED_EXPRESSION;PRIMARY_THRESHOLD0.3_FIXED;ANATOMY_NOT_VALIDATED","review_documents":[]}
    core_by_dataset={}
    if accept_final:
        e.stage(root/"comparison","EXECUTED_REVIEW_REQUIRED","final_produced_artifact_comparison")
        tables={name:e.table(root/"comparison"/name) for name in ["EFFECT_COMPARISONS.tsv","FINE_POPULATION_CROSSWALK.tsv","SKIPPED_COMPARISONS.tsv","CORE_RELIABILITY_SCOPE.tsv","LABEL_RELIABILITY_REVIEW.tsv"]}
        effects=tables["EFFECT_COMPARISONS.tsv"]
        # Retain every assigned pair, but flag literal program-name mismatches.
        kind={}
        for branch in BRANCHES:
            ix=e.table(Path(cfg["deseq2"]["finalized_root_by_branch"][branch])/"POPULATION_INDEX.tsv",dtype={"population_key":str})
            kind.update({(branch,r.population_key):r.analysis_kind for r in ix.itertuples(index=False)})
        effects["different_program_labels"]=effects.apply(lambda r:kind.get((r.left_branch,r.left_population_key)) in ["layer_program","expression_program"] and str(r.left_analysis_label)!=str(r.right_analysis_label),axis=1)
        effects["interpretation_scope"]=np.where(effects.different_program_labels,"DIFFERENT_PROGRAM_LABELS_POPULATION_DEFINITION_SENSITIVITY",np.where(effects.comparison_scope.eq("SAME_SHARED_TOTAL_AND_SAME_DESEQ2_ENGINE"),"SAME_COUNTS_SAME_ENGINE_REPRODUCIBILITY","PRODUCED_FIXED_POPULATION_EFFECT_COMPARISON_WITH_OVERLAP_AND_LABEL_LIMITS"))
        for name,frame in tables.items():review[name.replace(".tsv","").lower()+"_shard"]=shards.add("review",{"rows":records(frame),"source_table":name},e.hashes[str((root/"comparison"/name).resolve())])
        core=e.table(root/"comparison/CORE_WT_cKO_RELIABILITY_ALL_GENES.tsv.gz",dtype={"gene_id":str,"tool":str,"population":str})
        symbol_map=e.table(root/"deseq2/original_full/AUDIT_GENE_ID_SYMBOL_MANIFEST.tsv",dtype={"gene_id":str,"gene_symbol":str})
        if symbol_map.gene_id.duplicated().any():raise ValueError("Audited original gene symbol axis is not unique")
        symbols=symbol_map.set_index("gene_id").gene_symbol.replace("",np.nan)
        if not core.gene_id.isin(symbols.index).all():raise ValueError("Core reliability gene absent from original symbol manifest")
        core.insert(0,"gene",core.gene_id.map(symbols))
        for (tool,pop),z in core.groupby(["tool","population"],sort=False):
            key=shards.add("core_reliability",{"rows":records(z),"scope":"ORIGINAL_FULL_CORE_TOTAL;ALL3_PRIMARY_AND_ORIGINAL_PARTITION_LOO_CHECKS","definition":review["core_reliability_definition"]})
            for branch in BRANCHES:core_by_dataset[f"{tool}__{branch}__{pop}"]=key
    supplement=Path(cfg["study_root"])/"reviews/CA1_programme_agreement_v2_0_0"
    if (supplement/"REVIEW_PROVENANCE.json").exists():
        e.json(supplement/"REVIEW_PROVENANCE.json")
        programme_tables={}
        for name in ["programme_semantic_agreement.tsv","same_named_programme_overlap.tsv"]:
            frame=e.table(supplement/name)
            programme_tables[name]=frame
            review[name.replace(".tsv","")+"_shard"]=shards.add("ca1_program_review",{"rows":records(frame),"source_table":name},e.hashes[str((supplement/name).resolve())])
        a=programme_tables["programme_semantic_agreement.tsv"]
        b=programme_tables["same_named_programme_overlap.tsv"]
        main_a=a[a.branch.eq("original_full")&np.isclose(a.leading_WT_score_threshold,.3)]
        main_b=b[b.branch.eq("original_full")&np.isclose(b.leading_WT_score_threshold,.3)]
        review["ca1_primary_same_name_summary"]={"overall":records(main_a),"programmes":records(main_b)}
        review["ca1_same_name_guard"]="Use the primary0.3 same-name reciprocal-overlap table. A same-name program failing50% reciprocal overlap is not an equivalent cross-tool population. Different-label/Hungarian matches are population-definition sensitivity, not equivalent anatomy. Threshold checks do not alter primary memberships."
        p=supplement/"CA1_programme_agreement_and_threshold_sensitivity.pdf"
        if p.exists():e.bind(p);review["review_documents"].append({"title":"CA1 program agreement and threshold sensitivity","source_path":str(p.resolve()),"source_sha256":sha(p)})
    asset_map=Path(bundle_root)/"data/asset_map.json"
    if asset_map.exists():
        assets=e.json(asset_map)
        for asset in assets:
            if "/diagnostics/" not in asset.get("bundle_path",""):continue
            portable=(Path(bundle_root)/asset["bundle_path"]).resolve(strict=True)
            if not portable.is_relative_to(Path(bundle_root).resolve()) or sha(e.bind(portable))!=asset["sha256"]:raise ValueError("Portable review PDF asset binding failed")
            doc={"title":Path(asset["bundle_path"]).stem.replace("_"," "),"path":asset["bundle_path"],"source_path":asset["source_path"],"sha256":asset["sha256"]}
            review["review_documents"]=[d for d in review["review_documents"] if d.get("source_path")!=asset["source_path"]]+[doc]
    return review,core_by_dataset

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config",required=True,type=Path);ap.add_argument("--bundle-root",required=True,type=Path)
    ap.add_argument("--dataset-index",required=True,type=Path);ap.add_argument("--allow-partial",action="store_true")
    ap.add_argument("--rebuild",action="store_true");ap.add_argument("--check",action="store_true");args=ap.parse_args()
    cfg=json.loads(args.config.read_text());state=readiness(cfg)
    if args.check:print(json.dumps({**state,"dataset_index_exists":args.dataset_index.exists()},indent=2));return
    complete=state["complete_branches"]
    if not state["count_validation_complete"] or (not args.allow_partial and (complete!=BRANCHES or not state["comparison_complete"])):raise ValueError("Required native DE/count/comparison producers are incomplete")
    if not complete:raise ValueError("No closed DE producer is available")
    out=args.bundle_root.resolve();out.mkdir(parents=True,exist_ok=True)
    catalog_path=out/"DE_CATALOG.json";catalog_js=out/"data/de_catalog.js";receipt_path=out/"DE_CATALOG_COMPLETE.json"
    if any(p.exists() for p in [catalog_path,catalog_js,receipt_path]) and not args.rebuild:raise ValueError("Refusing existing catalog; use explicit --rebuild for our offline artifact")
    root=Path(cfg["results_root"]);module_path=root.parent/"workflow/diagnostics/compare_revision.py"
    spec=importlib.util.spec_from_file_location("frozen_native_result_reader",module_path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # Packaging preserves the source decimal-to-Float64 roundtrip. This wrapper
    # changes only table reading for display; the frozen native reader and
    # scientific results are untouched, and downloads retain exact source text.
    class RoundTripEvidence(module.Evidence):
        def table(self,path,**kwargs):
            kwargs.setdefault("float_precision","round_trip")
            return super().table(path,**kwargs)
    e=RoundTripEvidence();e.bind(args.config);e.bind(module_path);e.bind(__file__)
    ac=e.stage(root/"count_validation","COMPLETE")
    if ac.get("exact_integer_counts_PASS") is not True or ac.get("n_differing_entries")!=0 or ac.get("config_sha256")!=sha(args.config):raise ValueError("Independent count acceptance/config binding failed")
    final=state["comparison_complete"] and complete==BRANCHES
    if final:
        cr=e.json(root/"comparison/COMPLETE.json")
        if cr.get("config_sha256")!=sha(args.config) or cr.get("source_code_sha256")!=sha(module_path) or cr.get("exact_integer_counts_PASS") is not True or cr.get("same_shared_total_counts_PASS") is not True:raise ValueError("Final comparator acceptance/source binding failed")
    coords=coordinate_registry(e,args.dataset_index);shards=Shards(out)
    atom_index=e.table(Path(cfg["initial_annotation_dir"])/"atomic_population_index.tsv",dtype={"population_key":str,"cluster_id":str})
    atomic={(r.branch,r.population_key):r for r in atom_index.itertuples(index=False)}
    data=module.Results.__new__(module.Results);data.e=e;data.cache=OrderedDict();data.endpoint_audits=[];data.contracts=[];data.model_roots={};data.rosters={}
    datasets=[];endpoint_navigation=[];sources={};all_origins={};indices={};fit_sets={key:[] for key in coords}
    for branch in BRANCHES:
        finalized=Path(cfg["deseq2"]["finalized_root_by_branch"][branch]);e.stage(finalized,"EXECUTED_REVIEW_REQUIRED","finalize_populations",branch)
        ip=e.bind(finalized/"POPULATION_INDEX.tsv");atp=e.bind(finalized/"POPULATION_ATOMS.tsv")
        if ac.get("validated_sources",{}).get(str(ip))!=sha(ip) or ac.get("validated_sources",{}).get(str(atp))!=sha(atp):raise ValueError("Count acceptance does not bind consumed final populations/atoms")
        ix=e.table(ip,dtype={"population_key":str,"membership_sha256":str,"analysis_label":str,"category_family":str});indices[branch]=ix
        mi=streaming_members(e,finalized/"POPULATION_MEMBERSHIPS.tsv.gz",ix)
        atoms=e.table(atp,dtype=str).groupby("population_key",sort=False).atomic_population_key.apply(list).to_dict()
        sample=e.table(cfg["deseq2"]["sample_manifest_tsv_by_branch"][branch],dtype={"sample":str,"condition":str})
        if len(sample)!=12 or sample["sample"].duplicated().any() or sample.condition.value_counts().to_dict()!={"WT":6,"cKO":6}:raise ValueError("Incorrect all12 mouse roster")
        data.rosters[branch]=dict(zip(sample["sample"],sample.condition));data.model_roots[branch]=root/"deseq2"/branch
        if branch in complete:
            dr=e.stage(data.model_roots[branch],"EXECUTED_REVIEW_REQUIRED","deseq2",branch)
            master=e.table(data.model_roots[branch]/"SOURCE_AND_MASTER_HASHES.tsv",dtype=str)
            matches=master[master.path.str.endswith("/config/analysis_v2_0_0.json")]
            if len(matches)!=1 or matches.sha256.iloc[0]!=sha(args.config):raise ValueError("DE producer master binding differs")
            for guard in master.itertuples(index=False):
                if sha(e.bind(guard.path))!=guard.sha256:raise ValueError("DE producer source/master hash differs")
            registry=e.table(data.model_roots[branch]/"PSEUDOBULK_EXECUTION_AUDIT.tsv",dtype={"population_key":str})
            primary=registry[registry.endpoint.eq("PRIMARY")]
            if primary.population_key.duplicated().any() or set(primary.population_key)!=set(ix.population_key) or registry.status.eq("ERROR").any():raise ValueError("Incomplete/failed completed model registry")
        for tool in TOOLS:
            for pop in POPS:
                ident=f"{tool}__{branch}__{pop}";baseline=coords[ident];tot=ix[ix.tool.eq(tool)&ix.population.eq(pop)&ix.analysis_kind.eq("total")]
                if len(tot)!=1 or set(baseline["ids"])!=set(mi[tot.iloc[0].population_key]):raise ValueError("Coordinate baseline is not exact produced TOTAL membership")
                datasets.append({"id":ident,"tool":tool,"branch":branch,"population":pop,"cluster_scope":"original_frozen","n_nuclei":len(baseline["ids"]),"coordinate_membership_sha256":member_sha(baseline["ids"]),"fits_shard":f"fits__{ident}","producer_status":"COMPLETE" if branch in complete else "PENDING_NATIVE_PRODUCER"})
        for row in ix.itertuples(index=False):
            ident=f"{row.tool}__{branch}__{row.population}";lookup=coords[ident]["lookup"]
            atom_keys=sorted(atoms.get(row.population_key,[]))
            atom_rows=[atomic[(branch,a)] for a in atom_keys]
            if any(a.tool!=row.tool or a.population!=row.population for a in atom_rows):raise ValueError("Population atoms cross dataset identity")
            cluster_ids=sorted({str(a.cluster_id) for a in atom_rows})
            ids=mi[row.population_key]
            if not set(ids).issubset(lookup):raise ValueError("Fit membership escapes coordinate dataset")
            members=sorted(lookup[i] for i in ids)
            if sha(e.bind(row.counts_file))!=row.counts_sha256:raise ValueError("Produced model count source hash differs")
            labels=["PRIMARY"]
            if branch in complete and branch=="original_full" and str(row.core_total).lower() in ["true","1"]:labels += ["KNOWN_PARTITION_ADJUSTED","MAIN_PARTITION_ONLY"]+["LOO_omit_"+s for s in data.rosters[branch]]
            for label in labels:
                fid=f"{branch}::{row.tool}::{row.population_key}::{label}"
                fit={"id":fid,"dataset_id":ident,"tool":row.tool,"branch":branch,"population":row.population,"cluster_scope":"original_frozen","population_key":row.population_key,"endpoint":label,"analysis_kind":row.analysis_kind,"category_family":row.category_family,"analysis_label":row.analysis_label,"membership_sha256":row.membership_sha256,"n_nuclei":int(row.n_nuclei),"cluster_ids":cluster_ids,"atom_keys":atom_keys,"members":members,"counts_source_sha256":row.counts_sha256,"count_branch":branch,"numerator":"cKO","denominator":"WT","model_enabled":str(row.model_enabled).lower() in ["true","1"],"rows_shard":None,"download_shard":None,"stable_audit_shard":None,"coverage":[],"model_settings":{},"native_source_paths":{"counts":row.counts_file}}
                if branch not in complete:fit.update(status="PENDING_NATIVE_PRODUCER",reason="No completed DE producer receipt; partial endpoint files are not exposed")
                else:
                    result,ad=data.endpoint(branch,row,label);ep=data.model_roots[branch]/row.tool/row.population_key/label
                    fit.update(clean(ad));fit["model_n_nuclei"]=fit.get("n_nuclei");fit["n_nuclei"]=int(row.n_nuclei)
                    fit["model_endpoint"]=ad.get("endpoint");fit["endpoint"]=label
                    fit["native_source_paths"]["endpoint_audit"]=str(ep/"AUDIT/EXECUTION_AUDIT.tsv")
                    cov=e.table(ep/"AUDIT/MICE_ALL12_AUDIT.tsv",dtype={"sample":str,"condition":str})
                    if len(cov)!=12 or set(cov["sample"])!=set(data.rosters[branch]) or cov["sample"].duplicated().any() or (cov.condition!=cov["sample"].map(data.rosters[branch])).any():raise ValueError("Endpoint coverage animal/genotype identities differ")
                    if label=="PRIMARY":
                        minimum=1 if row.population=="Inhibitory" else 20
                        if not np.array_equal(module.truth(cov.included_in_model),cov.n_nuclei.ge(minimum)&cov.total_UMIs.gt(0)):raise ValueError("PRIMARY mouse inclusion does not follow declared coverage policy")
                    fit["coverage"]=records(cov);fit["model_settings"]=settings(e,ep/"AUDIT/MODEL_SETTINGS.tsv")
                    fit["source_coverage_shard"]=shards.download(e,ep/"AUDIT/MICE_ALL12_AUDIT.tsv","all12_mouse_audit")
                    pf=e.table(ep.parent/"AUDIT/FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv",dtype={"gene_id":str})
                    fit["gene_universe"]={"count_cutoff_per_mouse":int(pf.count_cutoff_per_mouse.iloc[0]),"required_mice":int(pf.required_mice.iloc[0]),"frozen_gene_count":int(module.truth(pf.passes_primary_prefilter).sum()),"policy":str(pf.universe_policy.iloc[0]),"audit_shard":shards.download(e,ep.parent/"AUDIT/FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv","frozen_gene_prefilter")}
                    warning=ep/"AUDIT/MODEL_WARNINGS.txt";fit_type=ep/"AUDIT/DISPERSION_FIT_TYPE.txt"
                    fit["warning_text"]=e.bind(warning).read_text().strip() if warning.exists() else ""
                    fit["actual_dispersion_fit_type"]=e.bind(fit_type).read_text().strip() if fit_type.exists() else None
                    if result is not None:
                        # Exact source text is separate from parsed display/IDs.
                        fit["download_shard"]=shards.download(e,ep/"ALL_RESULTS.txt","DESeq2_standard_seven_columns")
                        fit["stable_audit_shard"]=shards.download(e,ep/"AUDIT/ALL_RESULTS_STABLE_ID_AUDIT.tsv","DESeq2_stable_ID_audit")
                        rows=result[["gene_id",*SEVEN[1:]]].copy();rows.insert(1,"gene",result.display_symbol)
                        fit["rows_shard"]=shards.add("de_rows",{"rows":records(rows),"columns":["gene_id",*SEVEN],"stable_ID_binding":"Exact native export-row audit; duplicate/missing symbols retained","NA_policy":"null=unavailable; never nonsignificant"})
                        fit["display_download_filename"]=f"{branch}_{row.tool}_{row.population_key}_{label}_DESeq2.tsv"
                        fit["native_source_paths"].update(results=str(ep/"ALL_RESULTS.txt"),stable_audit=str(ep/"AUDIT/ALL_RESULTS_STABLE_ID_AUDIT.tsv"))
                    all_origins[ad.get("execution_origin","UNKNOWN")]=all_origins.get(ad.get("execution_origin","UNKNOWN"),0)+1
                fit["inference_caveats"]=["BH is within this endpoint, not across overlapping population families","Low-nucleus IN support is exploratory and does not validate cell identity","Same-count shared totals use one DESeq2 implementation; equality is not independent replication","Expression-associated layer/program labels do not establish physical anatomy","An arbitrary selection has inference only when count branch, atoms and exact membership match a produced fit"]
                fit_sets[ident].append(fit)
                endpoint_navigation.append({k:fit.get(k) for k in ["id","dataset_id","endpoint","analysis_kind","category_family","analysis_label","status","reason","n_WT","n_cKO","residual_df","n_endpoint_genes","n_pvalue","n_padj","DEG_FDR05","coverage_flag","warning_count","execution_origin"]})
    review,core_links=review_metadata(e,cfg,root,shards,coords,final,out)
    for dataset in datasets:
        if dataset["id"] in core_links:dataset["core_reliability_shard"]=core_links[dataset["id"]]
        payload={"fits":fit_sets[dataset["id"]],"dataset_id":dataset["id"],"coordinate_membership_sha256":dataset["coordinate_membership_sha256"],"producer_status":dataset["producer_status"]}
        version=hashlib.sha256(dumps(payload).encode()).hexdigest()[:16]
        dataset["fits_shard"]=shards.add("fits",payload,fixed_id=dataset["fits_shard"]+"__"+version)
    e.finish()
    catalog={"schema_version":"satb2-offline-1","status":"COMPLETE" if final else "PARTIAL_COMPLETED_PRODUCERS_ONLY","datasets":datasets,"shards":shards.index,"review_metadata":review,"endpoint_navigation":endpoint_navigation,"default_contrast":{"numerator":"cKO","denominator":"WT"},"exact_selection_contract":{"required":["dataset_identity","actual_count_branch","exact_atom_union","exact_sorted_coordinate_members"],"no_match_status":"NO MATCHING FIT","no_match_request":{"operation":"request_native_pseudobulk_fit","scientific_execution_performed":False,"required_fields":["tool","branch","population","cluster_scope","selected_cluster_ids","selected_atom_keys","selected_coordinate_members","membership_sha256","numerator","denominator"]}},"config_sha256":sha(args.config),"dataset_index_sha256":sha(args.dataset_index),"native_comparison_completed":final,"scientific_fitting_performed":False,"raw_nucleus_count_aggregation_performed":False}
    archive=None;old_receipt={};archived_orphans=[]
    if args.rebuild and catalog_path.exists():
        # Preserve the prototype outside the deliverable, then prune only files
        # owned and byte-bound by its own receipt. Scientific/base assets never
        # enter this operation, and all current source/result shards remain.
        archive=out.parent/"offline_dashboard_catalog_history"/(time.strftime("%Y%m%dT%H%M%SZ",time.gmtime())+"__"+sha(catalog_path)[:12])
        archive.mkdir(parents=True,exist_ok=False)
        for p in [catalog_path,catalog_js,receipt_path]:
            if p.exists():
                dst=archive/p.relative_to(out);dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(p.read_bytes())
        if receipt_path.exists():old_receipt=json.loads(receipt_path.read_text())
    catalog_path.write_text(dumps(catalog)+"\n");catalog_js.parent.mkdir(parents=True,exist_ok=True)
    catalog_js.write_text("window.SATB2_DE_CATALOG="+dumps(catalog)+";\n")
    e.finish()
    outputs={str(catalog_path.relative_to(out)):sha(catalog_path),str(catalog_js.relative_to(out)):sha(catalog_js)}
    outputs.update({v["path"]:v["sha256"] for v in shards.index.values()})
    if final and archive is not None:
        for relative,digest in old_receipt.get("output_sha256",{}).items():
            if not relative.startswith("data/shards/") or relative in outputs:continue
            p=(out/relative).resolve(strict=True)
            if not p.is_relative_to(out) or p.is_symlink() or sha(p)!=digest:raise ValueError("Owned prototype shard changed; refusing archive prune")
            dst=archive/relative;dst.parent.mkdir(parents=True,exist_ok=True)
            count=p.stat().st_size;p.rename(dst)
            if sha(dst)!=digest:raise ValueError("Prototype shard archive digest differs")
            archived_orphans.append({"path":relative,"sha256":digest,"bytes":count})
    receipt={"status":"COMPLETE" if final else "PARTIAL","stage":"offline_completed_DE_catalog","schema_version":"satb2-offline-1","config_sha256":sha(args.config),"source_code_sha256":sha(__file__),"native_reader_sha256":sha(module_path),"input_sha256":e.hashes,"output_sha256":outputs,"n_datasets":len(datasets),"n_endpoints":len(endpoint_navigation),"n_lazy_shards":len(shards.index),"execution_origins":all_origins,"scientific_fitting_performed":False,"raw_nucleus_count_aggregation_performed":False,"native_outputs_modified":False,"exact_source_downloads":True,"numeric_parsing":"source decimal-to-Float64 round_trip","duplicate_symbols_preserved":True,"unavailable_FDR_policy":"null;not_nonsignificant","native_comparator_completed":final,"previous_catalog_archive":str(archive) if archive else None,"owned_orphan_shards_archived":archived_orphans,"owned_orphan_bytes_archived":sum(v["bytes"] for v in archived_orphans)}
    receipt_path.write_text(dumps(receipt)+"\n")
    print(json.dumps({k:receipt[k] for k in ["status","n_datasets","n_endpoints","n_lazy_shards","scientific_fitting_performed"]},indent=2))

if __name__=="__main__":main()
