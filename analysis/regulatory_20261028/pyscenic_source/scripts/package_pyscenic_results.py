#!/usr/bin/env python3
"""Seal completed positive motif-supported pySCENIC/AUCell and mouse association outputs.
No fitting, scoring, averaging, differential expression or FDR adjustment in this adapter.
Run only in the immutable native workflow; writes a fresh atlas and FINAL_VALIDATION last.
"""
from __future__ import annotations
import argparse,base64,csv,gzip,hashlib,json,math,re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
SCHEMA="satb2-regulatory-atlas/v1"
ENDPOINTS={
 "CA1_Deep":("CA1_DL","P0244_seurat_CA1_ProS_LAYER_PROGRAM",8750,"CA1 deep-associated expression program"),
 "CA1_Superficial":("CA1_SL","P0247_seurat_CA1_ProS_LAYER_PROGRAM",11706,"CA1 superficial-associated expression program"),
 "IN_C1":("IN_C1","P0140_seurat_Inhibitory_CLUSTER",1718,"Inhibitory C1"),
 "IN_C6":("IN_C6","P0145_seurat_Inhibitory_CLUSTER",501,"Inhibitory C6"),
 "IN_C8":("IN_C8","P0147_seurat_Inhibitory_CLUSTER",468,"Inhibitory C8 · mixed identity review")}
REQUIRED_NATIVE=["CA1/prepared/TF_CANDIDATE_COVERAGE.tsv","IN/prepared/TF_CANDIDATE_COVERAGE.tsv","REGULON_ACTIVITY_PER_MOUSE.tsv","REGULON_COVERAGE.tsv","REGULON_TARGETS.tsv","LINEAGE_STATUS.tsv","PROVENANCE.json",
 "SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv","SATB2_NO_SELF_TARGET_COVERAGE.tsv"]
REQUIRED_MOUSE=["REGULON_ASSOCIATIONS.tsv","REGULON_SAME_PARTITION_ASSOCIATIONS.tsv","REGULON_LOO_MODELS.tsv","REGULON_LOO_STABILITY.tsv",
 "REGULON_MOUSE_BOOTSTRAP.tsv","REGULON_MODEL_ELIGIBILITY.tsv","REGULON_MOUSE_SCORES.tsv","MODEL_SUMMARY.tsv","SATB2_STATUS.tsv","PROVENANCE.json",
 "SATB2_NO_SELF_TARGET_ASSOCIATIONS.tsv","SATB2_NO_SELF_TARGET_SAME_PARTITION_ASSOCIATIONS.tsv","SATB2_NO_SELF_TARGET_LOO_MODELS.tsv",
 "SATB2_NO_SELF_TARGET_LOO_STABILITY.tsv","SATB2_NO_SELF_TARGET_MOUSE_BOOTSTRAP.tsv","SATB2_NO_SELF_TARGET_MOUSE_SCORES.tsv"]
def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def receipt(path, root=None):
    path = Path(path).resolve()
    return {"path": str(path.relative_to(Path(root).resolve())) if root else str(path),
            "bytes": path.stat().st_size, "sha256": digest(path)}


def raw_json(value, compact=False):
    return (json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True,
                       separators=(",", ":") if compact else None, indent=None if compact else 2)+"\n").encode()


def load(path):
    return json.loads(Path(path).read_text())


def table(path, columns=()):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)), f"Missing/duplicate TSV header: {path}")
        require(set(columns) <= set(reader.fieldnames), f"Missing native columns {set(columns)-set(reader.fieldnames)}: {path}")
        rows = list(reader)
        require(all(None not in row for row in rows), f"Malformed native TSV: {path}")
    return rows


def num(value, required=False):
    if value is None or str(value) in ("", "NA", "NaN", "nan", "null", "None"):
        require(not required, "Required native numeric value is missing")
        return None
    value = float(value)
    require(math.isfinite(value), "Native numeric value is nonfinite")
    return value


def integer(value):
    value = num(value, True)
    require(value >= 0 and value.is_integer(), "Native count is not a nonnegative integer")
    return int(value)


def boolean(value):
    if value in (True, "TRUE", "True", "true", "1", 1):
        return True
    require(value in (False, "FALSE", "False", "false", "0", 0), f"Native Boolean missing/invalid: {value}")
    return False


def same(a, b):
    return a == b or a is not None and b is not None and math.isclose(a, b, rel_tol=1e-13, abs_tol=1e-14)


def child(root, relative):
    root = Path(root).resolve()
    path = (root/relative).resolve()
    require(path.is_relative_to(root), f"Native receipt escapes its result root: {relative}")
    require(not (root/relative).is_symlink(), f"Symlink is not a native result: {relative}")
    return path


def check_record(record):
    require(isinstance(record, dict) and {"path", "bytes", "sha256"} <= set(record), "Input lacks a complete hash binding")
    path = Path(record["path"])
    require(path.is_absolute() and path.is_file() and not path.is_symlink(), f"Unavailable bound input: {path}")
    require(re.fullmatch(r"[a-f0-9]{64}", record["sha256"]) is not None, "Invalid input SHA256")
    actual = receipt(path)
    require(actual["bytes"] == record["bytes"] and actual["sha256"] == record["sha256"], f"Bound input changed: {path}")
    return actual


def sealed_step(folder, schema, required, cfg_record, source_path):
    require(not (folder/"FAILED.json").exists(), f"Native failure marker exists: {folder}")
    marker_path = folder/"COMPLETE.json"
    marker = load(marker_path)
    require(marker.get("schema") == schema and marker.get("status") == "COMPLETE", f"Native step is not complete: {folder}")
    require(check_record(marker["config"])["sha256"] == cfg_record["sha256"], "Native step configuration differs")
    source = check_record(marker["source"])
    require(source["sha256"] == digest(source_path), f"Executed source differs: {source_path}")
    bound = marker.get("output_files", marker.get("outputs", {}))
    require(bound and set(required) <= set(bound), f"Native step receipt does not cover required tables: {folder}")
    for relative, record in bound.items():
        path = child(folder, relative)
        require(path.is_file() and path.stat().st_size == record["bytes"] and digest(path) == record["sha256"], f"Native output changed: {path}")
    for record in marker.get("inputs", {}).values():
        if isinstance(record, dict) and {"path", "bytes", "sha256"} <= set(record):
            check_record(record)
    # No extra unsealed output may enter the archive or UI silently.
    found = {str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file()}
    require(found == set(bound)|{"COMPLETE.json"}, f"Unsealed native output files: {found-set(bound)-{'COMPLETE.json'}}")
    return marker


def frozen_cohort(cfg):
    extra = cfg["additional_hash_inputs"]
    audit = table(extra["frozen_membership_audit"]["path"], ["endpoint", "population_key", "retained_nuclei", "membership_sha256", "native_population_membership_sha256"])
    frozen = {r["endpoint"]: r for r in audit}
    require(len(frozen) == len(audit) == 5 and set(frozen) == set(ENDPOINTS), "Frozen cohort does not have exactly five endpoints")
    for ep, r in frozen.items():
        require(r["population_key"] == ENDPOINTS[ep][1] and integer(r["retained_nuclei"]) == ENDPOINTS[ep][2], f"Frozen endpoint identity/count differs: {ep}")
        require(re.fullmatch(r"[a-f0-9]{64}", r["membership_sha256"]) and r["membership_sha256"] == r["native_population_membership_sha256"], f"Membership hashes disagree: {ep}")
    sample_rows = table(cfg["samples"]["path"], ["sample", "condition"])
    samples = {r["sample"]: r["condition"] for r in sample_rows}
    require(len(samples) == len(sample_rows) == 12 and Counter(samples.values()) == {"WT": 6, "cKO": 6}, "Mouse roster is not six WT + six cKO")
    rows = table(extra["mouse_coverage"]["path"], ["endpoint", "sample", "condition", "n_nuclei"])
    coverage = {(r["endpoint"], r["sample"]): r for r in rows}
    require(len(coverage) == len(rows) == 60 and set(coverage) == {(ep, s) for ep in ENDPOINTS for s in samples}, "Frozen coverage does not contain 60 exact endpoint-mouse rows")
    for (ep, sample), r in coverage.items():
        require(r["condition"] == samples[sample] and integer(r["n_nuclei"]) > 0, "Frozen mouse genotype/nucleus coverage differs")
    for ep in ENDPOINTS:
        require(sum(integer(coverage[ep, s]["n_nuclei"]) for s in samples) == ENDPOINTS[ep][2], "Frozen mouse counts do not sum to endpoint")
    return frozen, samples, coverage



def lineage(ep):
    require(ep in ENDPOINTS,"Unknown frozen endpoint")
    return "CA1" if ep.startswith("CA1_") else "IN"

def verify_inputs(cfg):
    require(cfg.get("schema")=="satb2-pyscenic/v1","Unexpected pySCENIC configuration schema")
    require(cfg.get("min_regulon_targets")==10,"Protocol fixes ten measured regulon targets")
    require({"samples","prepared_complete","lineages","references","additional_hash_inputs"}<=set(cfg),"Missing bound scientific inputs")
    require(set(cfg["lineages"])=={"CA1","IN"},"Need exact CA1/IN input maps")
    for lin in ("CA1","IN"):
        require({"counts","genes","cells"}<=set(cfg["lineages"][lin]),"Lineage input map incomplete")
    require({"tf_names","motif_annotations","rankings"}<=set(cfg["references"]),"Missing exact matched motif references")
    require(len(cfg["references"]["rankings"])==2,"Require two approved gene-based mm10 rankings")
    require({"frozen_membership_audit","mouse_coverage"}<=set(cfg["additional_hash_inputs"]),"Missing frozen cohort evidence")
    found={}
    def walk(value,key):
        if isinstance(value,dict):
            if {"path","bytes","sha256"}<=set(value):found[key]=check_record(value)
            else:
                for k,v in value.items():walk(v,key+"/"+k)
        elif isinstance(value,list):
            for i,v in enumerate(value):walk(v,key+"/"+str(i))
    walk(cfg,"")
    require(found,"No immutable input records")
    return found

def probability(value):
    out=num(value)
    require(out is None or 0<=out<=1,"Native P/FDR or AUCell AUC outside[0,1]")
    return out

def method(lin):return "pyscenic_motif_"+lin

def feature(lin,reg):return method(lin)+"::"+reg

def node_id(ep,reg,symbol):return "ps_"+hashlib.sha256(raw_json([ep,reg,symbol],True)).hexdigest()[:24]

def identity(r):
    ep=r["endpoint"];lin=lineage(ep)
    require(r["lineage"]==lin and r["regulon_id"].startswith(lin+"::"),"Cross-lineage regulon alias")
    require(r["endpoint_id"]==ENDPOINTS[ep][0] and r["source_endpoint_id"]==ENDPOINTS[ep][1],"Changed exact endpoint identity")
    require(r["native_name"].endswith("(+)") and r["TF"],"Not a native positive motif-supported regulon")
    return ep,lin,r["regulon_id"]

def validate_score_maps(score_rows,cov,samples,source_coverage):
    values={}
    expected={(ep,reg,s) for (ep,reg),c in cov.items() if boolean(c["eligible"]) for s in samples}
    for r in score_rows:
        ep,lin,reg=identity(r);key=ep,reg,r["sample"]
        require(key not in values and (ep,reg) in cov and key in expected,"Duplicated or unmodelable native AUCell mouse score")
        c=cov[ep,reg];s=r["sample"]
        require(r["condition"]==samples[s] and integer(r["n_nuclei"])==integer(source_coverage[ep,s]["n_nuclei"]),"Wrong frozen mouse membership/count")
        require(r["score_unit"]=="native_AUCell_AUC" and boolean(r["eligible"]),"Wrong score units/eligibility")
        require(probability(r["aucell_mean"]) is not None,"Missing native mean AUCell")
        for k in ("native_name","TF","n_targets","n_background_genes"):
            require(r[k]==c[k],f"Native AUCell score/coverage disagreement: {k}")
        values[key]=r
    require(set(values)==expected,"Native per-mouse regulon scores lack exact twelve-mouse coverage")
    return values

def associations(path,cov,values,samples,source_table):
    rows=table(path,["lineage","endpoint","regulon_id","native_name","TF","tested","effect_cKO_minus_WT","pvalue","association_fdr","n_WT","n_cKO","fdr_family"])
    require(len(rows)==len(cov),"Native association/coverage row counts disagree")
    effects=[];seen=set()
    for r in rows:
        ep,lin,reg=identity(r);key=ep,reg
        require(key in cov and key not in seen,"Unknown or duplicate native regulon association");seen.add(key)
        tested=boolean(r["tested"])
        require(not tested or boolean(cov[key]["eligible"]),"Unavailable AUC association tested")
        require(not tested or (integer(r["n_WT"]),integer(r["n_cKO"]))==(6,6),"Primary association does not use all twelve mice")
        require(not tested or num(r["effect_cKO_minus_WT"],True) is not None and probability(r["pvalue"]) is not None and probability(r["association_fdr"]) is not None,"Tested association missing native effect/P/FDR")
        effect={**r,"method_id":method(lin),"endpoint_id":ENDPOINTS[ep][0],"feature_id":feature(lin,reg),"contrast":"cKO_vs_WT",
          "effect_name":"Mean native AUCell AUC difference (cKO − WT)","effect":num(r["effect_cKO_minus_WT"]),"statistic":num(r.get("statistic")),
          "pvalue":probability(r["pvalue"]),"padj":probability(r["association_fdr"]),"association_fdr":probability(r["association_fdr"]),"tested":tested,
          "n_WT":integer(r["n_WT"]),"n_cKO":integer(r["n_cKO"]),"source_fit_id":"native_regulon_mouse_lm_primary","source_table":source_table}
        for k in ("mean_WT","mean_cKO","effect_cKO_minus_WT","SE","CI_low","CI_high"):
            if k in r:effect[k]=num(r[k])
        effects.append(effect)
    require(seen==set(cov),"Native association keys incomplete")
    return effects

def validate_nucleus_axes(native,cfg,regulons):
    """Inspect original nucleus/regulon axes only; never recompute scores or means."""
    audits={}
    for lin in ("CA1","IN"):
        cells=table(cfg["lineages"][lin]["cells"]["path"],["nucleus_id","endpoint","sample","condition"])
        source={r["nucleus_id"]:r for r in cells};require(len(source)==len(cells),"Duplicate frozen lineage nucleus IDs")
        auc=native/lin/"aucell/nucleus_AUC.tsv.gz"
        nmeta=native/lin/"aucell/nucleus_axis_metadata.tsv.gz"
        require(nmeta.is_file(),"Missing exact nucleus AUCell metadata")
        meta=table(nmeta,["nucleus_id","endpoint","sample","condition"]); axis=[r["nucleus_id"] for r in meta]
        require(len(set(axis))==len(axis) and set(axis)==set(source),"Native nucleus metadata differs from exact frozen lineage")
        for r in meta:
            c=source[r["nucleus_id"]]
            require(all(r[k]==c[k] for k in ("endpoint","sample","condition")),"Native nucleus endpoint/mouse identity changed")
        with gzip.open(auc,"rt",encoding="utf-8",newline="") as h:
            reader=csv.reader(h,delimiter="\t");header=next(reader)
            require(header[0]=="nucleus_id" and len(header)==len(set(header)),"Invalid nucleus AUCell header")
            expected=set(regulons[lin]);require(set(header[1:])==expected,"Native AUC regulon axis differs from scored regulons")
            seen=[]
            for row in reader:
                require(len(row)==len(header),"Malformed native nucleus AUCell matrix")
                seen.append(row[0])
                for v in row[1:]:
                    require(probability(v) is not None,"Missing/invalid native nucleus AUC")
            require(seen==axis,"Native AUC nucleus order differs from matching metadata")
        audits[lin]={"n_nuclei":len(axis),"n_scored_regulons":len(expected),"exact_source_nucleus_axis_verified":True,
          "auc_table":str(auc.relative_to(native)),"metadata_table":str(nmeta.relative_to(native)),"no_score_recomputation_or_averaging":True}
    return audits

def package(config_path,results):
    config_path,results=Path(config_path).resolve(),Path(results).resolve();cfg_record,cfg=receipt(config_path),load(config_path)
    require(results==Path(cfg["run_root"]).resolve()/"results","Results differ from immutable run_root")
    require(not (results/"FINAL_VALIDATION.json").exists(),"Refusing to replace completed motif transport")
    out=results/"atlas";require(not out.is_symlink() and (not out.exists() or out.is_dir() and not any(out.iterdir())),"Refusing nonempty atlas output")
    inputs_before=verify_inputs(cfg);src=Path(__file__).resolve().parent
    executed={"scripts/"+n:digest(src/n) for n in ("run_pyscenic_native.py","analyze_regulon_mouse.R","package_pyscenic_results.py")}
    native,mouse=results/"pyscenic_native",results/"regulon_mouse"
    native_marker=sealed_step(native,"satb2-pyscenic-native-complete/v1",REQUIRED_NATIVE,cfg_record,src/"run_pyscenic_native.py")
    mouse_marker=sealed_step(mouse,"satb2-regulon-mouse-complete/v1",REQUIRED_MOUSE,cfg_record,src/"analyze_regulon_mouse.R")
    for k in ("positive_motif_regulons_only","full_genotype_neutral_RNA_background","no_DEG_selection","native_nucleus_AUCell","exact_nucleus_axes_verified","all12_mouse_roster","min_measured_targets_10","Satb2_self_target_sensitivity_separate","inputs_config_source_SHA256_unchanged"):
        require(native_marker.get("checks",{}).get(k) is True,f"Native motif check failed: {k}")
    for k in ("all12_mouse_roster","native_AUCell_means_preserved","global_primary_BH_both_lineages","separate_sensitivity_BH","no_activity_threshold_filter","no_nucleus_pseudoreplication","no_new_gene_DESeq2","positive_motif_regulons_only","inputs_config_source_SHA256_unchanged"):
        require(mouse_marker.get("checks",{}).get(k) is True,f"Mouse regulon check failed: {k}")
    frozen,samples,source_coverage=frozen_cohort(cfg);provenance=load(native/"PROVENANCE.json")
    require(provenance["config"]["sha256"]==cfg_record["sha256"] and provenance["source"]["sha256"]==executed["scripts/run_pyscenic_native.py"],"Native provenance differs from executed source/config")
    statuses=table(native/"LINEAGE_STATUS.tsv",["lineage","status","n_regulons","n_nuclei","n_genes"]);status={r["lineage"]:r for r in statuses}
    require(len(status)==len(statuses)==2 and set(status)=={"CA1","IN"},"Wrong native lineage status map")
    columns=["lineage","endpoint","endpoint_id","source_endpoint_id","regulon_id","native_name","TF","n_targets","n_background_genes","eligible","reason"]
    cov_rows=table(native/"REGULON_COVERAGE.tsv",columns);cov={};regs={};bylin=defaultdict(set)
    for r in cov_rows:
        ep,lin,reg=identity(r);require((ep,reg) not in cov,"Duplicate native regulon coverage");cov[ep,reg]=r
        require(integer(r["n_targets"])>0 and integer(r["n_background_genes"])>=integer(r["n_targets"]),"Invalid native regulon target/background count")
        require(boolean(r["eligible"])==(integer(r["n_targets"])>=10),"Native measured-target score eligibility differs from fixed >=10 policy")
        if reg in regs:require(all(regs[reg][k]==r[k] for k in ("lineage","native_name","TF","n_targets","n_background_genes","eligible")),"Native lineage regulon definition changes across endpoints")
        else:regs[reg]=r
        bylin[lin].add(reg)
    for lin in ("CA1","IN"):
        s=status[lin];require(s["status"] in ("COMPLETE_WITH_REGULONS","EMPTY_NO_SUPPORTED_REGULONS"),"Native lineage incomplete")
        require(integer(s["n_regulons"])==len(bylin[lin]),"Native regulon count changed")
        require((s["status"]=="EMPTY_NO_SUPPORTED_REGULONS")== (len(bylin[lin])==0),"Native empty-lineage state disagrees")
        require(integer(s["n_nuclei"])==sum(info[2] for ep,info in ENDPOINTS.items() if lineage(ep)==lin),"Native lineage nucleus count differs")
        require({(ep,reg) for ep in ENDPOINTS if lineage(ep)==lin for reg in bylin[lin]}=={key for key in cov if lineage(key[0])==lin},"Each lineage regulon needs each matching frozen endpoint")
    values=validate_score_maps(table(native/"REGULON_ACTIVITY_PER_MOUSE.tsv"),cov,samples,source_coverage)
    r_scores=table(mouse/"REGULON_MOUSE_SCORES.tsv");rmap={(r["endpoint"],r["regulon_id"],r["sample"]):r for r in r_scores}
    require(len(rmap)==len(r_scores) and set(rmap)==set(values),"Mouse association score axes differ from native AUCell")
    for key,r in rmap.items():
        a=values[key];require(same(num(a["aucell_mean"],True),num(r["aucell_mean"],True)) and a["condition"]==r["condition"] and a["n_nuclei"]==r["n_nuclei"],"R changed original native AUCell means/membership")
    scored_bylineage={lin:{reg for reg in bylin[lin] if boolean(regs[reg]["eligible"])} for lin in ("CA1","IN")}
    axis_audit=validate_nucleus_axes(native,cfg,scored_bylineage)
    targets=table(native/"REGULON_TARGETS.tsv",["lineage","regulon_id","TF","target","weight"]);edge_map=defaultdict(list);seen=set()
    for r in targets:
        reg=r["regulon_id"];require(reg in regs and r["lineage"]==regs[reg]["lineage"] and r["TF"]==regs[reg]["TF"],"Unknown native regulon target")
        require((reg,r["target"]) not in seen and r["target"],"Duplicate/empty native regulon target");seen.add((reg,r["target"]))
        require(num(r["weight"],True)>=0,"Invalid native gene2weight")
        if "self_target" in r:require(boolean(r["self_target"])==(r["TF"]==r["target"]),"Native TF-self target flag differs")
        edge_map[reg].append(r)
    require(all(len(edge_map[reg])==integer(r["n_targets"]) for reg,r in regs.items()),"Native target membership count differs from coverage")
    # The sensitivity retains its own graph/score definition and multiplicity family.
    no_cov_rows=table(native/"SATB2_NO_SELF_TARGET_COVERAGE.tsv",columns+["network_variant"]);no_cov={}
    for r in no_cov_rows:
        ep,lin,reg=identity(r);require(r["TF"]=="Satb2" and reg==lin+"::Satb2(+)_NO_SELF_TARGET" and (ep,reg) not in no_cov,"Wrong no-self identity")
        no_cov[ep,reg]=r
    require(len(no_cov)==5 and {ep for ep,_ in no_cov}==set(ENDPOINTS),"Need explicit all-five-endpoint Satb2 no-self coverage")
    no_values=validate_score_maps(table(native/"SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv"),no_cov,samples,source_coverage)
    no_r=table(mouse/"SATB2_NO_SELF_TARGET_MOUSE_SCORES.tsv");no_rmap={(r["endpoint"],r["regulon_id"],r["sample"]):r for r in no_r}
    require(len(no_rmap)==len(no_r) and set(no_rmap)==set(no_values),"No-self native/R score axes differ")
    for key,r in no_rmap.items():require(same(num(r["aucell_mean"],True),num(no_values[key]["aucell_mean"],True)),"R changed no-self AUCell mean")
    for (ep,reg),r in no_cov.items():
        native_reg=lineage(ep)+"::Satb2(+)";main=cov.get((ep,native_reg))
        if boolean(r["eligible"]):
            require(main is not None and integer(r["n_targets"])>=10 and integer(r["n_targets"])==integer(main["n_targets"])-1,"No-self target count changed beyond own gene removal")
            require(integer(r["n_background_genes"])==integer(main["n_background_genes"]),"No-self RNA background changed")
            require(any(x["target"]=="Satb2" for x in edge_map[native_reg]),"No-self sensitivity without native selfgene")
    no_effects=associations(mouse/"SATB2_NO_SELF_TARGET_ASSOCIATIONS.tsv",no_cov,no_values,samples,"regulon_mouse/SATB2_NO_SELF_TARGET_ASSOCIATIONS.tsv")
    payload={"schema":SCHEMA,"analysis_id":cfg["analysis_id"],"source_fingerprint":"","created_at":datetime.now(timezone.utc).isoformat(),"status":"complete",
      "methods":[],"endpoints":[],"samples":[{"id":s,"genotype":samples[s]} for s in sorted(samples,key=lambda s:(samples[s]!="WT",s))],
      "features":[],"values":[],"effects":[],"hubs":[],"nodes":[],"edges":[],"embeddings":[],"coverage":[],"downloads":[]}
    for lin in ("CA1","IN"):
        eligible_count=len(scored_bylineage[lin]);label="CA1" if lin=="CA1" else "Inhibitory"
        payload["methods"].append({"id":method(lin),"label":label+" positive motif regulons · pySCENIC","family":"motif_regulon","status":"complete" if eligible_count else "not_available",
          "score_label":"Native AUCell AUC mean across each mouse's exact nuclei","units":"AUCell AUC","activity_unit":"mouse","endpoint_ids":[v[0] for ep,v in ENDPOINTS.items() if lineage(ep)==lin],
          "network_edge_basis":"Native motif-pruned regulon gene2weight (GRNBoost2 importance), including pySCENIC-added TF selfgene; candidate orientation is not verified binding or binding strength",
          "interpretation":"Positive expression-derived modules with mouse mm10 cisTarget motif enrichment; native AUCell scored per nucleus, averaged only within exact endpoint and independent mouse. AUC is target-set enrichment, not TF RNA, direct binding, or causal activity. Genotype partly confounds technical partition. Definitions differ by lineage and must not be merged by TF name. Primary BH is joint across every tested matching endpoint × regulon, both lineages; no-self Satb2 is a separate sensitivity.",
          "eligibility_policy":"Full genotype-neutral RNA background; no DEG selection. At least10 measured motif-supported native targets including native TF selfgene; all12 independent mice for primary models. Unavailable activity stays null.",
          "resources":{"native_lineage_status":status[lin],"native_nucleus_axis_audit":axis_audit[lin],"references":cfg["references"],"native_parameters":provenance}})
    for ep,(ui,source_id,count,label) in ENDPOINTS.items():
        payload["endpoints"].append({"id":ui,"label":label,"family":lineage(ep),"source_endpoint_id":source_id,"membership_sha256":frozen[ep]["membership_sha256"],"n_nuclei":count,
          "identity_warning":"Molecular expression program; physical layer undetermined" if lineage(ep)=="CA1" else "Mixed atlas identity remains under review" if ep=="IN_C8" else "Molecular reference resemblance, not connectivity/morphology"})
    for reg,r in sorted(regs.items()):
        payload["features"].append({"id":feature(r["lineage"],reg),"method_id":method(r["lineage"]),"label":r["native_name"],"kind":"regulon","source_id":reg,"TF":r["TF"],"native_name":r["native_name"],"units":"Native AUCell AUC","n_targets":integer(r["n_targets"]),
          "self_target_caveat":"Native pySCENIC target membership may include own TF; retained as native membership, not self-binding proof","color":"#9c65a9"})
    for (ep,reg),r in sorted(cov.items()):
        lin=lineage(ep);eligible=boolean(r["eligible"]);nt=integer(r["n_targets"])
        payload["coverage"].append({**r,"method_id":method(lin),"endpoint_id":ENDPOINTS[ep][0],"feature_id":feature(lin,reg),"eligible":eligible,"n_targets":nt,"n_background_genes":integer(r["n_background_genes"]),"measured_targets":nt,"required_targets":10})
        for sample in samples:
            a=values.get((ep,reg,sample));row={"method_id":method(lin),"endpoint_id":ENDPOINTS[ep][0],"feature_id":feature(lin,reg),"sample_id":sample,"value":None,"eligible":eligible,
              "n_nuclei":integer(source_coverage[ep,sample]["n_nuclei"]),"coverage_n":nt,"coverage_total":nt,"source_fit_id":"native_nucleus_AUCell_endpoint_mouse_mean","missing_reason":"" if eligible else r["reason"],"normalization":"Native AUCell on complete genotype-neutral normalized RNA ranking; unchanged native AUC"}
            if a:row.update({"value":num(a["aucell_mean"],True),"aucell_mean":num(a["aucell_mean"],True),"genotype_association_tested":False})
            payload["values"].append(row)
        # Graphs show native full motif-pruned memberships, including low-target
        # unmodelable sets; no unmeasured edge score is imputed.
        for edge in edge_map[reg]:
            tf,target=edge["TF"],edge["target"]
            for symbol,kind in ((tf,"tf"),(target,"target")):
                ident=node_id(ep,reg,symbol);payload["nodes"].append({"id":ident,"method_id":method(lin),"endpoint_id":ENDPOINTS[ep][0],"feature_id":feature(lin,reg),"label":symbol,"source_id":symbol,"kind":"tf" if symbol==tf else kind})
            payload["edges"].append({"method_id":method(lin),"endpoint_id":ENDPOINTS[ep][0],"feature_id":feature(lin,reg),"source":node_id(ep,reg,tf),"target":node_id(ep,reg,target),"weight":num(edge["weight"],True),"sign":"positive","directed":True,"self_target":tf==target,
              "edge_basis":"native_positive_motif_pruned_regulon_gene2weight_not_binding_strength_not_verified_causation"})
    payload["nodes"]=list({x["id"]:x for x in payload["nodes"]}.values())
    payload["effects"]=associations(mouse/"REGULON_ASSOCIATIONS.tsv",cov,values,samples,"regulon_mouse/REGULON_ASSOCIATIONS.tsv")
    tf_symbols=set(Path(cfg["references"]["tf_names"]["path"]).read_text().split())
    satb2_candidates={}
    for lin in ("CA1","IN"):
        candidates=table(native/lin/"prepared/TF_CANDIDATE_COVERAGE.tsv",["lineage","TF","official_TF_candidate","assayed_symbol","eligible_RNA","reason","forced_into_TF_list"])
        require(all(r["lineage"]==lin and not boolean(r["forced_into_TF_list"]) for r in candidates),"A TF was forced into the native candidate list")
        sat=[r for r in candidates if r["TF"]=="Satb2"]
        require(len(sat)==1 and boolean(sat[0]["official_TF_candidate"])==("Satb2" in tf_symbols),"Native Satb2 official TF-list coverage disagrees")
        satb2_candidates[lin]=sat[0]
    payload["Satb2_status"]={"official_TF_list_membership": "Satb2" in tf_symbols,"TF_list_reference":cfg["references"]["tf_names"],"native_TF_candidate_coverage":satb2_candidates,
      "native_endpoint_status":table(mouse/"SATB2_STATUS.tsv"),"interpretation":"Official TF-list membership does not guarantee measured eligibility or a motif-supported regulon; unsupported activity is not substituted by Satb2 RNA or forced into inference."}
    native_files=[]
    for folder in (native,mouse):
        for path in sorted(folder.rglob("*")):
            require(not path.is_symlink(),"Symlink in native output tree")
            if path.is_file():native_files.append(receipt(path,results))
    for record in native_files:
        if record["path"].endswith((".tsv",".tsv.gz",".json",".csv",".csv.gz")):payload["downloads"].append({**record,"native_bytes_preserved":True,"label":record["path"]})
    payload["native_supplements"]={"Satb2_no_self_target":{"coverage":no_cov_rows,"effects":no_effects,"tables":[r for r in native_files if "SATB2_NO_SELF_TARGET" in r["path"]],"primary_full_regulons_unchanged":True,
      "interpretation":"Only native TF selfgene removed when genuine Satb2 motif regulon supports at least10 remaining targets; same RNA rankings/background. Own gene addition is pySCENIC bookkeeping, not direct self-binding evidence. Sensitivity has independent five-endpoint BH family, not independent validation."},
      "mouse_sensitivities":[r for r in native_files if any(s in r["path"] for s in ("SAME_PARTITION","LOO_","BOOTSTRAP"))]}
    inputs_after=verify_inputs(cfg);require(inputs_before==inputs_after and receipt(config_path)==cfg_record,"Source input/config changed while packaging")
    for r in native_files:require(receipt(child(results,r["path"]),results)==r,"Native output changed while packaging")
    require(executed=={"scripts/"+n:digest(src/n) for n in ("run_pyscenic_native.py","analyze_regulon_mouse.R","package_pyscenic_results.py")},"Executed source changed")
    fingerprint={"configuration_sha256":cfg_record["sha256"],"inputs":inputs_before,"native_outputs":native_files,"memberships":{ep:frozen[ep]["membership_sha256"] for ep in ENDPOINTS},"executed_source_sha256":executed}
    payload["source_fingerprint"]=hashlib.sha256(raw_json(fingerprint,True)).hexdigest()
    payload["provenance"]={"configuration_sha256":cfg_record["sha256"],"configuration":cfg,"bound_inputs":inputs_before,"executed_source_sha256":executed,"packager_sha256":digest(__file__),"native_output_receipts":native_files,
      "native_motif_provenance":provenance,"native_mouse_association_complete":mouse_marker,"native_nucleus_axis_audit":axis_audit,"no_scientific_fits_in_packager":True,
      "TF_activity":"NATIVE_NUCLEUS_AUCELL_AND_INDEPENDENT_MOUSE_ASSOCIATIONS_COMPLETE","motif_supported_regulons":"NATIVE_POSITIVE_PYSCENIC_COMPLETE","activity_UMAP":"NOT_PACKAGED_NO_BOUND_COORDINATE_AXIS","original_DESeq2_and_LR_results_unchanged":True}
    raw=raw_json(payload,True);compressed=gzip.compress(raw,compresslevel=9,mtime=0)
    outputs={"data.json":raw,"data.json.gz":compressed,"data.js":b"/* Completed native motif-supported regulon AUCell; no fits in transport. */\nwindow.SATB2_REGULATORY_DATA_GZIP="+json.dumps(base64.b64encode(compressed).decode()).encode()+b";\n","downloads.json":raw_json(payload["downloads"])}
    transport={n:{"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()} for n,b in outputs.items()}
    manifest={"schema":"satb2-regulatory-transport/v1","analysis_id":cfg["analysis_id"],"source_fingerprint":payload["source_fingerprint"],"status":"complete","completed_primary":"positive_motif_supported_regulons_native_nucleus_AUCell_and_independent_mouse_associations",
      "record_counts":{k:len(payload[k]) for k in ("methods","endpoints","samples","features","values","effects","nodes","edges","coverage","embeddings")},"transport":transport,"scientific_fits_in_packager":False,"original_sources_unchanged":True}
    out.mkdir(parents=True,exist_ok=True)
    for name,content in {**outputs,"manifest.json":raw_json(manifest)}.items():
        with (out/name).open("xb") as handle:handle.write(content)
    final={"schema":"satb2-pyscenic-regulatory-final-validation/v1","status":"COMPLETE","analysis_id":cfg["analysis_id"],"configuration_sha256":cfg_record["sha256"],"source_fingerprint":payload["source_fingerprint"],"executed_source_sha256":executed,
      "TF_activity_run":True,"motif_regulon_run":True,"positive_motif_regulons_only":True,"activity_UMAP_packaged":False,"twelve_mice_per_endpoint_regulon":True,"unique_mouse_keys":True,"raw_mouse_values_finite":True,
      "frozen_membership_hashes_verified":True,"native_score_receipts_verified":True,"exact_nucleus_axes_verified":True,"native_nucleus_axis_audit":axis_audit,"Satb2_self_target_sensitivity_separate":True,"Satb2_official_TF_list_membership":payload["Satb2_status"]["official_TF_list_membership"],"minimum_measured_targets":10,
      "all_configured_inputs_SHA256_before_after_verified":True,"original_analysis_files_unchanged":True,"no_scientific_fits_in_packager":True,"frozen_endpoint_counts":{ep:i[2] for ep,i in ENDPOINTS.items()},
      "lineage_status":status,"native_file_receipts":native_files,"package_file_receipts":[receipt(p,results) for p in sorted(out.iterdir()) if p.is_file()],"atlas":manifest}
    with (results/"FINAL_VALIDATION.json").open("xb") as handle:handle.write(raw_json(final))
    return final

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--config",required=True,type=Path);parser.add_argument("--results",required=True,type=Path);args=parser.parse_args()
    require(args.config.is_absolute() and args.results.is_absolute(),"Use absolute immutable paths")
    final=package(args.config,args.results);print(json.dumps({"status":final["status"],"analysis_id":final["analysis_id"],"counts":final["atlas"]["record_counts"]},indent=2))
if __name__=="__main__":main()
