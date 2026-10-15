#!/usr/bin/env python3
"""Package produced annotation/UMAP/marker tables; no RNA extraction or model fitting."""
import argparse, csv, gzip, hashlib, json, math, re, shutil, time
from collections import defaultdict
from pathlib import Path

SCHEMA = 'satb2-offline-1'
TOOLS = ['scanpy', 'seurat']
BRANCHES = ['original_full', 'paired_original', 'paired_cellbender']
POPS = ['allcells', 'Excitatory', 'Inhibitory', 'CA1_ProS']
KEY = ['tool', 'branch', 'population', 'cluster_scope', 'cluster_id']
N_FIELDS = ['nucleus_id','sample','condition','cluster_id','source_cluster_id','current_cluster_id','reviewed_label','broad_label','identity_label','reference_class_name','reference_subclass_name','reference_supertype_name','reference_cluster_name','reference_broad_prediction','reference_class_bootstrapping_probability','reference_subclass_bootstrapping_probability','reference_supertype_bootstrapping_probability','reference_cluster_bootstrapping_probability','reference_broad_high_support','reference_subclass_high_support','reference_identity_review_status','layer_program_label','expression_program_label','label_confidence','physical_layer_assignment','identity_review_status','uncertainty_reason','keep_for_analysis','keep_for_plot','exclusion_reason','program_evidence_status','source_count_branch']
NUMERIC = {'reference_class_bootstrapping_probability','reference_subclass_bootstrapping_probability','reference_supertype_bootstrapping_probability','reference_cluster_bootstrapping_probability','n_nuclei','original_cluster_n_nuclei','source_cluster_n_nuclei_before_filter','atlas_majority_fraction','atlas_broad_majority_fraction','atlas_broad_high_support_fraction','atlas_subclass_high_support_fraction','n_mice_represented','n_WT_program_evidence_nuclei','n_WT_program_evidence_mice','superficial_program_WT_expression_z','deep_program_WT_expression_z'}
BOOL = {'reference_broad_high_support','reference_subclass_high_support','keep_for_analysis','keep_for_plot','assignment_carried_from_original','atomic_cluster_ID_changed'}
QC = [('qc_original_nCount_RNA','Original full-library UMI count'),('qc_original_nFeature_RNA','Original detected genes'),('branch_source_nCount_RNA','Chosen branch full-library UMI count'),('qc_pct_mt','Original mitochondrial percentage'),('qc_pct_ribo','Original ribosomal percentage'),('qc_doublet_score','Original scDblFinder doublet score'),('qc_doublet_call','Original doublet call'),('qc_strict_sensitivity_pass','Original strict-QC sensitivity flag'),('qc_high_genes_review','High detected-genes review flag'),('qc_umi_fraction_removed','CellBender UMI fraction removed'),('qc_corrected_nFeature_RNA','CellBender corrected detected genes')]


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for z in iter(lambda:f.read(1024*1024),b''): h.update(z)
    return h.hexdigest()

def number(v):
    if v is None or v in ['', 'NA', 'NaN', 'nan', 'Inf', '-Inf']: return None
    x=float(v)
    return x if math.isfinite(x) else None

def boolean(v):
    if str(v).lower() not in ['true','false','1','0']: raise ValueError('Invalid boolean '+str(v))
    return str(v).lower() in ['true','1']

def typed(row):
    return {k:(boolean(v) if k in BOOL else number(v) if k in NUMERIC else None if v in ['NA',''] else v) for k,v in row.items()}

def ident(tool,branch,pop): return '__'.join([tool,branch,pop])

def rows(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else Path(p).open()) as f:
        yield from csv.DictReader(f,delimiter='\t')

class Builder:
    def __init__(self, cfg, out):
        self.cfg_path=Path(cfg).resolve(); self.cfg=json.loads(self.cfg_path.read_text()); self.root=Path(self.cfg['results_root']).resolve(); self.old=Path(self.cfg['previous_run_root']).resolve(); self.out=Path(out).resolve(); self.sources={}; self.shards={}; self.expected={}; self.assets=[]
        self.out.mkdir(parents=True,exist_ok=True); (self.out/'data/shards').mkdir(parents=True,exist_ok=True); (self.out/'staging').mkdir(exist_ok=True)
    def bind(self,p):
        p=Path(p).resolve(strict=True)
        if not (p.is_relative_to(self.root) or p.is_relative_to(self.old) or p.is_relative_to(Path(self.cfg['study_root']).resolve())): raise ValueError('Unapproved source path')
        s=sha(p)
        if str(p) in self.sources and self.sources[str(p)] != s: raise ValueError('Changed source '+str(p))
        
        if str(p) in self.expected and s != self.expected[str(p)]: raise ValueError('Producer artifact changed '+str(p))
        self.sources[str(p)]=s; return p
    def complete(self,name,status):
        p=self.bind(self.root/name/'COMPLETE.json'); d=json.loads(p.read_text())
        if d['status'] != status: raise ValueError('Incomplete producer '+name)
        if d.get('config_sha256') != sha(self.cfg_path): raise ValueError('Config binding differs '+name)
        for relative,h in d.get('output_sha256',{}).items():
            target=(self.root/name/relative).resolve()
            if not target.is_relative_to(self.root/name):raise ValueError('Producer receipt escapes root')
            self.expected[str(target)]=h
        for artifact in d.get('artifacts',[]):
            if isinstance(artifact,dict) and 'path' in artifact and 'sha256' in artifact:self.expected[str(Path(artifact['path']).resolve())]=artifact['sha256']
        return d
    def shard(self,id,payload):
        import base64
        raw=json.dumps(payload,ensure_ascii=False,separators=(',',':'),allow_nan=False).encode()
        packed=gzip.compress(raw,compresslevel=6,mtime=0); encoded=base64.b64encode(packed).decode()
        path='data/shards/'+id+'.js'; p=self.out/path
        p.write_text('window.SATB2_SHARDS=window.SATB2_SHARDS||{};window.SATB2_SHARDS['+json.dumps(id)+']={encoding:"gzip-base64",data:'+json.dumps(encoded)+'};\n')
        self.shards[id]={'path':path,'sha256':sha(p),'uncompressed_bytes':len(raw),'gzip_bytes':len(packed)}
    def copy_asset(self,p,relative):
        p=self.bind(p); out=self.out/relative;out.parent.mkdir(parents=True,exist_ok=True)
        if not out.exists() or sha(out)!=self.sources[str(p)]: shutil.copy2(p,out)
        entry={'source_path':str(p),'bundle_path':relative,'path':relative,'sha256':sha(out),'bytes':out.stat().st_size};self.assets.append(entry);return entry
    def metadata(self,tool,branch,pop):
        if tool=='scanpy': return self.old/'results/scanpy'/branch/pop/'nucleus_metadata.tsv'
        sub='log_global' if pop=='allcells' else 'annotate_subsets/'+{'Excitatory':'EX','Inhibitory':'IN','CA1_ProS':'CA1_PROS_CANDIDATE'}[pop]
        return self.old/'results/seurat'/branch/sub/'nucleus_metadata.tsv'
    def build(self):
        self.complete('annotation_final','ANNOTATION_FINAL_COMPLETE_REVIEW_REQUIRED'); self.complete('marker_summaries','COMPLETE'); self.complete('plots','COMPLETE')
        cv=self.complete('count_validation','COMPLETE')
        if cv.get('exact_integer_counts_PASS') is not True: raise ValueError('Count acceptance absent')
        genes=list(rows(self.bind(self.root/'annotation_initial/gene_panel.tsv')))
        if len({g['gene_id'] for g in genes})!=len(genes): raise ValueError('Gene IDs nonunique')
        cluster_groups=defaultdict(list)
        for z in rows(self.bind(self.root/'annotation_final/cluster_annotations.tsv')):
            if z['cluster_scope']=='original_frozen': cluster_groups[ident(z['tool'],z['branch'],z['population'])].append(typed(z))
        expected={ident(t,b,p) for t in TOOLS for b in BRANCHES for p in POPS}
        if set(cluster_groups)!=expected: raise ValueError('24 original frozen datasets absent')
        annodir=self.out/'staging/annotations';annodir.mkdir(exist_ok=True)
        handles={d:(annodir/(d+'.tsv')).open('w') for d in expected}
        writers={d:csv.DictWriter(h,fieldnames=N_FIELDS,delimiter='\t') for d,h in handles.items()}
        for w in writers.values():w.writeheader()
        for z in rows(self.bind(self.root/'annotation_final/nucleus_annotations.tsv.gz')):
            if z['cluster_scope']=='original_frozen' and boolean(z['keep_for_analysis']): writers[ident(z['tool'],z['branch'],z['population'])].writerow({k:z[k] for k in N_FIELDS})
        for h in handles.values():h.close()
        index=self.out/'staging/DATASET_COORDINATE_INDEX.tsv.gz'
        datasets=[];sample_map={};qc_cache={};books={}
        for z in rows(self.bind(self.root/'plots/PLOT_BOOK_INDEX.tsv')):
            if z['cluster_scope']=='original_frozen':books[ident(z['tool'],z['branch'],z['population'])]=z
        with gzip.open(index,'wt') as ix:
            iw=csv.DictWriter(ix,fieldnames=['dataset_id','tool','branch','population','nucleus_id','row_index'],delimiter='\t');iw.writeheader()
            for tool in TOOLS:
              for branch in BRANCHES:
                if branch not in qc_cache:
                    qkeys=['sample','condition','frozen_original_nCount_RNA','nCount_RNA','nFeature_RNA','branch_UMIs_verified','pct_mt','pct_ribo','score_new','call_new','passes_strict_low_sensitivity','flag_high_genes_review','umi_fraction_removed','cellbender_nFeature_RNA','sequencing_partition']
                    qc_cache[branch]={z['nucleus_id']:{k:z.get(k) for k in qkeys} for z in rows(self.bind(self.old/'results/shared'/branch/'nucleus_metadata.tsv'))}
                qc=qc_cache[branch]
                for pop in POPS:
                    d=ident(tool,branch,pop); ann=sorted(list(rows(annodir/(d+'.tsv'))),key=lambda x:x['nucleus_id'])
                    if len({z['nucleus_id'] for z in ann})!=len(ann):raise ValueError('Duplicate nucleus ID')
                    mkeys=['sample','condition']+(['UMAP_harmony_1','UMAP_harmony_2','UMAP_unintegrated_1','UMAP_unintegrated_2'] if tool=='scanpy' else [])
                    meta={z['nucleus_id']:{k:z[k] for k in mkeys} for z in rows(self.bind(self.metadata(tool,branch,pop)))}
                    coord={}
                    if tool=='seurat':
                        folder=self.metadata(tool,branch,pop).parent
                        for basis,prefix in [('harmony','UMAPH'),('unintegrated','UMAPU')]:
                            for z in rows(self.bind(folder/('umap_'+basis+'.tsv'))):coord.setdefault(z['nucleus_id'],{}).update({'x_'+basis:number(z[prefix+'_1']),'y_'+basis:number(z[prefix+'_2'])})
                    nuclei=[]
                    for i,z in enumerate(ann):
                        nid=z['nucleus_id'];m=meta[nid];q=qc[nid]
                        if (z['sample'],z['condition'])!=(q['sample'],q['condition']) or (z['sample'],z['condition'])!=(m['sample'],m['condition']):raise ValueError('Nucleus animal/genotype mismatch')
                        r=typed(z)
                        if tool=='scanpy':r.update({k:number(m[v]) for k,v in [('x_harmony','UMAP_harmony_1'),('y_harmony','UMAP_harmony_2'),('x_unintegrated','UMAP_unintegrated_1'),('y_unintegrated','UMAP_unintegrated_2')]})
                        else:r.update(coord[nid])
                        if any(r[k] is None for k in ['x_harmony','y_harmony','x_unintegrated','y_unintegrated']):raise ValueError('Missing UMAP coordinate')
                        r.update({'qc_original_nCount_RNA':number(q.get('frozen_original_nCount_RNA',q['nCount_RNA'])),'qc_original_nFeature_RNA':number(q['nFeature_RNA']),'branch_source_nCount_RNA':number(q.get('branch_UMIs_verified')),'qc_pct_mt':number(q['pct_mt']),'qc_pct_ribo':number(q['pct_ribo']),'qc_doublet_score':number(q['score_new']),'qc_doublet_call':q['call_new'],'qc_strict_sensitivity_pass':boolean(q['passes_strict_low_sensitivity']),'qc_high_genes_review':boolean(q['flag_high_genes_review']),'qc_umi_fraction_removed':number(q.get('umi_fraction_removed')),'qc_corrected_nFeature_RNA':number(q.get('cellbender_nFeature_RNA')),'sequencing_partition':q['sequencing_partition']})
                        nuclei.append(r);sample_map[z['sample']]=z['condition'];iw.writerow({'dataset_id':d,'tool':tool,'branch':branch,'population':pop,'nucleus_id':nid,'row_index':i})
                    count=sum(int(c['n_nuclei']) for c in cluster_groups[d])
                    if count!=len(nuclei):raise ValueError('Annotation cluster denominator mismatch')
                    self.shard('coordinates__'+d,{'nuclei':nuclei,'clusters':cluster_groups[d],'coordinate_measurement':'saved_native_UMAP; no refit after shared mask','reference_scope':'source-verified reused original MapMyCells; shared reference is not independent truth'})
                    book=books[d];asset=self.copy_asset(book['pdf'],'assets/pdf/main/'+d+'.pdf')
                    datasets.append({'id':d,'tool':tool,'branch':branch,'population':pop,'cluster_scope':'original_frozen','coordinates_shard':'coordinates__'+d,'markers_shard':'markers__'+d,'fits_shard':'fits__'+d,'n_nuclei':len(nuclei),'coordinate_fit_nuclei':int(book['coordinate_fit_nuclei']),'coordinate_refit_performed':False,'pdf_path':asset['path'],'baseline_definition':'final keep_for_analysis shared iterative mask; original Harmony r0.4 clusters'})
                    print(json.dumps({'dataset':d,'n_nuclei':len(nuclei),'coordinates_gzip_MB':round(self.shards['coordinates__'+d]['gzip_bytes']/1e6,2)}),flush=True)
                del qc_cache[branch]
        # Sum-based cKO arithmetic uses already produced pooled/WT aggregates only.
        pooled=rows(self.bind(self.root/'marker_summaries/gene_cluster_summary.tsv'));wt=rows(self.bind(self.root/'marker_summaries/gene_cluster_summary_WT.tsv'))
        completed=set();current_d=None;payload={'rows':[],'WT_rows':[],'cKO_rows':[]}
        def save_markers(d,payload):
            payload.update({'gene_overlay_measurement':'existing cluster mean log1p full-library-normalized expression; not per-nucleus RNA','normalization':'log1p(chosen-source count / complete chosen-source nucleus library * 10000), existing native sums','cKO_derivation':'difference of existing pooled and WT numerators/denominators; cKO z is unavailable and remains null'})
            self.shard('markers__'+d,payload);completed.add(d)
            print(json.dumps({'markers':d,'n_gene_cluster_rows':len(payload['rows'])}),flush=True)
        from itertools import zip_longest
        fields=['cluster_id','gene_id','gene_symbol','n_nuclei','sum_counts','n_detected','sum_log1p_normalized','mean_log1p_normalized','fraction_detected','mean_source_counts','scaled_expression_z','n_mice_with_nuclei','n_mice_eligible_for_expression_review','n_mice_expression_supported','pct_difference_vs_rest','mean_log1p_difference_vs_rest','keep_for_plot','measurement_status','scale_reference_branch','scale_reference_mean','scale_reference_sd','scale_definition']
        for a,w in zip_longest(pooled,wt):
            if a is None or w is None or any(a[k]!=w[k] for k in KEY+['gene_id']):raise ValueError('Pooled/WT aggregate row binding mismatch')
            d=ident(a['tool'],a['branch'],a['population'])
            if current_d is not None and current_d != d:
                save_markers(current_d,payload);payload={'rows':[],'WT_rows':[],'cKO_rows':[]}
            if d in completed:raise ValueError('Noncontiguous marker dataset rows')
            current_d=d
            def keep(z):
                out={}
                for k in fields:
                    if k in ['cluster_id','gene_id','gene_symbol','measurement_status','scale_reference_branch','scale_definition']:out[k]=z[k]
                    elif k=='keep_for_plot':out[k]=boolean(z[k])
                    else:out[k]=number(z[k])
                return out
            ar,wr=keep(a),keep(w); cr={k:ar[k] for k in ['cluster_id','gene_id','gene_symbol','keep_for_plot']}
            for k in ['n_nuclei','sum_counts','n_detected','sum_log1p_normalized','n_mice_with_nuclei','n_mice_eligible_for_expression_review','n_mice_expression_supported']:
                cr[k]=ar[k]-wr[k]
                if cr[k]<-1e-7:raise ValueError('Negative cKO difference from produced aggregates')
                if cr[k]<0:cr[k]=0
            n=cr['n_nuclei'];cr.update({'mean_log1p_normalized':cr['sum_log1p_normalized']/n if n else None,'fraction_detected':cr['n_detected']/n if n else None,'mean_source_counts':cr['sum_counts']/n if n else None,'scaled_expression_z':None,'pct_difference_vs_rest':None,'mean_log1p_difference_vs_rest':None,'measurement_status':'PRODUCED_POOLED_MINUS_PRODUCED_WT_AGGREGATES;NO_NEW_RNA_AGGREGATION'})
            cr.update({'scale_reference_branch':None,'scale_reference_mean':None,'scale_reference_sd':None,'scale_definition':'UNAVAILABLE; cKO z not computed'})
            payload['rows'].append(ar);payload['WT_rows'].append(wr);payload['cKO_rows'].append(cr)
        if current_d is not None:save_markers(current_d,payload)
        if completed != expected:raise ValueError('Marker datasets incomplete')
        memberships=list(rows(self.bind(self.root/'annotation_initial/panel_membership.tsv')));groups=defaultdict(list);provided=[]
        for z in memberships:
            if z['mapping_status']=='MAPPED_STABLE_ID':
                groups[z['group']].append(z['gene_id'])
                if '.xlsx:' in z['source']:provided.append(z['gene_id'])
        mapped=[z for z in memberships if z['mapping_status']=='MAPPED_STABLE_ID']
        panel_counts={'source_panel_rows':len(memberships),'mapped_rows':len(mapped),'unique_mapped_stable_IDs':len({z['gene_id'] for z in mapped}),'unmapped_rows':len(memberships)-len(mapped),'workbook_only_mapped_rows':sum('.xlsx:' in z['source'] for z in mapped),'workbook_only_unique_stable_IDs':len(set(provided))}
        gene_sets=[{'id':'provided_all','label':'Reviewed supplied/canonical panel union (139 stable genes)','gene_ids':sorted({z['gene_id'] for z in mapped}),'source':'147 review-panel rows;145 mapped rows;139 unique stableIDs including canonical diagnostics;2 unmapped AC rows remain audited'}, {'id':'workbook_only','label':'Mapped workbook genes only','gene_ids':sorted(set(provided)),'source':'exact workbook rows without added canonical diagnostic genes'}]
        gene_sets += [{'id':'panel_'+re.sub('[^A-Za-z0-9_-]','_',g),'label':g,'gene_ids':sorted(set(v)),'source':'supplied/canonical set membership; expression evidence does not establish anatomy or morphology'} for g,v in groups.items()]
        rep=defaultdict(list);rep_rows=[]
        for z in rows(self.bind(self.root/'marker_summaries/representative_markers.tsv')):
            if z['branch']=='original_full':
                rep[(z['tool'],z['population'])].append(z['gene_id']);rep_rows.append({k:z[k] for k in ['tool','population','cluster_id','gene_id','gene_symbol','representative_rank','selection_scope']})
        gene_sets += [{'id':'representative_'+t+'_'+p,'label':t+' '+p+' representative marker union','tool':t,'population':p,'gene_ids':sorted(set(v)),'source':'existing native representative selection; original groups frozen to paired counts'} for (t,p),v in rep.items()]
        downloads=[{'label':'Annotated primary UMAP and marker book: '+z['id'],'path':z['pdf_path'],'scope':'primary original_frozen'} for z in datasets]
        for p,label in [(Path(self.cfg['study_root'])/'reviews/CA1_programme_agreement_v2_0_0/CA1_programme_agreement_and_threshold_sensitivity.pdf','CA1 programme semantic agreement and WT-threshold sensitivity'),(Path(self.cfg['study_root'])/'reports/visual_guide_v4/CA1_technical_diagnostics.pdf','CA1 saved-input technical controls'),(Path(self.cfg['study_root'])/'reports/visual_guide_v4/count_branches_and_filtering.pdf','Count branches and shared analysis masks'),(Path(self.cfg['study_root'])/'reports/visual_guide_v4/crosswalks/CA1_tool_cluster_crosswalks.pdf','CA1 atomic-cluster crosswalks'),(Path(self.cfg['study_root'])/'reports/visual_guide_v4/programs/CA1_ProS_expression_programs_all_tools_branches.pdf','CA1 expression programme UMAPs'),(Path(self.cfg['study_root'])/'reports/visual_guide_v4/programs/Inhibitory_expression_programs_all_tools_branches.pdf','IN expression programme UMAPs')]:
            if p.exists():asset=self.copy_asset(p,'assets/pdf/diagnostics/'+p.name);downloads.append({'label':label,'path':asset['path'],'scope':'read-only completed metadata diagnostic; primary DE unchanged'})
        aliases={z['gene_symbol']:z['gene_id'] for z in genes}
        manifest={'schema_version':SCHEMA,'study':{'title':'SATB2 WT/cKO snRNA-seq: offline reviewed comparison','source_run':str(self.root.parent),'source_fingerprint':sha(self.cfg_path),'notes':['Twelve independent mice: six WT and six cKO.','Source-verified existing MapMyCells labels are a shared reference, not independent truth.','Expression-associated programmes remain provisional; physical layers UNDETERMINED.','Gene UMAP color is a produced cluster mean, not measured per-nucleus gene expression.','Baseline excludes small All/EX groups through the shared iterative mask; IN and CA1 have no 100-nucleus cutoff.','Offline edits do not fit new statistics; an exact source-bound membership is required to reuse a produced DESeq2 endpoint.']},'genes':genes,'populations':[{'id':p,'label':{'allcells':'All nuclei','Excitatory':'Excitatory','Inhibitory':'Inhibitory','CA1_ProS':'CA1-ProS candidates'}[p]} for p in POPS],'branches':[{'id':p,'label':{'original_full':'Original full cohort','paired_original':'Paired original counts','paired_cellbender':'Paired CellBender counts'}[p]} for p in BRANCHES],'samples':[{'sample':s,'condition':c} for s,c in sorted(sample_map.items())],'datasets':datasets,'shards':self.shards,'gene_sets':gene_sets,'representative_markers':rep_rows,'downloads':downloads,'qc_fields':[{'key':k,'label':l,'type':('category' if k=='qc_doublet_call' else 'boolean' if k in ['qc_strict_sensitivity_pass','qc_high_genes_review'] else 'number')} for k,l in QC],'defaults':{'population':'CA1_ProS','branch':'original_full','basis':'harmony','annotation':'reviewed_label','marker_ids':[aliases[s] for s in ['Calb1','Reln','Etv1','Zbtb20','Camk2a','Slc17a7','Satb2'] if s in aliases]},'panel_mapping_counts':panel_counts,'de_catalog_status':'separate source-bound catalog registration; not a new fit','measurement_contract':{'gene_overlay':'cluster_mean_only','source_nucleus_expression_extracted':False,'new_RNA_aggregation_performed':False,'new_scientific_fits_performed':False,'physical_layers_validated':False}}
        (self.out/'data/asset_map.json').write_text(json.dumps(self.assets,indent=2)+'\n')
        (self.out/'data/manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
        (self.out/'data/manifest.js').write_text('window.SATB2_DATA='+json.dumps(manifest,ensure_ascii=False,separators=(',',':'),allow_nan=False)+';\n')
        self.sources[str(self.cfg_path)]=sha(self.cfg_path)
        for p,h in self.sources.items():
            if sha(p)!=h:raise ValueError('Source changed while packaging '+p)
        # Only package-owned transient annotation splits are removed, not the fit index.
        for p in annodir.glob('*.tsv'):p.unlink()
        annodir.rmdir()
        receipt={'status':'COMPLETED_DATA_BUNDLE_BASE_REVIEW_REQUIRED','schema_version':SCHEMA,'config_sha256':sha(self.cfg_path),'source_code_sha256':sha(__file__),'input_sha256':self.sources,'n_datasets':len(datasets),'n_genes':len(genes),'panel_mapping_counts':panel_counts,'n_coordinate_rows':sum(z['n_nuclei'] for z in datasets),'coordinate_index_path':str(index),'coordinate_index_sha256':sha(index),'source_count_validation_PASS':True,'scientific_fitting_performed':False,'raw_nucleus_count_aggregation_performed':False,'cKO_summary_arithmetic_from_produced_aggregates':True,'per_nucleus_gene_values_available':False,'shard_compressed_bytes':sum(z['gzip_bytes'] for z in self.shards.values()),'output_sha256':{str(p.relative_to(self.out)):sha(p) for p in [self.out/'data/manifest.json',self.out/'data/manifest.js',index]}}
        (self.out/'DATA_BUNDLE_BASE_COMPLETE.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({k:receipt[k] for k in ['status','n_datasets','n_genes','n_coordinate_rows','shard_compressed_bytes']}),flush=True)

    def augment_assets(self):
        receipt_path=self.out/'DATA_BUNDLE_BASE_COMPLETE.json'
        receipt=json.loads(receipt_path.read_text())
        if receipt['status']!='COMPLETED_DATA_BUNDLE_BASE_REVIEW_REQUIRED' or receipt['config_sha256']!=sha(self.cfg_path):raise ValueError('Base data receipt mismatch')
        for relative,h in receipt['output_sha256'].items():
            if sha(self.out/relative)!=h:raise ValueError('Base packaged data changed '+relative)
        self.complete('plots','COMPLETE')
        manifest=json.loads((self.out/'data/manifest.json').read_text());self.assets=json.loads((self.out/'data/asset_map.json').read_text())
        paths={x['path'] for x in manifest['downloads']};nsecondary=0
        for z in rows(self.bind(self.root/'plots/PLOT_BOOK_INDEX.tsv')):
            if z['cluster_scope']=='fresh_branch':
                d=ident(z['tool'],z['branch'],z['population']);relative='assets/pdf/secondary_fresh/'+d+'.pdf'
                if relative not in paths:
                    a=self.copy_asset(z['pdf'],relative);manifest['downloads'].append({'label':'Secondary fresh-cluster view (download only): '+d,'path':a['path'],'scope':'fresh_branch secondary visualization; original frozen inference labels remain authoritative'});paths.add(relative)
                nsecondary+=1
        p=Path(self.cfg['study_root'])/'reports/ca1_technical_diagnostics_v3/CA1_technical_diagnostics.pdf';relative='assets/pdf/diagnostics/CA1_technical_diagnostics_v3.pdf'
        if p.exists() and relative not in paths:
            a=self.copy_asset(p,relative);manifest['downloads'].append({'label':'Full CA1 fixed-input technical controls and graph-conversion diagnostics','path':a['path'],'scope':'completed native saved-input diagnostic; primary partitions unchanged'})
        if nsecondary!=16:raise ValueError('Expected 16 secondary native books')
        # Keep one canonical source-binding entry per portable asset.
        self.assets=list({x['bundle_path']:x for x in self.assets}.values())
        (self.out/'data/asset_map.json').write_text(json.dumps(self.assets,indent=2)+'\n')
        (self.out/'data/manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
        (self.out/'data/manifest.js').write_text('window.SATB2_DATA='+json.dumps(manifest,ensure_ascii=False,separators=(',',':'),allow_nan=False)+';\n')
        receipt['asset_supplement_source_code_sha256']=sha(__file__);receipt['input_sha256'].update(self.sources);receipt['n_primary_pdf_books']=24;receipt['n_secondary_download_only_pdf_books']=16;receipt['n_asset_downloads']=len(manifest['downloads'])
        for name in ['data/manifest.json','data/manifest.js','data/asset_map.json']:receipt['output_sha256'][name]=sha(self.out/name)
        receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({'asset_supplement':'COMPLETE','primary_books':24,'secondary_download_only_books':16,'downloads':len(manifest['downloads'])}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--outdir',required=True);ap.add_argument('--assets-only',action='store_true');args=ap.parse_args();b=Builder(args.config,args.outdir)
    if not args.assets_only:b.build()
    b.augment_assets()

