"""Explicitly fabricated layout fixture. No study inputs or inference are used."""
import importlib.util
from pathlib import Path
import json
import numpy as np
import pandas as pd

path=Path(__file__).with_name('build_effect_figures.py');spec=importlib.util.spec_from_file_location('figure_fixture',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
out=path.parents[1]/'reviews/effect_figure_SYNTHETIC_LAYOUT_v3'
if out.exists() and any(out.iterdir()):raise RuntimeError('Preserve prior fixture; use a new version')
out.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(20261014);effects={};coverage={};robust={};statuses=[];core=[];selections=[];mouse=[];cross=[];comps=[]
samples=['wt_1122','wt_1123','wt_1178','wt_1220','wt_1226','wt_1228','cKO_1138','cKO_1139','cKO_1160','cKO_1165','cKO_933','cKO_991'];conditions=['WT']*6+['cKO']*6
for pop in m.POPS:
    gids=[f'FAKE_{pop}_{i:04d}' for i in range(500)]
    lfc=rng.normal(0,.65,len(gids));lfc[:4]=[1.8,1.6,1.4,1.2]
    padj=np.clip(10**(-abs(lfc)*4),0,1);padj[:4]=[1e-14,1e-12,1e-10,1e-8];padj[20:28]=np.nan
    robust[pop]=set(gids[:8]) if pop!='Inhibitory' else set()
    base=pd.DataFrame(dict(gene_id=gids,display_symbol=[f'FakeGene{i}' for i in range(500)],baseMean=100,log2FoldChange=lfc,lfcSE=.2,stat=lfc/.2,pvalue=padj/2,padj=padj))
    for branch in m.BRANCHES:
        z=base.copy()
        if branch=='paired_cellbender':z.log2FoldChange+=rng.normal(0,.12,len(z));z.loc[0:1,'padj']=np.nan
        effects[(pop,branch)]=z
        coverage[(pop,branch)]=pd.DataFrame(dict(sample=samples,condition=conditions,n_nuclei=rng.integers(400,1400,12)))
        for tool in ['scanpy','seurat']:
            for i in range(8):statuses.append(dict(tool=tool,branch=branch,population=pop,endpoint='PRIMARY',status='TESTED_REVIEW_REQUIRED' if i<6 else 'NOT_TESTED'))
    core.append(pd.DataFrame(dict(tool='scanpy',population=pop,gene_id=gids,PRIMARY_FDR05_all_three_count_branches_same_sign=[i<18 for i in range(500)],stringent_reliability_flag=[g in robust[pop] for g in gids])))
    for rank,gid in enumerate(sorted(robust[pop])[:2],1):
        row=base[base.gene_id.eq(gid)].iloc[0];selections.append(dict(population=pop,gene_id=gid,gene=row.display_symbol,selection_rank=rank,log2FoldChange=row.log2FoldChange,padj=row.padj))
        for sample,cond in zip(samples,conditions):mouse.append(dict(population=pop,gene_id=gid,sample=sample,condition=cond,normalized_aggregate_UMIs=float(rng.uniform(600,1100)*(2 if cond=='cKO' else 1))))
    for i in range(4):
        key=dict(population=pop,left_population_key=f'S_{pop}_{i}',right_population_key=f'R_{pop}_{i}')
        cross.append(dict(**key,branch='original_full',analysis_kind='cluster',left_label=str(i),right_label=str(i+2),jaccard=.35+i*.16))
        comps.append(dict(**key,comparison_scope='FROZEN_POPULATION_DEFINITION_COMPARISON_SAME_DESEQ2_ENGINE',left_branch='original_full',right_branch='original_full',jaccard=.35+i*.16,effect_Spearman=.72+i*.07,fine_population_interpretation_guard_pass=i>1,left_coverage_flag='LOW_NUCLEUS_COVERAGE_EXPLORATORY' if i==0 and pop=='Inhibitory' else 'STANDARD_DECLARED_COVERAGE',right_coverage_flag='STANDARD_DECLARED_COVERAGE'))
plt,PdfPages=m.setup_plotting();writer=m.FigureWriter(out,plt,PdfPages,synthetic=True)
m.draw_status(plt,writer,pd.DataFrame(statuses));m.draw_coverage(plt,writer,coverage);m.draw_volcano(plt,writer,effects,robust)
metrics=m.draw_concordance(plt,writer,effects,robust);m.draw_fine(plt,writer,pd.DataFrame(cross),pd.DataFrame(comps));m.draw_reliability(plt,writer,effects,pd.concat(core,ignore_index=True));m.draw_mouse_genes(plt,writer,pd.DataFrame(mouse),pd.DataFrame(selections));writer.close()
z=m.paired_effect_frame(effects[('allcells','paired_original')],effects[('allcells','paired_cellbender')],set())
assert z.loc['FAKE_allcells_0000','plot_class']=='BH hit; counterpart FDR unavailable'
assert len(writer.rows)==7 and len(metrics)==4
(out/'SYNTHETIC_TEST_RECEIPT.json').write_text(json.dumps(dict(status='PASS',scope='FABRICATED_LAYOUT_ONLY',actual_study_inputs_read=False,scientific_fitting_performed=False,raw_nucleus_aggregation_performed=False,figure_count=7,source_code_sha256=m.sha(path),output_sha256={p.name:m.sha(p) for p in out.iterdir() if p.is_file()}),indent=2)+'\n')
print(out)
