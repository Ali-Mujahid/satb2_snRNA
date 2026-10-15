"""Fabricated DataFrame-only contract checks; no study/model/count input files."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
import pandas as pd

source=Path(__file__).with_name("build_effect_figures.py")
spec=importlib.util.spec_from_file_location("effect_reader_fixture",source)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.BRANCHES=["original_full"]
samples=[f"FAKE_WT_{i}" for i in range(6)]+[f"FAKE_KO_{i}" for i in range(6)]
conditions=["WT"]*6+["cKO"]*6
rows=pd.DataFrame([
  dict(tool="scanpy",population="Inhibitory",analysis_kind="cluster",analysis_label="FAKE_LOW",population_key="FAKE_LOW",model_enabled=True),
  dict(tool="scanpy",population="Inhibitory",analysis_kind="cluster",analysis_label="FAKE_EMPTY",population_key="FAKE_EMPTY",model_enabled=False)])
ids=["FAKE_GENE_1","FAKE_GENE_2","FAKE_GENE_3"]
z=pd.DataFrame(dict(gene_id=ids,display_symbol=["FakeSymbol","FakeSymbol",np.nan],pvalue=[.01,.2,np.nan],padj=[.03,.25,np.nan]))
coverage={}
prefilters={}
audits={}
for key in rows.population_key:
    ns=[1,2,4,0,0,0]*2 if key=="FAKE_LOW" else [0]*12
    coverage[key]=pd.DataFrame(dict(sample=samples,condition=conditions,n_nuclei=ns,total_UMIs=[n*100 for n in ns],included_in_model=[n>0 for n in ns]))
    prefilters[key]=pd.DataFrame(dict(gene_id=ids,count_cutoff_per_mouse=5,required_mice=3,passes_primary_prefilter=[key=="FAKE_LOW"]*3))
    audits[key]=dict(status="TESTED_REVIEW_REQUIRED" if key=="FAKE_LOW" else "NOT_TESTED",n_WT=3 if key=="FAKE_LOW" else 0,n_cKO=3 if key=="FAKE_LOW" else 0,residual_df=4 if key=="FAKE_LOW" else 0,n_frozen_primary_genes=3 if key=="FAKE_LOW" else 0,min_included_mouse_nuclei=1,max_included_mouse_nuclei=4,coverage_flag="LOW_NUCLEUS_COVERAGE_EXPLORATORY" if key=="FAKE_LOW" else "EXCLUDED_BY_PRESPECIFIED_ATOMIC_MASK",warning_count=0,execution_origin="FABRICATED_ONLY")

class FakeEvidence:
    def table(self,path,**kwargs):
        key="FAKE_LOW" if "FAKE_LOW" in str(path) else "FAKE_EMPTY"
        if path.name=="MICE_ALL12_AUDIT.tsv":return coverage[key].copy()
        if path.name=="FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv":return prefilters[key].copy()
        raise AssertionError("Unexpected fixture read")
data=SimpleNamespace(indices={"original_full":rows},model_roots={"original_full":Path("/FAKE_NO_INPUTS")},rosters={"original_full":dict(zip(samples,conditions))},
    endpoint=lambda branch,row:(z.copy() if row.population_key=="FAKE_LOW" else None,audits[row.population_key].copy()))
module=SimpleNamespace(truth=m.truth)
review=m.review_all_primary_exports(FakeEvidence(),data,module)
low=review[review.population_key.eq("FAKE_LOW")].iloc[0]
assert low.included_mice_below5_nuclei==6 and low.zero_nucleus_mice==6
assert low.rows_with_duplicated_symbol==2 and low.missing_symbol_rows==1
assert low.n_pvalue_available==2 and low.n_FDR_available==2 and low.n_BH05==1
assert review.loc[review.population_key.eq("FAKE_EMPTY"),"n_modeled_genes"].iloc[0]==0
audits["FAKE_LOW"]["coverage_flag"]="STANDARD_DECLARED_COVERAGE"
try:m.review_all_primary_exports(FakeEvidence(),data,module)
except ValueError as ex:assert "coverage flag" in str(ex)
else:raise AssertionError("Missing exploratory IN flag was accepted")
audits["FAKE_LOW"]["coverage_flag"]="LOW_NUCLEUS_COVERAGE_EXPLORATORY"
prefilters["FAKE_LOW"]["required_mice"]=5
try:m.review_all_primary_exports(FakeEvidence(),data,module)
except ValueError as ex:assert "prefilter policy" in str(ex)
else:raise AssertionError("Wrong relaxed IN gene coverage was accepted")
prefilters["FAKE_LOW"]["required_mice"]=3
coverage["FAKE_EMPTY"].loc[0,"sample"]="FAKE_WRONG_MOUSE"
try:m.review_all_primary_exports(FakeEvidence(),data,module)
except ValueError as ex:assert "identities differ" in str(ex)
else:raise AssertionError("Wrong mouse identity in disabled endpoint was accepted")
print(json.dumps(dict(status="PASS",scope="FABRICATED_DATAFRAMES_ONLY",source_code_sha256=m.sha(source),actual_study_inputs_read=False,scientific_fitting_performed=False,raw_nucleus_aggregation_performed=False,guards=["IN_min1_all_nonempty","IN_below5_exploratory_flag","IN_gene_cutoff5_required3","empty_registry_preserved","missing_duplicate_symbol_and_FDR_coverage"]),indent=2))
