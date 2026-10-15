#!/usr/bin/env python3
"""Produced-table synthetic checks; no study inputs or scientific fitting."""
import unittest
import numpy as np
import pandas as pd
import compare_revision as c


def table(effects, fdr=None, p=None):
    n = len(effects)
    return pd.DataFrame({'gene_id':[f'g{i}'for i in range(n)],'baseMean':np.full(n,10.),'log2FoldChange':effects,
                         'lfcSE':np.ones(n),'stat':effects,'pvalue':p if p is not None else np.full(n,.001),
                         'padj':fdr if fdr is not None else np.full(n,.01)})


class ComparisonChecks(unittest.TestCase):
    def test_effect_NA_availability_never_false_nonsignificance(self):
        a=table([1.,2.,-3.],[.01,.06,np.nan],[.01,.02,np.nan])
        b=table([1.,-2.,-3.],[.02,np.nan,.01],[.01,np.nan,.4])
        z=c.effect_metrics(a,b)
        self.assertEqual(z['common_finite_FDR'],1)
        self.assertEqual(z['common_finite_pvalues'],1)
        self.assertAlmostEqual(z['effect_same_sign_fraction'],2/3)
        self.assertEqual(z['FDR05_Jaccard_on_common_available_rows'],1.)
        self.assertFalse(z['statistical_implementation_independence'])

    def test_rectangular_member_matching_and_support(self):
        left=pd.DataFrame({'population_key':['a','b']});right=pd.DataFrame({'population_key':['x','y','z']})
        m={'a':set(range(40)),'b':set(range(40,80)),'x':set(range(35))|set(range(40,45)),
           'y':set(range(35,40))|set(range(45,80)),'z':set(range(80,85))}
        z=c.crosswalk(left,right,m)
        self.assertEqual(sum(x['shared_nuclei']for x in z),70)
        self.assertTrue(all(x['match_pass']for x in z))
        self.assertTrue(all(x['jaccard']>.7 for x in z))

    def test_two_coarse_genes_do_not_certify_mixed_fine_labels(self):
        annotations=pd.DataFrame({'tool':['scanpy']*3,'branch':['original_full']*3,'population':['Inhibitory']*3,
                                 'cluster_scope':['original_frozen']*3,'cluster_id':['0','1','2'],
                                 'atlas_broad_majority_fraction':[1.]*3,'atlas_broad_high_support_fraction':[1.]*3,
                                 'atlas_majority_fraction':[.99,.51,.99],'atlas_subclass_high_support_fraction':[1.]*3,
                                 'identity_review_status':['REFERENCE_SUPPORTED','REFERENCE_MIXED_SUBCLASSES','FULL_COHORT_CONFLICTING_MARKER_SIGNATURE'],
                                 'n_nuclei':[400,400,400]})
        registry=annotations[['tool','branch','population','cluster_scope','cluster_id']].copy()
        registry['n_supported_representative_genes']=2;registry['marker_review_status']='SUPPORTED'
        z=c.label_reliability(annotations,registry)
        self.assertEqual(z.fine_label_evidence_supported.tolist(),[True,False,False])
        self.assertFalse(z.physical_layer_validated.any())

    def test_stringent_requires_all_branches_partitions_and_12LOO(self):
        row=pd.Series(dict(tool='scanpy',population='Inhibitory',analysis_kind='total',category_family='',analysis_label='TOTAL',population_key='P1'))
        class Fake:
            rosters={'original_full':{f'm{i}':'WT'if i<6 else'cKO'for i in range(12)}}
            families={'paired_original':{c.identity(row):row},'paired_cellbender':{c.identity(row):row}}
            def endpoint(self,branch,row,label='PRIMARY'):
                z=table([1.,2.,3.,4.,5.]);ad={'status':'TESTED_REVIEW_REQUIRED','warning_count':0}
                if branch=='paired_cellbender':z.loc[1,'log2FoldChange']=-2.;z.loc[4,'padj']=np.nan
                if label=='KNOWN_PARTITION_ADJUSTED':z.loc[3,'log2FoldChange']=-4.
                if label=='LOO_omit_m3':z.loc[2,'log2FoldChange']=np.nan
                return z,ad
        z,scope=c.core_reliability(Fake(),row)
        self.assertEqual(z.stringent_reliability_flag.tolist(),[True,False,False,False,False])
        self.assertEqual(z.loc[2,'n_LOO_effect_estimable'],11)
        self.assertEqual(z.loc[0,'min_LOO_log2FoldChange'],1.)
        self.assertEqual(scope['n_stringent'],1)


if __name__=='__main__':unittest.main()
