# Read-only scientific result review

Completed native output reviewed: `runs/20261022_liana_v1/results/communication/`.
No model was fitted or rerun by this review. Full validation receipt:
`native_output_scientific_review.json` — all 43 independent checks pass.

## What was verified

Across 64,936 native interaction records and 226,930 independent-mouse score
records, the four interaction keys are unique within scope/condition; native
means, expression products, geometric means and group-specific logFCs retain
their raw numerical meanings; 1,000-permutation CellPhoneDB P values lie on the
expected 0.001 grid. Every endpoint-specific DEG support flag was reconstructed
independently from the original stable-ID DESeq2 source tables and matches.
Per-genotype supporting-mouse and eligible-mouse counts were independently
recounted from their native source tables and match. C18 remains untested and
never contributes a gene-level DEG flag. IN-total likewise contributes no
resource-overlapping DEG flag, for the different reason that its four DEGs do
not occur in this protein LR resource.

Pooled/WT/cKO cluster-resolution expression-qualified records: 20,471 / 19,467 /
21,680. IN-total resolution: 1,108 / 1,040 / 1,170. These are candidate records
passing expression gates. Their differences are not numbers of significantly
changed interactions.

A plotting recovery that accepts 1e-9 numerical tolerance and clips detection
fractions to [0,1] is appropriate for tiny sparse-mean roundoff (reported maximum
excess 1.18e-13). No raw statistical result or fit needs to change for that repair.

## Gene-level evidence is mostly on the CA1 side

The inhibitory-total DEGs have no ligand/receptor overlap in mouseconsensus.
Only inhibitory clusters C1 and C6 have DEGs overlapping LR genes. Their LR
subunits are:

- C1: Cd200, Col16a1, Col4a1, Fstl5, Grm7, Nlgn1, Ptprz1.
- C6: Cdh7, Cntn3, Robo1, Sema6d, Slit2.

CA1 deep has 742 primary unambiguous DEG-supported genes. The original CA1
superficial export has 1,825 DEGs; two ambiguous symbols are correctly excluded
from primary interaction DEG support, leaving 1,823 unambiguous gene records.
These are existing DESeq2 results, not new differential fits.

With pooled magnitude_rank<=0.05, DEG support, at least 3 supporting WT and 3
supporting cKO mice, and an inhibitory endpoint marked REFERENCE_SUPPORTED,
there are 1,126 directed candidate rows: 1,102 with CA1-only gene-level DEG
support, 16 with inhibitory-only support, and 8 with both-side support. An
additional native CellPhoneDB specificity P<=0.05 filter retains 734: 724
CA1-only, 7 inhibitory-only, and 3 both-side. All 10 remaining rows with any
inhibitory-side DEG support involve C6. This optional specificity filter is not
a genotype-change test or a newly adjusted FDR.

The pooled primary-quadrant counts before optional CellPhoneDB specificity
filtering are in `scientific_candidate_counts.tsv`.

## A balanced 16-record short list

`curated_16_replicated_candidates.tsv` provides four diverse examples for each
CA1-program/direction combination, including exact source and target, ligand
and receptor, native ranks and expression, matched DESeq2 stable ID, signed
log2(cKO/WT), padj, and mouse support. Every selected record is expression
supported in all 6 WT and all 6 cKO mice. Selection is explicitly curated for
diversity; it is not a new statistical score or an exhaustive top-rank list.

Most useful observations to describe carefully:

1. Deep CA1 Slit1 -> C6 Robo1: Slit1 is lower in cKO (log2FC -0.30242,
   padj 0.0191792), as is C6 Robo1 (-0.458227, 3.38377e-05). This is a candidate
   with reduced RNA for both interacting components, not proof that signaling
   activity is reduced. CA1 Slit3 -> C6 Robo1 is also strongly expression
   supported, but its DEG evidence is the inhibitory receptor only.
2. C6 Cntn3 -> CA1 Ptprg: Cntn3 is lower (-0.538346, 0.00370056), while Ptprg
   is higher in deep (+0.396663, 0.0164891) and superficial (+0.490437,
   2.92705e-08) programs. Opposite component directions prevent assigning a
   single up/down activity label to this interaction.
3. CA1 Ntng1 -> C1 Lrrc4c: Ntng1 is higher in deep (+1.01493, 3.40261e-25)
   and superficial (+1.49741, 4.10692e-64) programs. The selected C1 receptor
   is not independently a DEG; the candidate has CA1-side evidence only.
4. C2 Efna5 -> CA1 Epha6: Epha6 is lower in deep (-0.327446, 0.000150768)
   and superficial (-0.500858, 7.16316e-11). The DEG evidence is again the
   CA1 endpoint only.
5. C15 Ntng1 -> superficial Lrrc4c: CA1 Lrrc4c is lower (-0.818584,
   1.19752e-29), whereas C15 Ntng1 is not a DEG. This must not be confused
   with the opposite CA1 Ntng1 -> IN receptor orientation above.

## Interpretation boundaries and lower-priority records

Nrxn/Nlgn pairs can have exceptionally high expression/magnitude while lacking
cell-group specificity in some endpoints (native CellPhoneDB P=1 and aggregate
specificity_rank=1). Preserve and display magnitude versus specificity
separately. Curated reversed Nrxn/Nlgn orientations can also have exactly the
same scores; they are not independent evidence of two physiological signaling
routes. Native zero permutation P values mean below empirical resolution,
not an exact zero probability.

Do not promote Rims1 -> Cacna1c from a resource-ranked row to a claim of a
secreted extracellular ligand. Rims1 is a cytoplasmic presynaptic active-zone
scaffold; full native exports should retain the row, while extracellular
interpretation should flag its localization. Curated primary references:
https://www.ncbi.nlm.nih.gov/datasets/gene/116837/
https://www.informatics.jax.org/homology/GOGraph/Rims1

C18 is DE-untested. C20 has an identity-review flag including immune-associated
signatures and is not used in the reference-supported shortlist. Mixed/provisional
inhibitory groups remain in complete expression results without being promoted
to confident fine identities. CA1 deep/superficial are molecular expression
programs, not demonstrated anatomical layer assignments.

All identified candidates remain transcript- and prior-supported hypotheses.
Neither rank, detection, native cell-label permutation P values, nor component
DESeq2 evidence establishes cell contact, protein availability, receptor
activation/inhibition, synaptic connectivity, or a statistically tested change
in interaction activity between genotypes.

## Additional, stronger consensus-specificity tier

The curated 16-record table above applies native CellPhoneDB specificity
P<=0.05 but **does not** apply aggregate specificity_rank<=0.05. These must not
be conflated. In particular, Slit3 -> C6 Robo1 has native CellPhoneDB P below
resolution but aggregate specificity_rank approximately 0.166; it is excluded
from the stronger tier below.

Starting from the same 1,126 pooled replicate/identity/DE/magnitude-supported
records, native CellPhoneDB P<=0.05 retains 734. Requiring aggregate
specificity_rank<=0.05 as well retains **324** directed records: 112 CA1 -> IN
and 212 IN -> CA1. Within the 324, 322 have CA1-only DEG evidence and two have
both-side DEG evidence; none has inhibitory-only DEG evidence. The two both-side
records are C6 Cntn3 -> deep Ptprg and C6 Cntn3 -> superficial Ptprg, with opposite
signed component DESeq2 effects as described above.

Exact counts by program/direction for all three nested tiers are exported in
`candidate_tier_counts.tsv`. The complete 324-record stronger-tier table is
`strong_consensus_specificity_candidates.tsv`. This is an explicitly defined
post-inference prioritization filter, not a genotype test or claim that excluded
candidates do not interact.

Some Nrxn/Nlgn orientations are genuinely well supported by the native
specificity metrics, including superficial Nlgn1 -> C3 Nrxn1 (aggregate
specificity_rank ~0.00272), while other groups are context-weak. The caveat is
endpoint-specific and does not dismiss the full gene family. C6 is described as
Sst-reference-supported; its fine backprojection resemblance is not proven
morphology or connectivity.
