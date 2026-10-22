# Supplemental complex-circle review — completed final report

Reviewed the actual saved PNGs in
`runs/20261022_liana_complete/results/report/gene_networks`.
No code, figures, native evidence, or analysis was changed. The checks below
read completed figure exports only; no model or statistical fit was run.

## Visual review

**Pass: no material clipping or labeling defect found in the three inspected
panels.** The gene/complex names, population sectors, arrowheads, circle/square
role legend, direction colors, titles and interpretation footnotes are readable.

| Inspected saved panel | Edge rows | Endpoint-role nodes | Fixed selected identities |
|---|---:|---:|---:|
| `IN_clusters/CA1_to_IN/WT/complex_circle.png` | 12 | 8 | 12 |
| `IN_clusters/IN_to_CA1/cKO/complex_circle.png` | 12 | 8 | 12 |
| `IN_total/bidirectional/pooled/complex_circle.png` | 12 | 19 | 12 |

- CA1→IN WT has blue arrows from circular CA1 ligand nodes to square
  inhibitory receptor nodes. IN→CA1 cKO reverses the sender family and uses
  orange arrows pointing toward CA1 square receptor nodes. The bidirectional
  pooled panel contains both source colors with consistent circle/square roles.
- Node names are associated with explicit population sectors; a molecule that
  occurs in different groups or roles remains a separate endpoint node. Review
  asterisks remain visible, including IN C12 and inhibitory total.
- Captions accurately describe all-subunit detection areas, native LIANA
  geometric-mean expression widths, common scales, fixed cross-condition
  selection, and the absence of a genotype interaction test. No arrow color is
  interpreted as activating or inhibiting signaling.
- The top 12 cluster-level *interaction identities* repeatedly include
  Nlgn1/Nrxn1 at different cell-pair endpoints. This follows the documented
  source×target×ligand×receptor ranking and is not twelve different molecular
  pairs. The broader dotplots/full tables provide additional molecular breadth;
  the circle selection should not be described as twelve distinct LR species.

## Saved-export checks

For all 18 panels, the interaction-row and node-row counts exactly match the
manifest. Every edge's source and target node identifiers occur in its node
table; endpoint-role node keys are unique; no panel exceeds 40 nodes; all PNG,
PDF and SVG assets exist. Node identities and coordinates are identical across
WT, cKO and pooled contexts within each scope/direction, as required for
comparison. The reviewed cluster-direction panels use the recorded shared
geometric-mean score maximum 4.321054; the reviewed inhibitory-total
bidirectional set uses 3.7425296.

The 10 primary figures inspected in the earlier
`figures_visual_review.md` were also checked by SHA-256 against their copies in
`runs/20261022_liana_complete/results/report/figures`: every pair is
byte-identical. That earlier visual review therefore covers the corresponding
assets in the completed final report as well.

This is representative visual review and saved-export consistency checking,
not a claim that every PDF page was separately rendered or that predictions
were experimentally validated.
