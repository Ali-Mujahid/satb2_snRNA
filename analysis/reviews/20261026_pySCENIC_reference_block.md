# pySCENIC reference readiness — resource blocked

Full mouse motif-supported regulon inference has not run. The user selected full pySCENIC plus CollecTRI TF activity and hdWGCNA, then confirmed that no cached mouse ranking and motif annotation files are available and instructed us to continue the available analyses.

The isolated pySCENIC0.12.1/Python3.9 environment passed import-only checks. Both supported official cisTarget hosts resolved to the same IPv4 endpoint and refused TCP connections before TLS or HTTP. Normal HEAD and ranged GET requests failed; no scientific reference file was downloaded. Registered remote compute targets were also unreachable. These observations are recorded in references/pyscenic_mm10/blocked_resource_manifest.json. They establish a resource-access blocker, not a failure of the biological dataset.

Required references are the mm10 RefSeq-r80 v10-cluster gene-based rankings for 500bp-up/100bp-down and 10kbp-up/10kbp-down windows, matching MGI v10nr-cluster motif annotations, and the matching mouse TF list. The installed environment alone cannot provide these resources. No human, mm9, region-based or unverified alternate reference is substituted.

The remaining steps, once matching references can be acquired and checked, are a separately reviewed immutable GRNBoost2 → cisTarget motif-pruning → AUCell workflow. The public protocol mouse TF list available from GitHub omits Satb2; SATB2 coverage must be assessed explicitly rather than forcing an unsupported regulon. Absence of an inferred SATB2 regulon would not exclude its architectural regulatory roles.

CollecTRI activity is curated-target enrichment, and hdWGCNA is RNA coexpression. Neither is presented as completed de novo motif-supported regulon inference. The atlas must display this branch as blocked until native results exist.
