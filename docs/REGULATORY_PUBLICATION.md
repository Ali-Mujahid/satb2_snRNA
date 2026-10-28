# COMPLETE-only regulatory publication transport

`scripts/build_regulatory_publication.py` packages completed hdWGCNA outputs. It
does not fit models, infer TF activity, run scientific workflows, push Git,
create tags/releases or deploy Pages. Do not invoke it against real outputs until
the approved native run has written a successful `FINAL_VALIDATION.json`.

Required gates include validated primary CA1 and IN networks, exact twelve-mouse
keys, finite native eigengenes, frozen memberships, native SHA256/MD5 receipts,
the executed immutable scientific source/config, and unchanged original bound
input hashes. Optional consensus limitations remain explicit. Curated TF activity
can be added only from its separately completed and verified native contract.
Motif-supported pySCENIC remains explicitly resource-blocked; nucleus activity
UMAPs are unavailable.

Use absolute paths and a **new nonexistent** deliverable directory:

```bash
python3 scripts/build_regulatory_publication.py \
  --results /absolute/approved/run/results \
  --source /absolute/immutable/scientific/source \
  --config /absolute/immutable/config.json \
  --outdir /absolute/new/publication_stage \
  --tag atlas-v3.7.0 \
  --stage-only
```

For the combined release, add **all three** optional arguments:

```bash
  --tf-results /absolute/approved/tf_run/results \
  --tf-source /absolute/immutable/tf/source \
  --tf-config /absolute/immutable/tf_config.json
```

The independent TF gate requires
`satb2-tf-regulatory-final-validation/v1`, status
`COMPLETE_WITH_BLOCKED_MOTIF`, exhaustive native/package SHA256 receipts,
executed source hashes, and unchanged bound inputs. Both components must retain
identical five frozen endpoint membership hashes/counts, the same twelve mice and
genotypes, and the same original raw-count and DESeq2 source bytes. Method,
feature and network-node IDs are disjoint; duplicate mouse-score keys are
rejected. Unmodelable TF activities retain explicit null/false eligibility rows.

The combined payload is generated only in the new publication stage. Native HD
and TF transports, files, source and configs remain byte-identical. Historical
HD `TF_NOT_RUN` provenance stays nested under its original component; the new
top provenance explicitly declares completed curated CollecTRI ULM activity and
blocked motif inference. Native Wald enrichment and self-edge/mouse sensitivity
supplements remain distinct from primary mouse-score inference.

An explicit copy-only export recovery is accepted without relaxing scientific
completion gates. Its configuration must bind the original native configuration
and a frozen exhaustive `prepared/` plus `hdwgcna/` inventory. Only `run_root`,
`workflow_version`, and `reuse_completed_native` may differ; every other input and
scientific parameter must remain identical. The builder independently checks
every original and recovered SHA256, all four completed checkpoint MD5 seals,
the original executed fitting source, and the no-refit copy receipt. The original
configuration, inventory and fitting-source receipts are included in the ZIP;
provenance distinguishes the original fits from downstream export recovery.

`--stage-only` prepares everything without modifying the repository site. Omit it
only after review to copy the verified staged regulatory subtree into
`site/regulatory/`. Existing nuclei/communication pages, their links and scientific
sources are not edited. A new output directory is required for each build.

The deliverable contains:

- `site/regulatory/`: compressed `data.js`, compact manifest, frontend, native
  figure gallery, table/archive download index, and exact PNG/PDF files that fit.
- `release_assets/`: one standalone HTML, exact native/source ZIP, and any native
  PNG/PDF overflow assets. Every required upload is listed in `RELEASE_ASSETS.json`.
- The native ZIP includes selected raw MTX/cells/genes, RDS state, every native
  PDF/PNG/SVG and table, scientific source/config, setup scripts/manifests and
  unchanged DESeq2/LR reference inputs. Original full H5AD/annotation inputs and
  installed package libraries remain external hash-bound dependencies.
- With TF inputs, the same ZIP additionally includes every TF result under
  `tf_component/results/`, exact scientific source/configuration under
  `tf_component/source/` and `tf_component/config/`, all TF bound reference
  inputs, and the new merged JSON under `combined_transport/`. The publication
  merge script is preserved under `publication_source/`. Both components' native
  figure collections appear in the gallery with separate source labels.
- `PUBLICATION_RECEIPT.json`, `ARCHIVE_CONTENTS.json`, native file receipts and
  portable HTML receipts preserve source identity and transport hashes.

New Pages files are limited to **30,000,000 bytes**, and the prospective complete
site to **950,000,000 bytes**. Raw JSON, duplicate gzip, RDS, MTX and native SVG
stay outside Pages in the release ZIP. Native PNGs are preferred for Pages
thumbnails and displayed with CSS; they are not resized or altered. Overflow
figures retain direct, versioned release-asset links. If the complete ZIP exceeds
the per-asset budget it is staged as ordered byte parts with reconstruction
instructions; the download index clearly requires every part.

Validate before any publication:

```bash
python3 scripts/build_regulatory_publication.py --verify-outdir /absolute/publication_stage
node tests/regulatory.cjs
PYTHONDONTWRITEBYTECODE=1 python3 tests/regulatory_publication.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/regulatory_verification.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/regulatory_recovery.py
python3 scripts/verify_site.py
python3 scripts/verify_communication_atlas.py
python3 scripts/verify_regulatory_atlas.py --site-dir /absolute/publication_stage/site/regulatory
```

Every ZIP member is streamed back and checked against its native SHA256 (which
also checks ZIP CRC). The archive contains `VALIDATE_EXTRACTED.py` for validation
after extraction. Native/source/config and external bound inputs are rehashed
after staging to confirm they were not changed. The standalone HTML embeds CSS,
fflate, jsPDF, data and app; navigation points to absolute HTTPS pages, and no
external rendering dependencies remain.

Upload and verify **all** assets listed in `RELEASE_ASSETS.json` before deploying
the Pages content. The builder does not upload them. The original payload's
relative scientific download paths are not exposed as working UI links; native
tables and figures are routed through the separate verified gallery/archive
index. Do not publish a waiting page or synthetic fixture. Synthetic transport
tests are isolated outside `site/` and outside real native run directories.
