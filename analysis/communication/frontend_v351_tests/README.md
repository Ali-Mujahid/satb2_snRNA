# Historical v3.5.1 frontend tests

These files preserve the original broad-candidate model and application tests from
the `atlas-v3.5.1` source/tag. They describe the earlier six-fit explorer, including
OR endpoint DEG support and all-candidate views. They are historical source records,
not tests executable against the corrected v3.6 DEG-only endpoints or model.

The corrected active suites are `tests/communication_deg_ui.cjs`,
`tests/communication_deg_only.cjs`, and `tests/communication_app.cjs`.
`tests/communication_model.cjs` remains a compatibility entry point to the active
strict DEG-only invariant suite.

For historical reproduction, check out the original `atlas-v3.5.1` tag and use the
tests and paths in that checkout. Original scientific analyses and saved results
remain preserved; the corrected atlas changes eligibility and presentation.
