# Frontend regression checks

Run with a modern Node.js runtime from this repository:

```sh
node tests/transport.cjs
node tests/catalog_pages.cjs
node tests/number_format.cjs
node tests/overlay_selection.cjs
node tests/category_persistence.cjs
node tests/saved_semantics.cjs
```

These checks use JavaScript VM fixtures and the already-produced local binary payloads. They do not launch a browser, call the network, extract nucleus RNA or fit a statistical model. Set `SATB2_SITE_DIR` only if testing a different site directory; exact text downloads are read from its sibling `data_downloads/` directory. The saved-data check verifies all four populations and existing fit-selection/count/NA safeguards, so it needs the complete processed-object repository checkout. Test receipts are written beside the scripts.

Passing these tests does not establish hosted CORS, actual browser rendering, browser file download behavior or deployment success. The recorded actual-browser checklist is a separate acceptance step.

The number-format regression covers tiny scientific-notation exponents without trimming exponent digits. The overlay-selection regression uses native select-value semantics to retain a selected Calb1 and keep current choices visible outside search results. The corrected portable archive is built separately by `scripts/build_corrected_offline_runtime.py`; original runtime and scientific payloads remain preserved.

The category-shortcut check verifies that all five category buttons save population/programme state before asynchronous loading and restore that choice after a simulated reload.
