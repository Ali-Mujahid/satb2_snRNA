# Frontend regression checks

Run with a modern Node.js runtime from this repository:

```sh
node tests/transport.cjs
node tests/catalog_pages.cjs
node tests/saved_semantics.cjs
```

These checks use JavaScript VM fixtures and the already-produced local binary payloads. They do not launch a browser, call the network, extract nucleus RNA or fit a statistical model. Set `SATB2_SITE_DIR` only if testing a different site directory; exact text downloads are read from its sibling `data_downloads/` directory. The saved-data check verifies all four populations and existing fit-selection/count/NA safeguards, so it needs the complete processed-object repository checkout. Test receipts are written beside the scripts.

Passing these tests does not establish hosted CORS, actual browser rendering, browser file download behavior or deployment success. The recorded actual-browser checklist is a separate acceptance step.
