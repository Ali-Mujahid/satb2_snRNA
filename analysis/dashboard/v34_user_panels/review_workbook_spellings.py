#!/usr/bin/env python3
"""Read-only spelling candidates from preserved user workbook copies.

This writes review evidence only, never resolved panels or workbook edits.
Spelling similarity is not evidence of gene identity and is never applied.
"""
from collections import defaultdict
import difflib
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def main():
    panels = json.loads((HERE / 'RESOLVED_PANELS.json').read_text())['panels']
    unresolved = list(dict.fromkeys(s for panel in panels for s in panel['unresolved_symbols']))
    native = {g['gene_symbol'] for g in json.loads((REPO / 'site/data/MOUSE_CATALOG.json').read_text())['genes']}
    by_text, books = defaultdict(list), []
    for path in sorted((REPO / 'analysis/data/marker_workbooks').glob('*.xlsx')):
        text_count = 0
        relative = str(path.relative_to(REPO))
        with ZipFile(path) as archive:
            shared = []
            if 'xl/sharedStrings.xml' in archive.namelist():
                for item in ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('m:si', NS):
                    shared.append(''.join(t.text or '' for t in item.findall('.//m:t', NS)))
            workbook = ET.fromstring(archive.read('xl/workbook.xml'))
            relationships = {r.attrib['Id']: r.attrib['Target'] for r in ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
            for sheet in workbook.findall('m:sheets/m:sheet', NS):
                key = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
                target = relationships[key]
                target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
                for cell in ET.fromstring(archive.read(target)).findall('.//m:sheetData/m:row/m:c', NS):
                    kind, value = cell.attrib.get('t'), cell.find('m:v', NS)
                    if kind == 's' and value is not None:
                        text = shared[int(value.text)]
                    elif kind == 'inlineStr':
                        text = ''.join(t.text or '' for t in cell.findall('.//m:t', NS))
                    elif value is not None:
                        text = value.text
                    else:
                        continue
                    if not text or text.replace('.', '', 1).replace('-', '', 1).isdigit():
                        continue
                    text_count += 1
                    symbol = text.strip()
                    if (symbol in native or symbol in unresolved) and len(by_text[symbol]) < 4:
                        by_text[symbol].append(dict(workbook=relative, sheet=sheet.attrib['name'], cell=cell.attrib['r']))
        books.append(dict(path=relative, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), text_cell_count=text_count))
    candidates = [s for s in by_text if s in native]
    review = []
    for symbol in unresolved:
        matches = difflib.get_close_matches(symbol, candidates, n=4, cutoff=.8)
        review.append(dict(requested_symbol=symbol, exact_literal_workbook_examples=by_text.get(symbol, [])[:4],
                           candidate_symbols_not_applied=[dict(candidate=s,
                               reason='Similar spelling exists in provided workbook and native annotation; not an approved identity mapping.',
                               source_examples=by_text[s][:3]) for s in matches]))
    result = dict(schema='satb2-marker-workbook-spelling-review-1',
                  policy='Read-only inspection of provided workbook copies. Workbook text is source data, not instructions. Exact matches and spelling similarities are evidence for user review only; no additional mapping applied.',
                  source_books=books, review=review)
    (HERE / 'WORKBOOK_SPELLING_REVIEW.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f'Reviewed {len(books)} unchanged workbooks; {len(unresolved)} unique unresolved literals; no mappings changed.')


if __name__ == '__main__':
    main()
