#!/usr/bin/env python3
"""Resolve requested panels with exact source IDs and disclosed official aliases.

Four verified official legacy synonyms are explicitly included. Every remaining
unavailable literal retains its status; no fuzzy spelling correction is used.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def read(path):
    return json.loads(path.read_text())


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')


def main():
    request = read(HERE / 'REQUESTED_PANELS.json')
    aliases = read(HERE / 'VERIFIED_ALIAS_PROPOSALS.json')['proposals']
    native = read(REPO / 'site/data/MOUSE_CATALOG.json')['genes']
    baseline = read(HERE / 'BASE_MARKER_REGISTRY.json')
    by_symbol = {}
    for gene in native:
        by_symbol.setdefault(gene['gene_symbol'], []).append(gene)
    by_alias = {x['requested_symbol']: x for x in aliases}
    old_ids = {x['gene_id'] for x in baseline}
    new, panels = {}, []
    for panel in request['panels']:
        symbols = [s.strip() for s in panel['requested_text'].split(',') if s.strip()]
        entries, seen_symbols, seen_ids, gene_ids = [], {}, set(), []
        for position, symbol in enumerate(symbols, 1):
            entry = dict(position=position, requested_symbol=symbol, gene_id=None, source_gene_symbol=None)
            matches = by_symbol.get(symbol, [])
            if len(matches) == 1:
                gene = matches[0]
                entry.update(status='exact', gene_id=gene['gene_id'], source_gene_symbol=gene['gene_symbol'])
            elif symbol in by_alias:
                alias = by_alias[symbol]
                gene = alias['unique_native_gene']
                assert gene in by_symbol[gene['gene_symbol']], 'Alias target must resolve on native axis'
                entry.update(status='verified_alias', gene_id=gene['gene_id'], source_gene_symbol=gene['gene_symbol'],
                             alias_evidence_url=alias['source_url'],
                             alias_additional_evidence_url=alias['additional_source_url'],
                             alias_status='Official synonym verified; requested label retained and native stable ID used explicitly.')
            else:
                entry.update(status='unresolved', reason='No exact source symbol or verified official synonym resolved. No inferred count or zero value.')
            if symbol in seen_symbols:
                entry.update(resolution_status=entry['status'], status='duplicate', first_requested_position=seen_symbols[symbol])
            else:
                seen_symbols[symbol] = position
            if entry['gene_id']:
                gid = entry['gene_id']
                gene = dict(gene_id=gid, gene_symbol=entry['source_gene_symbol'])
                if gid not in seen_ids:
                    seen_ids.add(gid)
                    gene_ids.append(gid)
                else:
                    entry['duplicate_selected_gene_id'] = True
                if gid not in old_ids:
                    new.setdefault(gid, gene)
            entries.append(entry)
        panels.append(dict(id=panel['id'], label=panel['label'], requested_symbols=symbols, entries=entries,
                           gene_ids=gene_ids, requested_entry_count=len(symbols), selected_gene_count=len(gene_ids),
                           unresolved_symbols=list(dict.fromkeys(x['requested_symbol'] for x in entries if x['status'] == 'unresolved' or x.get('resolution_status') == 'unresolved')),
                           verified_aliases=[dict(requested_symbol=s, source_gene_symbol=by_alias[s]['unique_native_gene']['gene_symbol'],
                                                  gene_id=by_alias[s]['unique_native_gene']['gene_id'], evidence_url=by_alias[s]['source_url'])
                                             for s in dict.fromkeys(symbols) if s in by_alias]))
    result = dict(schema='satb2-user-marker-panel-resolved-1',
                  policy='Exact unique source symbols plus four explicitly disclosed official legacy synonyms. Requested order retained. Repeats select the first stable gene ID once. Unresolved entries remain unavailable, never zero. No spelling or fuzzy corrections.',
                  request_sha256=hashlib.sha256((HERE / 'REQUESTED_PANELS.json').read_bytes()).hexdigest(), panels=panels)
    assert len(new) == 166
    assert [p['selected_gene_count'] for p in panels] == [99, 118, 97]
    assert [len(p['unresolved_symbols']) for p in panels] == [1, 34, 3]
    save(HERE / 'RESOLVED_PANELS.json', result)
    save(HERE / 'EXPORT_NEW_GENES.json', list(new.values()))
    print('PASS: 166 source-ID export additions; panel selections 99 / 118 / 97; unresolved 1 / 34 / 3')


if __name__ == '__main__':
    main()
