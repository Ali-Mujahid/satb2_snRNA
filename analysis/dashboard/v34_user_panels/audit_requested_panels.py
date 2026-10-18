#!/usr/bin/env python3
"""Audit literal user marker requests against the unchanged native gene axis.

Exact, unique source-symbol matches only are selected. Alias suggestions are
recorded separately and never used as automatic replacements. Requested order,
spelling and duplicate positions remain in the audit.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def main():
    request_path = HERE / 'REQUESTED_PANELS.json'
    source_path = REPO / 'site/data/MOUSE_CATALOG.json'
    baseline_path = HERE / 'BASE_MARKER_REGISTRY.json'
    if not baseline_path.exists():
        save(baseline_path, read(REPO / 'site/data/manifest.json')['genes'])
    requests = read(request_path)
    source = read(source_path)['genes']
    baseline = read(baseline_path)
    by_symbol, by_lower = defaultdict(list), defaultdict(list)
    for gene in source:
        by_symbol[gene['gene_symbol']].append(gene)
        by_lower[gene['gene_symbol'].casefold()].append(gene)
    old_ids = {gene['gene_id'] for gene in baseline}
    new_by_id, all_selected_by_id, panels, rows = {}, {}, [], []
    for panel in requests['panels']:
        literal = [text.strip() for text in panel['requested_text'].split(',') if text.strip()]
        occurrences = Counter(literal)
        first_position, selected, unresolved, ambiguous, repeated = {}, [], [], [], []
        for position, symbol in enumerate(literal, 1):
            first = first_position.setdefault(symbol, position)
            matches = by_symbol.get(symbol, [])
            status = ('EXACT_UNIQUE' if len(matches) == 1 else
                      'AMBIGUOUS_SOURCE_SYMBOL' if matches else 'ABSENT_LITERAL_SOURCE_SYMBOL')
            item = dict(position=position, requested_symbol=symbol,
                        status=status, first_position=first,
                        repeat_of_position=first if position != first else None,
                        source_genes=matches,
                        case_insensitive_candidates=by_lower.get(symbol.casefold(), []) if not matches else [])
            rows.append(dict(panel_id=panel['id'], **item))
            if position != first:
                continue
            if len(matches) == 1:
                gene = matches[0]
                selected.append(gene)
                all_selected_by_id.setdefault(gene['gene_id'], gene)
                if gene['gene_id'] not in old_ids:
                    new_by_id.setdefault(gene['gene_id'], gene)
            elif matches:
                ambiguous.append(item)
            else:
                unresolved.append(symbol)
            if occurrences[symbol] > 1:
                repeated.append(dict(requested_symbol=symbol, positions=[i + 1 for i, g in enumerate(literal) if g == symbol]))
        panels.append(dict(id=panel['id'], label=panel['label'], requested_entry_count=len(literal),
                           unique_literal_count=len(first_position), selected_gene_count=len(selected),
                           selected_genes=selected, selected_gene_ids=[g['gene_id'] for g in selected],
                           source_resolved_symbols_in_requested_order=[g['gene_symbol'] for g in selected],
                           new_gene_count=sum(g['gene_id'] not in old_ids for g in selected),
                           absent_literal_count=len(unresolved), absent_literal_symbols=unresolved,
                           ambiguous_symbols=ambiguous, repeated_requested_symbols=repeated))
    result = dict(schema='satb2-user-marker-panel-literal-audit-1',
                  policy='Select exact unique source symbols only, preserving first requested occurrence. No aliases, spelling corrections, case-insensitive matches or inferred zeros applied.',
                  requested_panels_sha256=sha(request_path), source_catalog_sha256=sha(source_path),
                  baseline_marker_registry_sha256=sha(baseline_path),
                  source_gene_count=len(source), baseline_marker_gene_count=len(baseline),
                  new_exact_gene_count=len(new_by_id), selected_unique_gene_count=len(all_selected_by_id),
                  panels=panels, entries=rows)
    save(HERE / 'LITERAL_GENE_AUDIT.json', result)
    save(HERE / 'EXACT_NEW_GENES.json', list(new_by_id.values()))
    with (HERE / 'LITERAL_GENE_AUDIT.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['panel_id', 'position', 'requested_symbol', 'status', 'first_position', 'repeat_of_position', 'source_gene_ids', 'source_gene_symbols', 'case_insensitive_candidates'])
        writer.writeheader()
        for row in rows:
            writer.writerow({**{k: row[k] for k in writer.fieldnames[:6]},
                             'source_gene_ids': ';'.join(g['gene_id'] for g in row['source_genes']),
                             'source_gene_symbols': ';'.join(g['gene_symbol'] for g in row['source_genes']),
                             'case_insensitive_candidates': ';'.join(g['gene_symbol'] for g in row['case_insensitive_candidates'])})
    lines = ['# Exact marker-panel source audit', '',
             'Literal requested symbols were matched to the unchanged full native gene axis. Missing symbols are unavailable in that annotation, not measured zeros. No spelling correction or alias was applied. Repeated symbols select the first occurrence once and retain every original position in the JSON/CSV audit.', '',
             '| Panel | Requested entries | Unique literals | Exact unique genes | New vs 889 | Unresolved literals |',
             '|---|---:|---:|---:|---:|---:|']
    for panel in panels:
        lines.append(f"| {panel['label']} | {panel['requested_entry_count']} | {panel['unique_literal_count']} | {panel['selected_gene_count']} | {panel['new_gene_count']} | {panel['absent_literal_count']} |")
    lines.extend(['', f'Across all panels, {len(new_by_id)} new exact genes and {len(all_selected_by_id)} unique resolved genes were found.', ''])
    for panel in panels:
        lines += [f"## {panel['label']}", '', 'Unresolved literal symbols: ' + ', '.join(panel['absent_literal_symbols']) + '.', '']
        if panel['repeated_requested_symbols']:
            lines += ['Repeated entries: ' + '; '.join(f"{x['requested_symbol']} at {', '.join(map(str, x['positions']))}" for x in panel['repeated_requested_symbols']) + '.', '']
    (HERE / 'LITERAL_GENE_AUDIT.md').write_text('\n'.join(lines))
    print(json.dumps({k: result[k] for k in ['source_gene_count', 'baseline_marker_gene_count', 'new_exact_gene_count', 'selected_unique_gene_count']}))
    for panel in panels:
        print(panel['id'], panel['selected_gene_count'], 'exact;', panel['absent_literal_count'], 'unresolved')


if __name__ == '__main__':
    main()
