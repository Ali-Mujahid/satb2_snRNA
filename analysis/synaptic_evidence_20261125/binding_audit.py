"""Attach published binding and compartment evidence to the frozen RNA candidates."""
import csv
import json
from collections import Counter
from pathlib import Path


def write_tsv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def annotate_pairs(pairs, output: Path):
    audit = json.loads((output / 'binding_audit.json').read_text())
    definitions = {frozenset(p['genes']): p for p in audit['pairs']}
    assert len(definitions) == 10 and len(pairs) == 14
    annotated = []
    for original in pairs:
        p = dict(original)
        a = definitions[frozenset((p['ligand'], p['receptor']))]
        ca_source = p['source'].startswith('CA1_')
        ca_gene = p['ligand'] if ca_source else p['receptor']
        in_gene = p['receptor'] if ca_source else p['ligand']
        p.update(
            CA1_sublayer=p['source'] if ca_source else p['target'], CA1_gene=ca_gene,
            IN_cluster=p['target'] if ca_source else p['source'], IN_gene=in_gene,
            CA1_log2FC=p['ligand_log2FC'] if ca_source else p['receptor_log2FC'],
            IN_log2FC=p['receptor_log2FC'] if ca_source else p['ligand_log2FC'],
            CA1_FDR=p['ligand_FDR'] if ca_source else p['receptor_FDR'],
            IN_FDR=p['receptor_FDR'] if ca_source else p['ligand_FDR'],
            binding_status=a['status'], binding_verdict=a['title'],
            binding_assay=a['assay'], binding_context=a['experimental_context'],
            molecular_form=a['molecular_form'], interpretation=a['interpretation'],
            binding_sources=' | '.join(r['url'] for r in a['references']),
            CA1_compartment=a['locations'][ca_gene], IN_compartment=a['locations'][in_gene],
            contact_validation='Predicted P60 population relationship; cell-specific contact untested',
            corrected_source=p['source'], corrected_target=p['target'],
            corrected_source_gene=p['ligand'], corrected_target_gene=p['receptor'],
            molecular_mode='signal', synaptic_source='', synaptic_target='',
            synaptic_direction='Unresolved', synaptic_effect='Unresolved',
        )
        if a['status'] in ('unverified', 'indirect'):
            p['molecular_mode'] = 'unresolved'
        if frozenset(a['genes']) == frozenset(('Nrxn1', 'Nlgn1')):
            p.update(corrected_source='CA1_Superficial', corrected_target='IN_C1',
                     corrected_source_gene='Nrxn1', corrected_target_gene='Nlgn1',
                     molecular_mode='adhesion', synaptic_source='CA1_Superficial',
                     synaptic_target='IN_C1', synaptic_direction='CA1 superficial → IN C1 (predicted)',
                     synaptic_effect='Excitatory glutamatergic connection predicted by canonical protein localization')
        elif ca_gene == 'Ptprg':
            p['molecular_mode'] = 'adhesion'
        p['corrected_pair'] = p['corrected_source_gene'] + (' ↔ ' if p['molecular_mode'] == 'adhesion' else ' → ' if p['molecular_mode'] == 'signal' else ' ··· ') + p['corrected_target_gene']
        annotated.append(p)
    assert Counter(p['binding_status'] for p in annotated) == dict(direct=10, screen=1, indirect=1, unverified=2)
    assert sum(bool(p['synaptic_source']) for p in annotated) == 2
    keys = ['CA1_sublayer','CA1_gene','CA1_log2FC','CA1_FDR','IN_cluster','IN_gene','IN_log2FC','IN_FDR',
            'source','target','ligand','receptor','binding_status','binding_verdict','corrected_pair',
            'corrected_source','corrected_target','molecular_mode','CA1_compartment','IN_compartment',
            'synaptic_direction','synaptic_effect','binding_assay','binding_context','molecular_form',
            'interpretation','contact_validation','binding_sources']
    write_tsv(output / 'binding_audit_14_rows.tsv', [{k: p[k] for k in keys} for p in annotated])
    sources = [dict(pair=' / '.join(a['genes']), status=a['status'], citation=r['label'], url=r['url'])
               for a in audit['pairs'] for r in a['references']]
    write_tsv(output / 'binding_audit_sources.tsv', sources)
    return annotated, audit
