#!/usr/bin/env python3
"""Port marker-order UI to a derivative runtime while preserving its own shard loader.

The source is the hosted frontend in site/. Only named display methods, the marker
panel HTML, CSS and pure display helpers are copied. No scientific payload is read
or modified. The caller owns the derivative folder and any archive/manifest updates.
"""
from pathlib import Path
import argparse, hashlib, json, re

FUNCTIONS = ['defaults','markerRange','renderDots','render','exportDotValues','exportState','importState','exportPDF','bindControls','init']
NEW_API = 'markerOrderResult,markerGeneIds,orderedMarkerGroups,markerColumnGroups,selectMarkerPreset,applyMarkerOrder,moveMarker,markerOrderNote'

def named_function(source,name):
    match=re.search(r'^  (?:async )?function '+re.escape(name)+r'\(',source,re.M)
    if not match: raise ValueError('Missing function '+name)
    end=re.search(r'^  (?:async )?function |^  /\* MARKER_ORDERING_UI_START|^  const api=',source[match.end():],re.M)
    if not end: raise ValueError('Missing boundary after '+name)
    return match.start(),match.end()+end.start(),source[match.start():match.end()+end.start()]

def loader(source):
    return source[source.index('  async function loadShard(id) {'):source.index('  function descriptor(tool)')]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1]/'site');parser.add_argument('--target',type=Path,required=True);args=parser.parse_args()
    source=args.source.resolve();target=args.target.resolve()
    assert source!=target, 'Use an owned derivative folder as the target'
    before=(target/'app.js').read_text();s=before;reference=(source/'app.js').read_text();before_loader=loader(s)
    for name in FUNCTIONS:
        a,b,_=named_function(s,name);replacement=named_function(reference,name)[2];s=s[:a]+replacement+s[b:]
    start=reference.index('  /* MARKER_ORDERING_UI_START');stop=reference.index('  /* MARKER_ORDERING_UI_END */',start)+len('  /* MARKER_ORDERING_UI_END */')
    block=reference[start:stop]+'\n'
    if '  /* MARKER_ORDERING_UI_START' in s:
        a=s.index('  /* MARKER_ORDERING_UI_START');b=s.index('  /* MARKER_ORDERING_UI_END */',a)+len('  /* MARKER_ORDERING_UI_END */')
        s=s[:a]+block.rstrip('\n')+s[b:]
    else:s=s.replace('  function numericRange(tool) {',block+'  function numericRange(tool) {')
    if NEW_API not in s:s=s.replace('umapIndices,pageGenes,suggestionFor,get state()', 'umapIndices,pageGenes,suggestionFor,'+NEW_API+',get state()')
    assert NEW_API in s and loader(s)==before_loader, 'Port must retain target classic/native data loader exactly'
    (target/'app.js').write_text(s)
    old=(target/'index.html').read_text();host=(source/'index.html').read_text()
    section=lambda text: re.search(r'<section id="tab-markers".*?(?=<section id="tab-de")',text,re.S)
    m=section(old);new=section(host);assert m and new
    index=(old[:m.start()]+new.group()+old[m.end():]).replace('Search 878 available genes','Search all packaged genes').replace('including all 878 panel genes','including all packaged panel genes')
    (target/'index.html').write_text(index)
    for name in ['styles.css','v2_features.js']:(target/name).write_bytes((source/name).read_bytes())
    print(json.dumps({'status':'PASS','target_classic_loader_preserved':True,'changed_files':['app.js','index.html','styles.css','v2_features.js'],'source_frontend_sha256':hashlib.sha256(reference.encode()).hexdigest(),'target_frontend_sha256':hashlib.sha256(s.encode()).hexdigest(),'scientific_payloads_modified':False}))

if __name__=='__main__':main()
