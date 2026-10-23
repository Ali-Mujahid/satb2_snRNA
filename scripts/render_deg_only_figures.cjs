#!/usr/bin/env node
'use strict';
/* Re-render the corrected saved native DEG-only records. No biological model is fitted. */
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),zlib=require('node:zlib'),crypto=require('node:crypto'),cp=require('node:child_process');
const repo=path.resolve(__dirname,'..');
const runtime=process.env.CODEX_RUNTIME_DEPENDENCIES || '/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies';
const sharp=require(require.resolve('sharp',{paths:[path.join(runtime,'node/node_modules')]}));
const python=process.env.CCC_FIGURE_PYTHON || path.join(runtime,'python/bin/python3');
const out=path.resolve(process.argv[2] || path.join(repo,'../SATB2_CCC_20261022/deliverables/DEG_ONLY_20261023/figures'));
const input=path.join(repo,'site/communication/data/deg_only_data.js');
const sources=['site/communication/data/deg_only_data.js','site/communication/model.js','site/communication/charts.js','site/communication/networks.js','scripts/render_deg_only_figures.cjs'];
const hash=buf=>crypto.createHash('sha256').update(buf).digest('hex');
const context={window:{}};
for(const file of ['site/communication/data/deg_only_data.js','site/communication/model.js','site/communication/charts.js','site/communication/networks.js'])vm.runInNewContext(fs.readFileSync(path.join(repo,file),'utf8'),context,{filename:file});
const M=context.window.CCCModel,P=context.window.CCCPlots,N=context.window.CCCNetworks;
const packed=JSON.parse(zlib.gunzipSync(Buffer.from(context.window.CCC_DATA_GZIP,'base64'))),data=M.decode(packed),dataset=data.datasets['IN_clusters|pooled'];
if(!dataset)throw Error('Corrected pooled IN-clusters dataset is unavailable.');
const options={maxRank:1,tier:'DEG_only',bothDEG:true,ligandDirection:'any',receptorDirection:'any',directionPair:'any',search:'',sources:[],targets:[],condition:'pooled',showDEGDirections:true,degUpColor:'#c63d4f',degDownColor:'#397eb5',degMixedColor:'#8b6bb1',dotMetric:'magnitude_rank',heatMetric:'count',palette:'viridis',scaleMode:'auto',reverseScale:false,opacity:.8,fontSize:12,topN:120,pairLimit:80,order:'rank',dotRadius:9,edgeWidth:5,arrowSize:10,curvature:.35,arcWidth:17,arcGap:3,ribbonScale:1,nodeScale:1,edgeColorMode:'source',arcColorMode:'node',geneEdges:30,nodeColorMode:'palette'};
const directions=[['CA1_to_IN','CA1_to_inhibitory'],['IN_to_CA1','inhibitory_to_CA1'],['both','bidirectional']];
const views=[['dotplot',P.dotplot],['heatmap',P.heatmap],['circos',N.circos],['circle',N.circle],['geneCircle',N.geneCircle]];
const original=JSON.stringify(dataset.rows);
function describeDirections(rows){const out={};for(const r of rows)out[r.deg_direction_pair]=(out[r.deg_direction_pair]||0)+1;return out;}
function addProvenance(svg,description){
  const safe=x=>String(x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'}[c]));
  return svg.replace('</svg>','<metadata id="strict-deg-provenance">'+safe(JSON.stringify(description))+'</metadata></svg>');
}
async function main(){
  fs.mkdirSync(out,{recursive:true});
  const expected=15,outputFiles=[];
  for(const [direction,label] of directions){
    const state={...options,direction},rows=M.filter(dataset.rows,state),eligiblePairs=M.eligible(dataset.eligible_pairs,state);
    if(rows.some(r=>!M.strictDEG(r)))throw Error('Nonmatching DEG endpoint escaped the corrected biological selection.');
    const nodes=data.nodes.map(n=>({...n,label:n.id==='CA1_Deep'?'CA1 Deep':n.id==='CA1_Superficial'?'CA1 Sup.':n.id.replace(/^IN_/,'')}));
    const ctx={nodes,eligiblePairs,condition:'pooled',scope:'IN_clusters',direction,options:state};
    const directionCounts=describeDirections(rows),sourceFile='strict_DEG_only_'+label+'_source_rows.tsv';
    fs.writeFileSync(path.join(out,sourceFile),M.tsv(rows,packed.columns));
    for(const [view,render] of views){
      const stem='strict_DEG_only_pooled_'+label+'_'+view,provenance={analysis:'Strict matching-endpoint DEG-only projection of completed native LIANA results; no refits',populations:data.nodes.map(n=>n.id),direction,condition:'pooled',records:rows.length,deg_rule:'Each ligand subunit is DESeq2 padj < 0.05 in its sender population AND each receptor subunit is DESeq2 padj < 0.05 in its receiver population; unambiguous symbol mapping.',genotype_effect:'DESeq2 log2(cKO / WT) at endpoint gene level only; native lr_logfc remains group-versus-rest specificity; no interaction genotype test.',direction_counts:directionCounts,rank_filter:'magnitude_rank <= 1',source_data_file:sourceFile,source_sha256:hash(fs.readFileSync(input)),display_options:state};
      const svg=addProvenance(render(rows,ctx),provenance);
      if(/(?:NaN|Infinity|undefined)/.test(svg))throw Error('Invalid finite SVG geometry: '+stem);
      if(!svg.includes('data-deg-summary="true"'))throw Error('Signed endpoint evidence absent: '+stem);
      fs.writeFileSync(path.join(out,stem+'.svg'),svg);
      const png=await sharp(Buffer.from(svg),{density:192}).png().toBuffer();
      fs.writeFileSync(path.join(out,stem+'.png'),png);
      const dims=await sharp(png).metadata();
      outputFiles.push({stem,view,direction,direction_label:label,native_rows:rows.length,molecular_pairs:new Set(rows.map(M.lrKey)).size,deg_direction_counts:directionCounts,source_data_file:sourceFile,width_px:dims.width,height_px:dims.height,logical_width:dims.width/2,logical_height:dims.height/2});
    }
  }
  if(outputFiles.length!==expected)throw Error('Expected exactly 15 complete scientific views.');
  if(JSON.stringify(dataset.rows)!==original)throw Error('Saved native rows were mutated.');
  fs.writeFileSync(path.join(out,'figure_manifest.json'),JSON.stringify({status:'RENDERED_SVG_PNG_PENDING_PDF',expected_views:expected,options,source_hashes:Object.fromEntries(sources.map(p=>[p,hash(fs.readFileSync(path.join(repo,p)))])),figures:outputFiles},null,2)+'\n');
  const py=String.raw`
import json,hashlib,subprocess,concurrent.futures
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from pypdf import PdfReader
from PIL import Image,ImageDraw,ImageFont
out=Path(__import__('sys').argv[1])
manifest=json.loads((out/'figure_manifest.json').read_text())
for fig in manifest['figures']:
    stem=fig['stem']; png=out/(stem+'.png'); pdf=out/(stem+'.pdf')
    w,h=fig['logical_width']*.75,fig['logical_height']*.75
    c=canvas.Canvas(str(pdf),pagesize=(w,h),pageCompression=1,invariant=1)
    c.setTitle('Strict DEG-only CA1-inhibitory communication: '+fig['direction_label']+' '+fig['view'])
    c.setAuthor('SATB2 snRNA atlas')
    c.setSubject('Corrected matching-endpoint DESeq2 padj < 0.05; pooled native LIANA candidates; endpoint gene changes only, no interaction-level genotype test')
    c.drawImage(ImageReader(str(png)),0,0,width=w,height=h)
    c.showPage();c.save()
    p=PdfReader(pdf)
    assert len(p.pages)==1 and abs(float(p.pages[0].mediabox.width)-w)<.01 and abs(float(p.pages[0].mediabox.height)-h)<.01
    for ext in ['svg','png','pdf']:
        target=out/(stem+'.'+ext)
        fig[ext]={'file':target.name,'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
font_path='/System/Library/Fonts/Supplemental/Arial.ttf'
try: font=ImageFont.truetype(font_path,18)
except OSError: font=ImageFont.load_default()
cellw,cellh=800,750
sheet=Image.new('RGB',(cellw*3,cellh*5),'#dfe7ec')
draw=ImageDraw.Draw(sheet)
for i,fig in enumerate(manifest['figures']):
    col=i//5;row=i%5;im=Image.open(out/(fig['stem']+'.png')).convert('RGB')
    im.thumbnail((cellw-26,cellh-54))
    x=col*cellw+(cellw-im.width)//2;y=row*cellh+40+(cellh-54-im.height)//2
    sheet.paste(im,(x,y))
    draw.text((col*cellw+13,row*cellh+10),fig['direction_label']+' / '+fig['view']+' / '+str(fig['native_rows'])+' rows',font=font,fill='#233f50')
sheet.save(out/'CONTACT_SHEET.png')
tmp=out/'tmp'/'pdfs';tmp.mkdir(parents=True,exist_ok=True)
def verify_pdf(fig):
    stem=fig['stem'];target=tmp/stem
    subprocess.run(['pdftoppm','-png','-singlefile','-scale-to-x','1200','-scale-to-y','-1',str(out/(stem+'.pdf')),str(target)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    image=Image.open(target.with_suffix('.png')).convert('RGB')
    assert image.width==1200 and abs(image.height/image.width-fig['height_px']/fig['width_px'])<.002
    extrema=image.getextrema();assert all(lo<245 and hi==255 for lo,hi in extrema)
    return {'file':stem+'.pdf','pages':1,'roundtrip_dimensions':[image.width,image.height],'roundtrip_sha256':hashlib.sha256(target.with_suffix('.png').read_bytes()).hexdigest(),'status':'PASS'}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    manifest['pdf_roundtrip_checks']=list(pool.map(verify_pdf,manifest['figures']))
roundtrip=Image.new('RGB',(cellw*3,cellh*5),'#dfe7ec');draw=ImageDraw.Draw(roundtrip)
for i,fig in enumerate(manifest['figures']):
    col=i//5;row=i%5;im=Image.open(tmp/(fig['stem']+'.png')).convert('RGB');im.thumbnail((cellw-26,cellh-54))
    roundtrip.paste(im,(col*cellw+(cellw-im.width)//2,row*cellh+40+(cellh-54-im.height)//2))
    draw.text((col*cellw+13,row*cellh+10),fig['direction_label']+' / '+fig['view']+' / PDF',font=font,fill='#233f50')
roundtrip.save(out/'PDF_ROUNDTRIP_CONTACT_SHEET.png')
manifest['status']='PASS';manifest['pdf_count']=len(manifest['figures']);manifest['pdf_encoding']='Rendered at 192 dpi with Unicode-aware SVG fonts, then embedded losslessly into one-page PDFs. Matching vector SVGs are also provided.'
manifest['contact_sheet']='CONTACT_SHEET.png'
(out/'figure_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
`;
  cp.execFileSync(python,['-c',py,out],{stdio:'inherit'});
  fs.writeFileSync(path.join(out,'README.md'),'# Corrected strict DEG-only communication figures\n\n15 saved views: five chart types for CA1 to inhibitory, inhibitory to CA1, and the bidirectional union. Pooled native scores are reused without model refits.\n\nOnly five requested populations are active: CA1 Deep, CA1 Superficial, inhibitory C1, C6, and C8. Every ligand/receptor complex subunit must have unambiguous matching-population DESeq2 adjusted P < 0.05 on both endpoints. No rank prioritization is applied beyond native magnitude_rank <= 1.\n\nPooled view: 16 population-specific LR records; CA1 to inhibitory: 6; inhibitory to CA1: 10. C8 has no retained interactions and remains visible as a zero-edge selected population in the population-level plots. This does not demonstrate absent biological communication.\n\nUp/down refers to endpoint DESeq2 log2(cKO / WT). Dot left/right rims indicate ligand/receptor direction; gene-node outlines and arrow suffixes indicate endpoint gene direction. Interior dot color still shows native LIANA magnitude rank. Network arrows indicate listed ligand sender to receptor receiver, not activation or inhibition. Native lr_logfc is group-versus-rest specificity and is not a genotype fold change.\n\nThe gene-circle figures show all qualifying records (display limit 30 exceeds the 16-record maximum). All five selected groups are retained in heatmap, population circle, and circos views. SVG tooltips and source TSVs retain the signed effect evidence.\n\nPDFs embed the 192-dpi PNG losslessly to preserve Unicode direction arrows; SVGs retain vector lines and text. No browser was used for these renders. Numerical source rows, complete file hashes, plot defaults, record counts, and input checksums are included in figure_manifest.json.\n');
  const final=JSON.parse(fs.readFileSync(path.join(out,'figure_manifest.json'),'utf8'));
  console.log(JSON.stringify({status:final.status,output:out,pdf_count:final.pdf_count,views:final.figures.map(f=>({view:f.view,direction:f.direction,rows:f.native_rows})),manifest_sha256:hash(fs.readFileSync(path.join(out,'figure_manifest.json')))},null,2));
}
main().catch(error=>{console.error(error.stack || error);process.exitCode=1;});
