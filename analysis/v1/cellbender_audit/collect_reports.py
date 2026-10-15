from pathlib import Path
import csv, hashlib, json, re, subprocess
from pypdf import PdfReader
from PIL import Image, ImageDraw

ROOT=Path('/Users/mujahid357/Desktop/SATB2_LiSC_archive_transfer_20261012.BmzvXa/cellbender_out')
OUT=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender')
RENDER=OUT/'pdf_render'
RENDER.mkdir(parents=True,exist_ok=True)
records=[]
inventory=[]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
for directory in sorted(x for x in ROOT.iterdir() if x.is_dir()):
    sample=directory.name
    metrics=dict(csv.reader((directory/f'{sample}_metrics.csv').open()))
    logpath=directory/f'{sample}.log'
    log=logpath.read_text(errors='replace')
    record={'sample':sample,**{k:float(v) for k,v in metrics.items()},'completed':'Completed remove-background.' in log,'checkpoint_loaded':'Checkpoint loaded successfully.' in log,'checkpoint_mismatch_rejected':'Workflow hash does not match that of checkpoint.' in log,'no_saved_checkpoint':'No saved checkpoint.' in log,'html_report_error':'Unable to create report.' in log,'command':log.splitlines()[1],'source_log_path':str(logpath)}
    records.append(record)
    for suffix,kind in [('.h5','full_corrected'),('_filtered.h5','filtered_corrected'),('.log','training_log'),('_metrics.csv','metrics'),('.pdf','summary_pdf'),('_cell_barcodes.csv','called_barcodes')]:
        path=directory/f'{sample}{suffix}'
        inventory.append({'sample':sample,'role':kind,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
    pdf=directory/f'{sample}.pdf'
    reader=PdfReader(pdf)
    (RENDER/f'{sample}_text.txt').write_text('\n'.join(p.extract_text() for p in reader.pages))
    subprocess.run(['pdftoppm','-f','1','-singlefile','-r','90','-png',str(pdf),str(RENDER/sample)],check=True,capture_output=True)
with (OUT/'source_artifact_inventory.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(inventory[0]),delimiter='\t');w.writeheader();w.writerows(inventory)
(OUT/'recorded_metrics_and_log_flags.json').write_text(json.dumps(records,indent=2))
with (OUT/'recorded_metrics_and_log_flags.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
for offset in range(0,len(records),4):
    group=records[offset:offset+4]
    canv=Image.new('RGB',(1800,1420),'white')
    draw=ImageDraw.Draw(canv)
    for pos,record in enumerate(group):
        im=Image.open(RENDER/f'{record["sample"]}.png').convert('RGB');im.thumbnail((880,660))
        x=(pos%2)*900;y=(pos//2)*710
        draw.text((x+15,y+10),record['sample'],fill='black')
        canv.paste(im,(x+10,y+35))
    canv.save(RENDER/f'contact_{offset//4+1}.png')
print('collected',len(records),'samples,',len(inventory),'hashed artifacts')
