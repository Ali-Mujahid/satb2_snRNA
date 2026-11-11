import os,sys,hashlib,json,inspect,importlib.util
from pathlib import Path
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
import numpy as np,anndata as ad,pandas as pd,torch
from scene import SCENE
import scene.train as official_train
import scene.model as scene_model
source=Path(sys.argv[1]);inp=Path(sys.argv[2]);out=Path(sys.argv[3]);expected_input=sys.argv[4];expected_train=sys.argv[5]
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
assert digest(inp)==expected_input,'input hash mismatch'
assert digest(source)==expected_train,'frozen training source mismatch'
assert digest(Path(official_train.__file__))==expected_train,'installed training implementation differs from frozen source'
assert not out.exists(),'Refusing to overwrite output';out.mkdir(parents=True); snap=out/'checkpoints';snap.mkdir()
A=ad.read_h5ad(inp);assert A.shape==(27703,2518);assert A.obs['sample'].nunique()==12 and A.obs['condition'].value_counts().to_dict()=={'WT':14984,'cKO':12719};assert np.all(A.X.data>=0) and np.all(A.X.data==np.floor(A.X.data))
# Instrument an exact frozen copy of the SCENE 0.1.2 trainer, adding only coordinate snapshots.
src=source.read_text(); needle='    optimizer_updates = 0\n    for epoch in range(epochs):';assert src.count(needle)==1
src=src.replace(needle,'''    optimizer_updates = 0
    np.savez_compressed(_capture_dir / "epoch_0000.npz", cells=model.Z_cells.detach().cpu().numpy(), genes=model.Z_genes.detach().cpu().numpy())
    for epoch in range(epochs):''')
needle='        should_validate = (\n';assert src.count(needle)==1
src=src.replace(needle,'''        if epoch == 0 or (epoch + 1) % 5 == 0 or epoch + 1 == epochs:
            np.savez_compressed(_capture_dir / f"epoch_{epoch+1:04d}.npz", cells=model.Z_cells.detach().cpu().numpy(), genes=model.Z_genes.detach().cpu().numpy())

        should_validate = (
''')
ns={'__name__':'scene._capture_trainer','__package__':'scene','_capture_dir':snap};exec(compile(src,str(source),'exec'),ns);scene_model._train_prepared=ns['_train_prepared']
print('NATIVE 3D SCENE start:',A.shape,'mice',A.obs['sample'].nunique(),flush=True)
model=SCENE(latent_dim=3,likelihood='ZIP'); model.fit(A,layer=None,batch_keys=['sample'],epochs=1200,seed=42,validate=True,val_interval=20,print_every=20,cell_batch_size=512,loader_workers=0,device='cpu');model.save(out/'native_model');A.write_h5ad(out/'fitted.h5ad',compression='gzip');pd.DataFrame(model.training_history['train']).to_csv(out/'training.tsv',sep='\t',index=False);pd.DataFrame(model.training_history['validation']).to_csv(out/'validation.tsv',sep='\t',index=False)
assert len(list(snap.glob('epoch_*.npz')))==242
assert digest(inp)==expected_input and digest(source)==expected_train
(out/'COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','n_cells':A.n_obs,'n_genes':A.n_vars,'epochs':1200,'latent_dim':3,'capture_every':5,'checkpoints':242,'seed':42,'batch_key':'sample','likelihood':'ZIP','input_sha256':expected_input,'trainer_sha256':expected_train,'validation':'10% held-out positive edges plus matched negatives','interpretation':'SCENE native Euclidean 3D fit; no anatomical or genotype-DE inference'},indent=2))
print('FIT COMPLETE',flush=True)
