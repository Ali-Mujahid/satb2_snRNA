from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
import imageio.v2 as imageio,numpy as np
P=Path('/Users/mujahid357/Documents/codex_new/SATB2_SCENE_SEURAT_20261111'); W,H=1600,900
font='/System/Library/Fonts/Supplemental/Arial.ttf'; bold='/System/Library/Fonts/Supplemental/Arial Bold.ttf'; title=ImageFont.truetype(bold,32); small=ImageFont.truetype(font,20)
def slide(file,heading,sub):
 bg=Image.new('RGB',(W,H),'white'); im=Image.open(P/'figures'/file).convert('RGB'); im.thumbnail((1480,730),Image.Resampling.LANCZOS); bg.paste(im,((W-im.width)//2,100+(730-im.height)//2)); d=ImageDraw.Draw(bg);d.text((45,25),heading,font=title,fill='#163a4d');d.text((46,860),sub,font=small,fill='#4d6470');return bg
slides=[slide('scene_joint_clusters.png','SCENE · CA1 / ProS Seurat clusters','Actual SCENE joint UMAP; colors indicate saved Seurat clusters.'),slide('scene_joint_layers.png','SCENE · superficial and deep layer-associated programs','Blue: superficial-associated · orange: deep-associated · gray: mixed / weak.'),slide('cluster_grouped_30_dotplot.png','Seurat · 30 marker candidates grouped by cluster','Dot size: detected fraction · color: relative expression; pooled WT and cKO.')]
writer=imageio.get_writer(str(P/'SCENE_CA1_Seurat.mp4'),fps=24,codec='libx264',quality=8,macro_block_size=16)
# Fade between actual static model projections; no inferred temporal cell movement.
for i,img in enumerate(slides):
 for frame in range(120):
  if i and frame<24: frame_img=Image.blend(slides[i-1],img,frame/24)
  else: frame_img=img
  writer.append_data(np.asarray(frame_img))
writer.close();print('MP4 complete',flush=True)
