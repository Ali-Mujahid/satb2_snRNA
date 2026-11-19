/* Export the current C6 viewport and its category key using the positions drawn by the atlas. */
window.SATB2_C6_FIGURE = async function saveC6Figure(scene, format, download) {
  const ns='http://www.w3.org/2000/svg';
  function element(name, attributes, parent, content) {
    const el=document.createElementNS(ns,name);
    for(const [key,value] of Object.entries(attributes))el.setAttribute(key,String(value));
    if(content!==undefined)el.textContent=content;
    if(parent)parent.append(el);
    return el;
  }
  const width=Math.max(620,scene.width+32);
  const plotTop=54, legendTop=plotTop+scene.height+25, line=19;
  const height=legendTop+scene.legend.length*line+42;
  const svg=element('svg',{width,height,viewBox:`0 0 ${width} ${height}`,'font-family':'Arial, sans-serif'});
  element('title',{},svg,'Inhibitory neurons with Seurat C6 subclusters');
  element('rect',{width,height,fill:'white'},svg);
  element('text',{x:16,y:24,'font-size':18,'font-weight':700},svg,'Inhibitory neurons with Seurat C6 subclusters');
  element('text',{x:16,y:44,'font-size':13},svg,scene.heading);
  const plot=element('svg',{x:16,y:plotTop,width:scene.width,height:scene.height,overflow:'hidden'},svg);
  if(scene.guide.length)element('polyline',{points:scene.guide.map(p=>p.x+','+p.y).join(' '),stroke:'#e4edf1','stroke-width':18,fill:'none','stroke-linecap':'round'},plot);
  for(const point of scene.points)
    element('rect',{x:point.x,y:point.y,width:scene.size,height:scene.size,fill:point.color,'fill-opacity':scene.opacity},plot);
  for(const t of scene.labels)
    element('text',{x:t.x,y:t.y,'text-anchor':'middle','font-size':scene.font,'font-weight':700,fill:'#253e50',stroke:'white','stroke-width':3,'paint-order':'stroke'},plot,t.text);
  scene.legend.forEach((item,index)=>{
    const y=legendTop+index*line;
    element('rect',{x:18,y:y-10,width:10,height:10,fill:item.color},svg);
    element('text',{x:38,y,'font-size':13},svg,item.label);
  });
  element('text',{x:16,y:height-12,'font-size':12},svg,'C6 names describe provisional RNA marker programs. Positions are embedding coordinates.');
  const source=new XMLSerializer().serializeToString(svg);
  if(format==='svg'){download('Seurat_inhibitory_C6.svg',source,'image/svg+xml');return;}
  const url=URL.createObjectURL(new Blob([source],{type:'image/svg+xml'}));
  try {
    const image=new Image();image.src=url;await image.decode();
    const canvas=document.createElement('canvas');
    canvas.width=width*2;canvas.height=height*2;canvas.getContext('2d').drawImage(image,0,0,canvas.width,canvas.height);
    const png=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));
    if(!png)throw Error('PNG creation failed');
    download('Seurat_inhibitory_C6.png',png,'image/png');
  } finally {URL.revokeObjectURL(url);}
};
