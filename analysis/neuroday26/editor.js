/* The artwork is the document. Controls and selection handles live outside the SVG. */
(() => {
"use strict";
const $ = id => document.getElementById(id);
const NS = "http://www.w3.org/2000/svg", XL = "http://www.w3.org/1999/xlink";
const LEAVES = "text,path,circle,rect,ellipse,line,polygon,polyline,image,use";
const OBJECTS = LEAVES + ",g,a,svg";
const GROUPS = "g,a,svg";
const assetBase = window.POSTER_ASSET_BASE || document.baseURI;
const storageKey = location.pathname.replace(/index\.html$/,"").replace(/\/$/,"");
const svgNode = (name, attrs={}) => {
  const node = document.createElementNS(NS,name);
  for (const [key,val] of Object.entries(attrs)) node.setAttribute(key,val);
  return node;
};
let poster, original, selected=null, mode="object", zoom=1, drag=null;
let undoStack=[], redoStack=[], cachedLayers=[], db=null, dirty=false, timer=0;
let saveNumber=0, saveChain=Promise.resolve(), sourceHash="", lastDraftID=null, storageBlocked=false;
let generation=0, selectionBox=null, imageTarget=null;
const short = s => s.replace(/\s+/g," ").trim().slice(0,92);
const round = v => Math.round(v*1000)/1000;
const finite = v => v!=="" && Number.isFinite(Number(v));
const group = el => el && el!==poster && el.matches(GROUPS);
function shown(id,value) { $(id).hidden=!value; }
function notice(message) { $("status").textContent=message; }
const label = el => short(el.getAttribute("data-label") || (el.localName==="text" && el.textContent) ||
  el.querySelector(":scope > title")?.textContent || el.id?.replaceAll("-"," ") || el.localName);
function role(el) {return el.dataset.role || (group(el)?"group":el.localName);}
function numeric(el,name,fallback) {
  const str=el.getAttribute(name); return finite(str)&&str!==null ? +str : fallback;
}
function styleValue(el,name) {return getComputedStyle(el).getPropertyValue(name).trim();}
function validLink(value) {
  if(!value) return "";
  if(/^(?:\.{1,2}\/)?[a-zA-Z0-9_-]+(?:[./?#][a-zA-Z0-9_./?&=%#-]*)?$/.test(value))return value;
  try { const url=new URL(value); return ["https:","http:","mailto:"].includes(url.protocol)?url.href:null; }
  catch { return null; }
}
function prepared(source, importing=false) {
  const doc=new DOMParser().parseFromString(source,"image/svg+xml");
  const root=doc.documentElement;
  if(root.localName!=="svg" || root.namespaceURI!==NS || doc.querySelector("parsererror"))
    throw Error("Choose a valid SVG document.");
  if(importing) {
    // Accept static local artwork. Remove executable elements and external embedded resources.
    const allow=new Set(("svg g a defs title desc metadata symbol use rect circle ellipse line polygon polyline path image " +
      "text tspan textPath linearGradient radialGradient stop clipPath mask pattern marker").split(" "));
    for(const el of [root,...root.querySelectorAll("*")]) {
      if(!allow.has(el.localName)||el.namespaceURI!==NS){el.remove();continue;}
      if(el.localName==="metadata")el.textContent=el.textContent;
      for(const a of [...el.attributes]) {
        const v=a.value.trim(), key=a.localName.toLowerCase();
        if(key.startsWith("on") || ["base","target","tabindex"].includes(key)){el.removeAttributeNode(a);continue;}
        if(key==="href") {
          const permitted=v.startsWith("#") || (el.localName==="a" && validLink(v)) ||
            (el.localName==="image" && /^data:image\/(png|jpeg|webp);base64,[a-z\d+/=\s]+$/i.test(v));
          if(!permitted) el.removeAttributeNode(a);
        } else if((/url\s*\(/i.test(v) && !/^(?:url\(\s*['"]?#[\w.:-]+['"]?\s*\))(?:\s+(?:none|#[\da-f]+))?$/i.test(v) && key!=="style") ||
          (key==="style" && /url\s*\(|@import|expression|\\|\/\*/i.test(v))) el.removeAttributeNode(a);
      }
    }
  }
  if(!root.hasAttribute("viewBox")) throw Error("The SVG needs a viewBox for poster editing.");
  const vb=root.getAttribute("viewBox").split(/[\s,]+/).map(Number);
  if(vb.length!==4||vb.some(v=>!Number.isFinite(v))||vb[2]<=0||vb[3]<=0) throw Error("Invalid SVG viewBox.");
  return document.importNode(root,true);
}
function mount(root) {
  if(drag) cancelDrag();
  selected=null; poster=root;
  $("artwork").replaceChildren(root);
  select(null); rebuildLayers(); setZoom(zoom); generation++;
  $("downloadSvg").disabled=false;$("downloadPdf").disabled=false;
}
function cleanClone() {return poster.cloneNode(true);}
function serialize(node=poster) {return new XMLSerializer().serializeToString(node);}
function wireDB() {
  return new Promise((resolve,reject)=>{
    const req=indexedDB.open("satb2-neuroday-editor",1);
    req.onupgradeneeded=()=>req.result.createObjectStore("drafts");
    req.onsuccess=()=>resolve(req.result); req.onerror=()=>reject(req.error);
  });
}
function fromDB(key) {
  return new Promise((resolve,reject)=>{
    const req=db.transaction("drafts").objectStore("drafts").get(key);
    req.onsuccess=()=>resolve(req.result);req.onerror=()=>reject(req.error);
  });
}
function writeDraft(draft) {
  return new Promise((resolve,reject)=>{
    const tx=db.transaction("drafts","readwrite"), store=tx.objectStore("drafts");
    const req=store.get(storageKey);
    let conflict=false;
    req.onsuccess=()=>{
      const latest=req.result?.id ?? null;
      if(latest!==lastDraftID){conflict=true;tx.abort();return;}
      store.put(draft,storageKey);
    };
    tx.oncomplete=()=>{lastDraftID=draft.id;resolve();};
    tx.onabort=()=>reject(Error(conflict?"Another tab saved a newer draft. Download your SVG to keep these edits.":"Could not save. Download your SVG to keep your edits."));
    tx.onerror=()=>reject(tx.error);
  });
}
function saveSoon() {
  dirty=true;saveNumber++;clearTimeout(timer);
  if(!db||storageBlocked){notice("Draft in this tab · download to keep");return;}
  notice("Saving locally…");
  timer=setTimeout(flush,320);
}
function flush() {
  clearTimeout(timer);
  if(!dirty||!db||storageBlocked)return saveChain;
  const seq=saveNumber, content=serialize(), gen=generation;
  saveChain=saveChain.then(async()=>{
    if(seq!==saveNumber||gen!==generation)return;
    const id=crypto.randomUUID?.() || String(Date.now()+Math.random());
    await writeDraft({id,hash:sourceHash,svg:content});
    if(seq===saveNumber){dirty=false;notice("Saved in this browser");}
  }).catch(e=>{storageBlocked=true;notice(e.message||"Local save unavailable. Download your SVG.");});
  return saveChain;
}
document.addEventListener("visibilitychange",()=>{if(document.hidden)flush();});
window.addEventListener("beforeunload",e=>{if(dirty){flush();e.preventDefault();e.returnValue="";}});
function buttons(){ $("undo").disabled=!undoStack.length; $("redo").disabled=!redoStack.length; }
function adopt(action,alreadyApplied=false) {
  if(!alreadyApplied)action.forward();
  undoStack.push(action);if(undoStack.length>65)undoStack.shift();
  redoStack=[];buttons();afterChange(action.structure);
}
function afterChange(structure=false){
  if(selected&&!poster.contains(selected))selected=null;
  if(structure)rebuildLayers();
  refresh();saveSoon();
}
function attributeAction(changes) {
  const ops=changes.map(([el,key,value])=>({el,key,next:value==null?null:String(value),prev:el.getAttribute(key)})).filter(o=>o.prev!==o.next);
  if(!ops.length)return null;
  const set=(o,v)=>{if(v===null)o.el.removeAttribute(o.key);else o.el.setAttribute(o.key,v);};
  return {forward:()=>ops.forEach(o=>set(o,o.next)),back:()=>[...ops].reverse().forEach(o=>set(o,o.prev))};
}
function attrs(changes) {const action=attributeAction(changes);if(action)adopt(action);}
function undo(redo=false) {
  const source=redo?redoStack:undoStack, dest=redo?undoStack:redoStack;
  const a=source.pop();if(!a)return;(redo?a.forward:a.back)();dest.push(a);buttons();afterChange(a.structure);
}
function eligible(el) {
  return el && el!==poster && poster.contains(el) && el.matches(OBJECTS) &&
   !el.closest("defs,clipPath,mask,marker,pattern,symbol");
}
function containing(el) {
  let p=el?.parentElement;
  while(p&&p!==poster){if(group(p))return p;p=p.parentElement;}
  return null;
}
function clicked(target) {
  if(!target||!poster.contains(target))return null;
  let el=target.closest("tspan")?.closest("text") || target.closest(OBJECTS);
  if(!eligible(el))return null;
  if(mode==="group") {
    const semantic=el.closest('[data-role="figure"],[data-role="section"],[data-role="dotplot"],[data-role="button"],g[data-label],a[data-label]');
    el=semantic|| (group(el)?el:containing(el))||el;
  }
  return el;
}
function select(el) {
  selected=eligible(el)?el:null;
  refresh();
}
function domMatrix(m){return new DOMMatrix([m.a,m.b,m.c,m.d,m.e,m.f]);}
function world(el) {return domMatrix(poster.getScreenCTM()).inverse().multiply(domMatrix(el.getScreenCTM()));}
function bounds(el) {
  const box=el.getBBox(),matrix=world(el);
  const corners=[[box.x,box.y],[box.x+box.width,box.y],[box.x,box.y+box.height],[box.x+box.width,box.y+box.height]]
    .map(([x,y])=>new DOMPoint(x,y).matrixTransform(matrix));
  const xs=corners.map(p=>p.x),ys=corners.map(p=>p.y);
  return {x:Math.min(...xs),y:Math.min(...ys),w:Math.max(...xs)-Math.min(...xs),h:Math.max(...ys)-Math.min(...ys)};
}
function controlsBox(){
  const sel=$("selection");if(!selected||!poster.contains(selected)){sel.hidden=true;return;}
  const box=selected.getBoundingClientRect(), paper=$("paper").getBoundingClientRect();
  if(!Number.isFinite(box.x)){sel.hidden=true;return;}
  sel.hidden=false;sel.style.left=(box.left-paper.left)+"px";sel.style.top=(box.top-paper.top)+"px";
  sel.style.width=Math.max(2,box.width)+"px";sel.style.height=Math.max(2,box.height)+"px";
}
function hexColor(value,fallback="#24383a"){
  if(/^#[0-9a-f]{6}$/i.test(value))return value;
  if(/^#[0-9a-f]{3}$/i.test(value))return "#"+value.slice(1).split("").map(x=>x+x).join("");
  const a=value?.match(/rgba?\(\s*(\d+(?:\.\d+)?)[, ]+\s*(\d+(?:\.\d+)?)[, ]+\s*(\d+(?:\.\d+)?)/i);
  return a?"#"+a.slice(1,4).map(x=>Math.min(255,Math.round(+x)).toString(16).padStart(2,"0")).join(""):fallback;
}
function textNode(el){return el.matches("text")?el:el.matches('[data-role="button"]')?el.querySelector("text"):null;}
function linkNode(el){return el.closest("a") || (el.dataset.role==="button"&&el.querySelector("a")) || null;}
function marker(el){return el?.dataset.role==="marker-dot";}
function shapeKind(el){return el.getAttribute("data-editor-shape")||el.localName;}
function dotMeta(el,key,fallback) {
  const owner=el.closest("[data-"+key+"]");return owner?.getAttribute("data-"+key)??fallback;
}
function setVal(id,val){$(id).value=val??"";}
function pointValue(el){return el?.dataset.role==="xy-point";}
function refresh(){
  shown("empty",!selected);shown("props",!!selected);controlsBox();
  for(const b of $("layersList").querySelectorAll(".active"))b.classList.remove("active");
  if(!selected) return;
  const activeIndex=cachedLayers.indexOf(selected);
  if(activeIndex>=0)$("layersList").querySelector('[data-index="'+activeIndex+'"]')?.classList.add("active");
  $("kind").textContent=role(selected);$("name").textContent=label(selected);
  $("parent").disabled=!containing(selected);
  const b=bounds(selected); selectionBox=b;
  for(const [id,val] of [["posX",b.x],["posY",b.y],["width",b.w],["height",b.h]])setVal(id,round(val));
  $("width").disabled=b.w<1e-5;$("height").disabled=b.h<1e-5;
  const t=textNode(selected);shown("textSect",!!t);
  if(t){setVal("copy",textRuns(t).map(n=>n.data).join(""));setVal("font",styleValue(t,"font-family").split(",")[0].replace(/['"]/g,""));
    setVal("fontSize",round(parseFloat(styleValue(t,"font-size"))));setVal("weight",styleValue(t,"font-weight"));}
  const a=linkNode(selected);shown("linkSect",!!a);if(a)setVal("url",a.getAttribute("href")||a.getAttributeNS(XL,"href")||"");
  shown("pointSect",pointValue(selected));shown("imageSect",selected.localName==="image");
  if(pointValue(selected)){setVal("pointX",selected.dataset.xValue);setVal("pointY",selected.dataset.yValue);}
  const isDot=marker(selected),dots=group(selected)?selected.querySelectorAll('[data-role="marker-dot"]'):[];
  shown("dotSect",isDot);shown("plotSect",!!dots.length);if(isDot){setVal("detection",selected.dataset.detection);setVal("expression",selected.dataset.expression);setVal("maxRadius",dotMeta(selected,"radius-max","7"));}
  if(dots.length){setVal("plotRadius",dotMeta(dots[0],"radius-max","7"));setVal("plotLow",hexColor(dotMeta(dots[0],"color-low","#f1f1eb")));setVal("plotHigh",hexColor(dotMeta(dots[0],"color-high","#08695f")));}
  const sk=shapeKind(selected);shown("shapeSect",["circle","ellipse","rect","diamond","triangle"].includes(sk));setVal("shape",sk);
  const sample=group(selected)?[...selected.querySelectorAll(LEAVES)].find(n=>styleValue(n,"fill")!=="none")||selected:selected;
  setVal("fill",hexColor(styleValue(sample,"fill")));setVal("stroke",hexColor(styleValue(sample,"stroke"),"#254849"));
  setVal("strokeWidth",round(parseFloat(styleValue(sample,"stroke-width"))||0));
  setVal("opacity",round(parseFloat(styleValue(selected,"opacity"))*100));
  const code=selected.matches("path")?"d":selected.matches("polygon,polyline")?"points":null;
  shown("geometrySect",!!code);if(code)setVal("geometry",selected.getAttribute(code));
}
function rebuildLayers(){
  if(!poster)return;
  // Avoid constructing thousands of DOM buttons; all editable primitives remain reachable on the canvas.
  cachedLayers=[...poster.querySelectorAll(OBJECTS)].filter(eligible);
  renderLayers();
}
function renderLayers(){
  const q=$("layerSearch").value.trim().toLowerCase(), host=$("layersList");
  host.replaceChildren();let shownCount=0, total=0;
  for(let i=0;i<cachedLayers.length;i++){
    const el=cachedLayers[i];
    if(!q && !(group(el)&&(el.hasAttribute("data-label")||el.hasAttribute("id")))) continue;
    const name=label(el);
    if(q&&!name.toLowerCase().includes(q))continue;
    total++;if(shownCount>=160)continue;
    const b=document.createElement("button");b.type="button";b.textContent=name;b.title=role(el)+" · "+name;
    b.dataset.index=i;b.className=el.dataset.role==="section"?"section":"nested";if(el===selected)b.classList.add("active");
    host.append(b);shownCount++;
  }
  if(!shownCount || total>shownCount){const p=document.createElement("p");p.className="layer-note";p.textContent=shownCount?"Refine your search to see more objects.":q?"No matching objects.":"Use the canvas or search for text.";host.append(p);}
}
function setLayers(open){shown("layers",open);$("layersToggle").setAttribute("aria-expanded",String(open));requestAnimationFrame(setSize);}
$("layersToggle").onclick=()=>setLayers($("layers").hidden);$("closeLayers").onclick=()=>setLayers(false);$("browse").onclick=()=>setLayers(true);
$("layersList").onclick=e=>{const b=e.target.closest("button[data-index]");if(!b)return;select(cachedLayers[+b.dataset.index]);
  if(selected)selected.scrollIntoView({block:"nearest",inline:"nearest",behavior:"instant"});controlsBox();};
$("layerSearch").oninput=renderLayers;
function matrixText(m){return "matrix("+[m.a,m.b,m.c,m.d,m.e,m.f].map(v=>+v.toFixed(8)).join(" ")+")";}
function transformState(el) {
  const parent=el.parentElement, p=world(parent);
  return {parent:p,local:p.inverse().multiply(world(el)),start:el.getAttribute("transform"),box:bounds(el)};
}
function transformed(state,operation){return matrixText(state.parent.inverse().multiply(operation).multiply(state.parent).multiply(state.local));}
function scaleAt(x,y,sx,sy){return new DOMMatrix().translate(x,y).scale(sx,sy).translate(-x,-y);}
function moveNumber(id) {
  if(!selected || !finite($(id).value))return refresh();
  const val=Number($(id).value),state=transformState(selected),b=state.box;
  let m;
  if(id==="posX"||id==="posY")m=new DOMMatrix().translate(id==="posX"?val-b.x:0,id==="posY"?val-b.y:0);
  else {
    if(val<=0)return refresh();
    const factor=val/(id==="width"?b.w:b.h);if(!Number.isFinite(factor)||factor>1e5)return refresh();
    m=scaleAt(b.x,b.y,id==="width"||$("ratio").checked?factor:1,id==="height"||$("ratio").checked?factor:1);
  }
  attrs([[selected,"transform",transformed(state,m)]]);
}
for(const id of ["posX","posY","width","height"])$(id).onchange=()=>moveNumber(id);
function textRuns(el) {
  const out=[];
  function walk(parent){for(const child of parent.childNodes){
    if(child.nodeType===3)out.push(child);
    else if(child.nodeType===1&&["tspan","textPath"].includes(child.localName))walk(child);
  }}
  walk(el);return out;
}
function textAction(el,value) {
  let runs=textRuns(el);if(!runs.length){el.append(document.createTextNode(""));runs=textRuns(el);}
  const orig=runs.map(n=>n.data),prev=orig.join("");if(value===prev)return;
  let start=0,end=0;
  while(start<prev.length && start<value.length && prev[start]===value[start])start++;
  while(end<prev.length-start && end<value.length-start && prev[prev.length-1-end]===value[value.length-1-end])end++;
  const insert=value.slice(start,value.length-end);
  let pos=0,inserted=false;
  const next=orig.map((s,i)=>{
    const left=pos;pos+=s.length;
    const prefix=s.slice(0,Math.max(0,Math.min(s.length,start-left)));
    const suffix=s.slice(Math.max(0,Math.min(s.length,prev.length-end-left)));
    const getsInsert=!inserted && start>=left && (start<pos||i===orig.length-1);
    if(getsInsert)inserted=true;
    return prefix+(getsInsert?insert:"")+suffix;
  });
  const set=vals=>vals.forEach((v,i)=>runs[i].data=v);
  return {forward:()=>set(next),back:()=>set(orig),structure:true};
}
$("copy").onchange=()=>{if(selected&&textNode(selected)){const action=textAction(textNode(selected),$("copy").value.replace(/[\r\n]+/g," "));if(action)adopt(action);}};
function styleChanges(el,prop,val) {
  // Inline CSS inside imported diagrams takes precedence over presentation attributes.
  const copy=el.cloneNode(false);if(copy.style.getPropertyValue(prop))copy.style.removeProperty(prop);
  const out=[];
  if(copy.getAttribute("style")!==el.getAttribute("style"))out.push([el,"style",copy.getAttribute("style")||null]);
  out.push([el,prop,val]);
  return out;
}
for(const [id,prop] of [["font","font-family"],["fontSize","font-size"],["weight","font-weight"]]) {
  $(id).onchange=()=>{const el=selected&&textNode(selected),v=$(id).value.trim();if(!el||!v||id==="fontSize"&&(!finite(v)||+v<=0))return refresh();attrs(styleChanges(el,prop,v));};
}
function paint(prop,value) {
  if(!selected)return;
  const nodes=group(selected)?[...selected.querySelectorAll(LEAVES)].filter(n=>
    prop==="fill"?styleValue(n,"fill")!=="none"&&!n.matches("line,image"):prop==="stroke-width"?styleValue(n,"stroke")!=="none":styleValue(n,"stroke")!=="none"): [selected];
  if(!nodes.length){if(group(selected)&&prop==="stroke"){for(const n of selected.querySelectorAll(LEAVES))if(!n.matches("image"))nodes.push(n);}else if(!group(selected))nodes.push(selected);}
  const changes=nodes.flatMap(el=>styleChanges(el,prop,value));attrs(changes);
}
$("fill").onchange=()=>paint("fill",$("fill").value);$("stroke").onchange=()=>paint("stroke",$("stroke").value);
$("noFill").onclick=()=>paint("fill","none");$("noStroke").onclick=()=>paint("stroke","none");
$("strokeWidth").onchange=()=>{if(finite($("strokeWidth").value)&&+$("strokeWidth").value>=0)paint("stroke-width",$("strokeWidth").value);else refresh();};
$("opacity").onchange=()=>{if(selected&&finite($("opacity").value))attrs(styleChanges(selected,"opacity",Math.max(0,Math.min(100,+$("opacity").value))/100));else refresh();};
$("url").onchange=()=>{
  if(!selected)return;const a=linkNode(selected),value=validLink($("url").value.trim());
  if(!a||value===null){notice("Enter a full https, http or mailto destination.");refresh();return;}
  const old=a.getAttribute("href"),oldXL=a.getAttributeNS(XL,"href");
  const set=(u,x)=>{if(u===null)a.removeAttribute("href");else a.setAttribute("href",u);if(x===null)a.removeAttributeNS(XL,"href");else a.setAttributeNS(XL,"xlink:href",x);};
  adopt({forward:()=>set(value||null,null),back:()=>set(old,oldXL)});
};
function mix(a,b,t){
  a=hexColor(a);b=hexColor(b);t=Math.min(1,Math.max(0,t));
  return "#"+[1,3,5].map(i=>Math.round(parseInt(a.slice(i,i+2),16)*(1-t)+parseInt(b.slice(i,i+2),16)*t).toString(16).padStart(2,"0")).join("");
}
function dotColor(el,value,lowOverride,highOverride) {
  const low=lowOverride||dotMeta(el,"color-low","#f1f1eb"),high=highOverride||dotMeta(el,"color-high","#08695f");
  const mid=dotMeta(el,"color-mid",""),lo=+dotMeta(el,"expression-min","0"),hi=+dotMeta(el,"expression-max","1");
  if(hi<=lo)return high;
  const split=+dotMeta(el,"expression-mid",lo<0&&hi>0?"0":String((lo+hi)/2));
  if(mid&&split>lo&&split<hi)return value<=split?mix(low,mid,(value-lo)/(split-lo)):mix(mid,high,(value-split)/(hi-split));
  return mix(low,high,(value-lo)/(hi-lo));
}
function symbolGeometry(shape,cx,cy,w,h) {
  const x=cx-w/2,y=cy-h/2;
  if(shape==="circle")return {cx,cy,r:Math.min(w,h)/2};
  if(shape==="ellipse")return {cx,cy,rx:w/2,ry:h/2};
  if(shape==="rect")return {x,y,width:w,height:h};
  const pts=shape==="diamond"?[[cx,y],[x+w,cy],[cx,y+h],[x,cy]]:[[cx,y],[x+w,y+h],[x,y+h]];
  return {points:pts.map(pt=>pt.map(round).join(",")).join(" ")};
}
function markerCenter(el){
  if(el.hasAttribute("data-center-x"))return [+el.getAttribute("data-center-x"),+el.getAttribute("data-center-y")];
  if(el.matches("circle,ellipse"))return [numeric(el,"cx",0),numeric(el,"cy",0)];
  const b=el.getBBox();return [b.x+b.width/2,b.y+b.height/2];
}
function sizeDot(el,det,rad) {
  const [cx,cy]=markerCenter(el),diam=2*rad*Math.sqrt(Math.min(100,Math.max(0,det))/100);
  return Object.entries(symbolGeometry(shapeKind(el),cx,cy,diam,diam)).map(([k,v])=>[el,k,v]);
}
function combined(actions){
  const list=actions.filter(Boolean);if(!list.length)return;
  adopt({forward:()=>list.forEach(a=>a.forward()),back:()=>[...list].reverse().forEach(a=>a.back()),structure:list.some(a=>a.structure)});
}
function linkedValue(el,value){
  const id=el.dataset.valueLabel?.replace(/^#/,"");
  if(!id)return null;
  const n=poster.querySelector("#"+CSS.escape(id));
  return n?.matches("text")?textAction(n,(value>=0?"+":"−")+Math.abs(value).toFixed(2)):null;
}
function toolTip(el,content){
  const n=el.querySelector(":scope > title");if(!n)return null;
  const before=n.textContent;if(before===content)return null;
  return {forward:()=>{n.textContent=content;},back:()=>{n.textContent=before;}};
}
function dataValue(id) {
  if(!marker(selected)||!finite($(id).value))return refresh();const value=+$(id).value,el=selected;
  if(id==="expression") {
    const lo=+dotMeta(el,"expression-min",0),hi=+dotMeta(el,"expression-max",1);
    const clamped=Math.max(lo,Math.min(hi,value));
    combined([attributeAction([[el,"data-expression",value],...styleChanges(el,"fill",dotColor(el,clamped))]),
      linkedValue(el,value),toolTip(el,label(el)+" · detected "+el.dataset.detection+"% · expression "+value)]);
  } else {
    const det=id==="detection"?Math.min(100,Math.max(0,value)):+el.dataset.detection;
    const rad=id==="maxRadius"?Math.max(.1,value):+dotMeta(el,"radius-max",7);
    combined([attributeAction([[el,"data-detection",det],[el,"data-radius-max",rad],...sizeDot(el,det,rad)]),
      toolTip(el,label(el)+" · detected "+det+"% · expression "+el.dataset.expression)]);
  }
}
for(const [id,axis] of [["pointX","x"],["pointY","y"]])$(id).onchange=()=>{
  const el=selected;if(!pointValue(el)||!finite($(id).value))return refresh();
  const value=+$(id).value,origin=+el.getAttribute("data-"+axis+"-origin"),scale=+el.getAttribute("data-"+axis+"-scale");
  const x=axis==="x"?value:+el.dataset.xValue,y=axis==="y"?value:+el.dataset.yValue;
  const entries=[[el,"data-"+axis+"-value",value],[el,axis==="x"?"cx":"cy",origin+scale*value]];
  const color=y>x?el.dataset.colorY:el.dataset.colorX;if(color)entries.push(...styleChanges(el,"fill",color));
  combined([attributeAction(entries),toolTip(el,label(el)+" · x "+x+" · y "+y)]);
};
for(const id of ["detection","expression","maxRadius"])$(id).onchange=()=>dataValue(id);
for(const id of ["plotRadius","plotLow","plotHigh"])$(id).onchange=()=>{
  if(!group(selected))return; const arr=[];
  for(const el of selected.querySelectorAll('[data-role="marker-dot"]')) {
    const v=$(id).value;
    if(id==="plotRadius"){
      if(!finite(v)||+v<=0)return refresh();arr.push([el,"data-radius-max",v],...sizeDot(el,+el.dataset.detection,+v));
    }else{
      const key=id==="plotLow"?"data-color-low":"data-color-high";
      arr.push([el,key,v],...styleChanges(el,"fill",dotColor(el,+el.dataset.expression,id==="plotLow"?v:null,id==="plotHigh"?v:null)));
    }
  }
  attrs(arr);
};
$("shape").onchange=()=>{
  const el=selected;if(!el)return;
  const desired=$("shape").value,b=el.getBBox(),cx=b.x+b.width/2,cy=b.y+b.height/2;
  const tag=["diamond","triangle"].includes(desired)?"polygon":desired;
  const node=svgNode(tag),geom=new Set(["x","y","cx","cy","r","rx","ry","width","height","points"]);
  for(const a of el.attributes)if(!geom.has(a.name))node.setAttributeNS(a.namespaceURI,a.name,a.value);
  node.setAttribute("data-editor-shape",desired);
  if(marker(el)){node.setAttribute("data-center-x",cx);node.setAttribute("data-center-y",cy);}
  for(const [k,v] of Object.entries(symbolGeometry(desired,cx,cy,b.width,b.height)))node.setAttribute(k,v);
  for(const child of el.childNodes)node.append(child.cloneNode(true));
  adopt({forward(){el.replaceWith(node);selected=node;},back(){node.replaceWith(el);selected=el;},structure:true});
};
$("geometry").onchange=()=>{
  if(!selected)return;const k=selected.matches("path")?"d":"points",v=$("geometry").value.trim();
  if(k==="d"&&v&&!/^[mMzZlLhHvVcCsSqQtTaA0-9eE.,+\s-]+$/.test(v) || k==="points"&&v&&!/^[0-9eE.,+\s-]+$/.test(v)){notice("Enter valid SVG path coordinates.");return refresh();}
  attrs([[selected,k,v]]);
};
function pointerWorld(e){return new DOMPoint(e.clientX,e.clientY).matrixTransform(poster.getScreenCTM().inverse());}
function startDrag(e,corner=null) {
  if(e.button!==0||!poster||drag)return;
  let el=corner===null?clicked(e.target):selected;
  if(e.altKey&&corner===null){
    const stack=[...new Set(document.elementsFromPoint(e.clientX,e.clientY).filter(n=>poster.contains(n)).map(clicked).filter(Boolean))];
    if(stack.length)el=stack[(stack.indexOf(selected)+1)%stack.length];
  }
  if(!el){if(!$("paper").contains(e.target)||e.target===$("artwork")||e.target===poster)select(null);return;}
  select(el);
  const state=transformState(el);
  drag={el,state,corner,id:e.pointerId,origin:pointerWorld(e),client:[e.clientX,e.clientY], moved:false};
  if(e.pointerType!=="touch")e.preventDefault();
  $("stage").setPointerCapture(e.pointerId);
}
function duringDrag(e) {
  if(!drag||e.pointerId!==drag.id)return;
  const {el,state,corner,origin}=drag,p=pointerWorld(e);
  if(!drag.moved && Math.hypot(e.clientX-drag.client[0],e.clientY-drag.client[1])<4)return;
  drag.moved=true;$("stage").classList.add("is-dragging");
  let m;const dx=p.x-origin.x,dy=p.y-origin.y;
  if(corner===null)m=new DOMMatrix().translate(dx,dy);
  else{
    const b=state.box,l=corner===0||corner===3,t=corner===0||corner===1;
    let sx=b.w>0?Math.max(.015,1+(l?-dx:dx)/b.w):1,sy=b.h>0?Math.max(.015,1+(t?-dy:dy)/b.h):1;
    if($("ratio").checked&&!e.shiftKey){const f=Math.abs(sx-1)>Math.abs(sy-1)?sx:sy;sx=f;sy=f;}
    m=scaleAt(l?b.x+b.w:b.x,t?b.y+b.h:b.y,sx,sy);
  }
  el.setAttribute("transform",transformed(state,m));controlsBox();
}
function cancelDrag() {
  if(!drag)return;
  const {el,state,id}=drag;drag=null;
  if(state.start==null)el.removeAttribute("transform");else el.setAttribute("transform",state.start);
  if($("stage").hasPointerCapture(id))$("stage").releasePointerCapture(id);
  $("stage").classList.remove("is-dragging");refresh();
}
function endDrag(e) {
  if(!drag||drag.id!==e.pointerId)return;const {el,state,id,moved}=drag;const val=el.getAttribute("transform");drag=null;
  if($("stage").hasPointerCapture(id))$("stage").releasePointerCapture(id);
  $("stage").classList.remove("is-dragging");
  if(moved&&val!==state.start)adopt({forward:()=>el.setAttribute("transform",val),back:()=>state.start===null?el.removeAttribute("transform"):el.setAttribute("transform",state.start)},true);
}
$("artwork").addEventListener("pointerdown",e=>startDrag(e));
$("selection").querySelectorAll("i").forEach((h,i)=>h.onpointerdown=e=>startDrag(e,i));
$("stage").onpointermove=duringDrag;$("stage").onpointerup=endDrag;$("stage").onpointercancel=cancelDrag;
$("artwork").onclick=e=>{if(e.target.closest("a"))e.preventDefault();};
$("artwork").ondblclick=e=>{
  if(!poster)return;const el=e.target.closest("text")||clicked(e.target);if(!el)return;
  select(el);if(textNode(el)){$("copy").focus();$("copy").select();}
};
function setMode(next){mode=next;$("modeObject").setAttribute("aria-pressed",String(mode==="object"));$("modeGroup").setAttribute("aria-pressed",String(mode==="group"));}
$("modeObject").onclick=()=>setMode("object");$("modeGroup").onclick=()=>setMode("group");
$("parent").onclick=()=>select(containing(selected));
function applyTransform(el,m){const st=transformState(el);attrs([[el,"transform",transformed(st,m)]]);}
function positionOfPlacement(){
  const rect=$("stage").getBoundingClientRect(),paper=poster.getBoundingClientRect();
  const x=(Math.max(rect.left,paper.left)+Math.min(rect.right,paper.right))/2;
  const y=(Math.max(rect.top,paper.top)+Math.min(rect.bottom,paper.bottom))/2;
  return new DOMPoint(x,y).matrixTransform(poster.getScreenCTM().inverse());
}
function addNode(el,parent=poster,next=null) {
  const prev=selected;
  adopt({forward(){parent.insertBefore(el,next?.parentNode===parent?next:null);selected=el;},back(){el.remove();selected=prev;},structure:true});
}
function addElement(type) {
  if(!poster||!type)return;$("addType").value="";
  if(type==="image"){imageTarget=null;$("imageUpload").click();return;}
  const p=positionOfPlacement(),x=round(p.x),y=round(p.y);let n;
  if(type==="text"){
    n=svgNode("text",{x,y,fill:"#17383c","font-family":"Arial, sans-serif","font-size":22});n.textContent="New text";
  }else if(type==="button"){
    n=svgNode("a",{"data-role":"button","data-label":"Link button",href:"https://ali-mujahid.github.io/satb2_snRNA/",transform:"translate("+x+" "+y+")"});
    n.append(svgNode("rect",{width:145,height:36,rx:5,fill:"#11675f"}));
    const t=svgNode("text",{x:72.5,y:24,fill:"#ffffff","text-anchor":"middle","font-family":"Arial, sans-serif","font-size":17});t.textContent="Open atlas";n.append(t);
  }else if(type==="line")n=svgNode("line",{x1:x,y1:y,x2:x+100,y2:y,stroke:"#17383c","stroke-width":2});
  else n=type==="rect"?svgNode("rect",{x,y,width:85,height:45,fill:"#80b5aa"}):svgNode("ellipse",{cx:x,cy:y,rx:40,ry:23,fill:"#80b5aa"});
  addNode(n);
  if(type==="text"){$("copy").focus();$("copy").select();}
}
$("addType").onchange=e=>addElement(e.target.value);
$("replaceImage").onclick=()=>{imageTarget=selected;$("imageUpload").click();};
$("imageUpload").onchange=async e=>{
  const oldImage=imageTarget;imageTarget=null;
  const file=e.target.files[0];e.target.value="";if(!file)return;
  if(!["image/png","image/jpeg","image/webp"].includes(file.type)){notice("Choose a PNG, JPEG or WebP image.");return;}
  try{
    const url=await new Promise((ok,err)=>{const r=new FileReader();r.onload=()=>ok(r.result);r.onerror=err;r.readAsDataURL(file);});
    const image=new Image();image.src=url;await image.decode();
    if(oldImage&&poster.contains(oldImage)){
      const oldHref=oldImage.getAttribute("href"),oldXL=oldImage.getAttributeNS(XL,"href"),oldLabel=oldImage.getAttribute("data-label");
      const set=(u,x,name)=>{if(u===null)oldImage.removeAttribute("href");else oldImage.setAttribute("href",u);if(x===null)oldImage.removeAttributeNS(XL,"href");else oldImage.setAttributeNS(XL,"xlink:href",x);if(name===null)oldImage.removeAttribute("data-label");else oldImage.setAttribute("data-label",name);};
      adopt({forward:()=>set(url,null,file.name),back:()=>set(oldHref,oldXL,oldLabel),structure:true});return;
    }
    const p=positionOfPlacement(),w=Math.min(300,poster.viewBox.baseVal.width/3);
    addNode(svgNode("image",{x:round(p.x),y:round(p.y),width:w,height:w*image.height/image.width,href:url,"data-label":file.name}));
  }catch{notice("Could not open that image.");}
};
function removeSelection(){
  const n=selected;if(!n)return;const p=n.parentNode,before=n.nextSibling;
  adopt({forward(){n.remove();selected=null;},back(){p.insertBefore(n,before?.parentNode===p?before:null);selected=n;},structure:true});
}
function cloneSelection(){
  const n=selected;if(!n)return;const clone=n.cloneNode(true);
  const suffix="-copy-"+(crypto.randomUUID?.()||Math.random().toString(36).slice(2,10)).slice(0,8);
  const ids=new Map();
  for(const x of [clone,...clone.querySelectorAll("[id]")])if(x.id){ids.set(x.id,x.id+suffix);x.id+=suffix;}
  for(const x of [clone,...clone.querySelectorAll("*")])for(const a of [...x.attributes]){
    if(a.name==="id")continue;
    if(a.localName==="href"&&a.value.startsWith("#")&&ids.has(a.value.slice(1)))a.value="#"+ids.get(a.value.slice(1));
    else if(a.value.includes("url(#"))a.value=a.value.replace(/url\(#([^)]+)\)/g,(m,id)=>ids.has(id)?"url(#"+ids.get(id)+")":m);
  }
  const state=transformState(n);clone.setAttribute("transform",transformed(state,new DOMMatrix().translate(14,14)));
  addNode(clone,n.parentNode,n.nextSibling);
}
$("duplicate").onclick=cloneSelection;$("remove").onclick=removeSelection;
$("undo").onclick=()=>undo();$("redo").onclick=()=>undo(true);
document.addEventListener("keydown",e=>{
  if(e.target.closest("input,textarea,select,dialog,[contenteditable=true]"))return;
  if(!poster)return;
  const cmd=e.ctrlKey||e.metaKey,key=e.key.toLowerCase();
  if(cmd&&key==="z"){e.preventDefault();undo(e.shiftKey);return;}
  if(cmd&&key==="y"){e.preventDefault();undo(true);return;}
  if(cmd&&key==="d"&&selected){e.preventDefault();cloneSelection();return;}
  if(key==="escape"){if(drag)cancelDrag();else select(null);return;}
  if(!selected||cmd)return;
  if(["backspace","delete"].includes(key)){e.preventDefault();removeSelection();}
  if(key.startsWith("arrow")){
    const dx=(key==="arrowright")-(key==="arrowleft"),dy=(key==="arrowdown")-(key==="arrowup");
    if(dx||dy){e.preventDefault();applyTransform(selected,new DOMMatrix().translate(dx*(e.shiftKey?10:1),dy*(e.shiftKey?10:1)));}
  }
});
function baseWidth(){
  const stage=$("stage"),cs=getComputedStyle(stage);return Math.max(110,stage.clientWidth-parseFloat(cs.paddingLeft)-parseFloat(cs.paddingRight));
}
function setSize(){if(!poster)return;$("paper").style.width=round(baseWidth()*zoom)+"px";controlsBox();}
function setZoom(n){zoom=Math.min(8,Math.max(.15,n));$("zoomLevel").textContent=Math.round(zoom*100)+"%";setSize();}
$("zoomIn").onclick=()=>setZoom(zoom*1.25);$("zoomOut").onclick=()=>setZoom(zoom/1.25);
$("zoomLevel").onclick=()=>setZoom(1);$("fit").onclick=()=>{
  if(!poster)return;const v=poster.viewBox.baseVal,s=$("stage"),cs=getComputedStyle(s);
  setZoom(Math.min(1,(s.clientHeight-parseFloat(cs.paddingTop)-parseFloat(cs.paddingBottom))*v.width/v.height/baseWidth()));s.scrollTo(0,0);
};
new ResizeObserver(setSize).observe($("stage"));
function saveFile(blob,name) {
  const url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download=name;a.click();
  setTimeout(()=>URL.revokeObjectURL(url),30000);
}
$("downloadSvg").onclick=()=>{if(poster)saveFile(new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n',serialize(cleanClone())],{type:"image/svg+xml;charset=utf-8"}),"SATB2-Neuroday-edited.svg");};
function printFallback(){
  $("printArea")?.remove();const area=document.createElement("div");area.id="printArea";area.style.display="none";area.append(cleanClone());document.body.append(area);
  const done=()=>area.remove();window.addEventListener("afterprint",done,{once:true});window.print();
}
let pdfScript;
function loadPdf(){
  if(typeof window.downloadPosterPDF==="function")return Promise.resolve(true);
  if(!pdfScript) pdfScript=new Promise(resolve=>{
    const tag=document.createElement("script");tag.src=new URL("./pdf-export.js",assetBase).href;tag.onload=()=>resolve(typeof window.downloadPosterPDF==="function");tag.onerror=()=>resolve(false);document.head.append(tag);
  });
  return pdfScript;
}
$("downloadPdf").onclick=async()=>{
  const btn=$("downloadPdf");if(!poster)return;btn.disabled=true;const caption=btn.textContent;btn.textContent="Preparing PDF…";
  try{
    if(await loadPdf())await window.downloadPosterPDF(cleanClone());else printFallback();
  }catch{notice("PDF export needs the browser print dialog.");printFallback();}
  finally{btn.disabled=false;btn.textContent=caption;}
};
$("openSvg").onclick=()=>$("svgUpload").click();
$("svgUpload").onchange=async e=>{
  const file=e.target.files[0];e.target.value="";if(!file)return;
  try{
    const newDoc=prepared(await file.text(),true),prior=poster;
    adopt({forward(){mount(newDoc);},back(){mount(prior);},structure:true});
  }catch(err){notice(err.message||"Could not open the SVG.");}
};
$("restore").onclick=()=>{if(poster)$("confirmReset").showModal();};$("cancelReset").onclick=()=>$("confirmReset").close();
$("confirmRestore").onclick=()=>{
  $("confirmReset").close();mount(prepared(original));undoStack=[];redoStack=[];buttons();saveSoon();
};
async function init(){
  try{
    const res=await fetch(new URL("./poster.svg",assetBase).href,{cache:"no-cache"});if(!res.ok)throw Error("The poster SVG could not be loaded.");
    original=await res.text();const canonical=prepared(original);
    try{
      if(crypto.subtle)sourceHash=Array.from(new Uint8Array(await crypto.subtle.digest("SHA-256",new TextEncoder().encode(original))),n=>n.toString(16).padStart(2,"0")).join("");
      db=await wireDB();const draft=await fromDB(storageKey);
      if(draft?.id)lastDraftID=draft.id;
      if(draft?.svg){try{mount(prepared(draft.svg,true));notice(draft.hash&&sourceHash&&draft.hash!==sourceHash?"Local draft from an earlier edition":"Local draft restored");return;}
        catch {notice("Saved draft could not be opened.");}}
    }catch{db=null;}
    mount(canonical);notice(db?"Original poster":"Download SVG to keep edits");
  }catch(error){$("paper").hidden=true;shown("loadError",true);$("loadError").textContent=error.message+" Reload the page to try again.";notice("Poster unavailable");}
}
init();
})();
