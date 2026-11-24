/* Vector A0 PDF export for the neuroday poster editor. Everything runs locally in the browser. */
(() => {
  "use strict";
  const assetBase = new URL("vendor/", document.currentScript?.src || document.baseURI);
  const asset = name => new URL(name, assetBase).href;
  const page = { width: 841 * 72 / 25.4, height: 1189 * 72 / 25.4 };
  const typefaces = new Map();
  let renderer;

  function script(name) {
    return new Promise((resolve, reject) => {
      const tag = document.createElement("script");
      tag.src = asset(name);
      tag.crossOrigin = "anonymous";
      tag.onload = resolve;
      tag.onerror = () => {
        tag.remove();
        reject(new Error("Could not load the PDF renderer. Check your connection and retry."));
      };
      document.head.append(tag);
    });
  }
  function loadRenderer() {
    if (!renderer) {
      renderer = (async () => {
        if (!window.jspdf?.jsPDF) await script("jspdf.umd.min.js");
        if (!window.svg2pdf?.svg2pdf) await script("svg2pdf.umd.min.js");
        if (!window.jspdf?.jsPDF || !window.svg2pdf?.svg2pdf) throw new Error("PDF renderer initialization failed.");
      })().catch(error => { renderer = undefined; throw error; });
    }
    return renderer;
  }
  const faceNames = {
    PosterSans: "LiberationSans",
    PosterSerif: "LiberationSerif",
    PosterMono: "LiberationMono",
    PosterSymbols: "DejaVuSans"
  };
  const variants = { normal:"Regular", bold:"Bold", italic:"Italic", bolditalic:"BoldItalic" };

  function familyOf(name) {
    const first = name.split(",")[0].toLowerCase();
    if (/courier|monospace|mono\b/.test(first)) return "PosterMono";
    if (/georgia|times|serif/.test(first) && !first.includes("sans")) return "PosterSerif";
    return "PosterSans";
  }
  function loadFont(family, style) {
    const key = family + "-" + style;
    if (typefaces.has(key)) return typefaces.get(key);
    const file = family === "PosterSymbols" ? "DejaVuSans.ttf" : faceNames[family] + "-" + variants[style] + ".ttf";
    const loading = (async () => {
      const response = await fetch(asset(file), {mode:"cors", credentials:"omit"});
      if (!response.ok) throw new Error("Could not load a PDF typeface. Please retry.");
      const bytes = new Uint8Array(await response.arrayBuffer());
      // Load the same outlines in the browser and the PDF. svg2pdf measures text
      // with the DOM when placing tspans and centered/right-aligned labels.
      const browserFace = new FontFace(family, bytes, {
        weight: style.includes("bold") ? "700" : "400",
        style: style.includes("italic") ? "italic" : "normal"
      });
      document.fonts.add(await browserFace.load());
      const parts = [];
      for (let i=0; i<bytes.length; i+=16384) {
        parts.push(String.fromCharCode(...bytes.subarray(i, i+16384)));
      }
      return { file, base64:btoa(parts.join("")), family, style };
    })().catch(error => { typefaces.delete(key); throw error; });
    typefaces.set(key, loading);
    return loading;
  }

  function prepareFonts(svg) {
    // Read all inherited styles before changing any of them. This also handles
    // a bold or italic parent with a single specially styled scientific symbol.
    const tasks = new Map();
    const records = Array.from(svg.querySelectorAll("text,tspan"), node => {
      const css = getComputedStyle(node);
      const special = /[⁰⁻↗]/.test(node.textContent);
      const family = special ? "PosterSymbols" : familyOf(css.fontFamily);
      const bold = !special && (css.fontWeight === "bold" || Number(css.fontWeight) >= 600);
      const italic = !special && /italic|oblique/.test(css.fontStyle);
      const style = bold ? (italic ? "bolditalic" : "bold") : (italic ? "italic" : "normal");
      tasks.set(family + "-" + style, () => loadFont(family, style));
      return {node, family, bold, italic};
    });
    for (const {node, family, bold, italic} of records) {
      node.setAttribute("font-family",family);
      node.setAttribute("font-weight",bold ? "bold" : "normal");
      node.setAttribute("font-style",italic ? "italic" : "normal");
      node.style.fontFamily=family;
      node.style.fontWeight=bold ? "bold" : "normal";
      node.style.fontStyle=italic ? "italic" : "normal";
    }
    return Promise.all([...tasks.values()].map(start => start()));
  }

  function normalizeVectorFeatures(svg) {
    // Object-bounding-box coordinates are fractions in the PDF shading model.
    // Keep gradient proportions identical for both "%" and unitless authoring.
    for (const g of svg.querySelectorAll("linearGradient,radialGradient")) {
      if (g.getAttribute("gradientUnits") === "userSpaceOnUse") continue;
      for (const attr of ["x1","y1","x2","y2","cx","cy","r","fx","fy"]) {
        const value = g.getAttribute(attr);
        if (value?.endsWith("%")) g.setAttribute(attr,String(parseFloat(value)/100));
      }
    }
    // svg2pdf does not consume SVG letter-spacing. Keep each run's edited
    // visual length and anchor, without splitting text into individual glyphs.
    for (const t of Array.from(svg.querySelectorAll("text"))) {
      const css=getComputedStyle(t);
      const spacing=parseFloat(css.letterSpacing);
      if (spacing && t.getNumberOfChars()) {
        const width=t.getComputedTextLength();
        t.style.letterSpacing="0px";
        t.setAttribute("letter-spacing","0");
        const untracked=t.getComputedTextLength();
        const anchor=t.x.baseVal.numberOfItems ? t.x.baseVal.getItem(0).value : 0;
        if (width>0 && untracked>0) {
          const prior=t.getAttribute("transform") || "";
          t.setAttribute("transform", prior+" translate("+anchor+" 0) scale("+(width/untracked)+" 1) translate("+(-anchor)+" 0)");
        }
      }
      // Paint-order:stroke puts the halo behind the fill in SVG. Use two
      // PDF-native text draws so short white-outlined UMAP labels stay legible.
      if (css.paintOrder?.startsWith("stroke") && css.stroke!=="none") {
        const outline=t.cloneNode(true);
        outline.removeAttribute("id");
        outline.style.fill="none";outline.setAttribute("fill","none");
        t.before(outline);
        t.style.stroke="none";t.setAttribute("stroke","none");
      }
    }
    for (const node of svg.querySelectorAll("[font-weight]")) {
      const w=node.getAttribute("font-weight");
      if (/^\d+$/.test(w)) node.setAttribute("font-weight",Number(w)>=600?"bold":"normal");
    }
  }

  // svg2pdf draws the linked artwork but does not generate PDF link annotations.
  function addLinkAnnotations(pdf, svg) {
    const viewport=svg.getBoundingClientRect();
    if(!viewport.width || !viewport.height)return;
    const xScale=page.width/viewport.width, yScale=page.height/viewport.height;
    for(const anchor of svg.querySelectorAll("a")) {
      const href=anchor.getAttribute("href") || anchor.getAttributeNS("http://www.w3.org/1999/xlink","href");
      if(!href || href.startsWith("#"))continue;
      let url;
      try {url=new URL(href,window.POSTER_ASSET_BASE || document.baseURI);}
      catch {continue;}
      if(!["https:","http:","mailto:"].includes(url.protocol))continue;
      const bounds=anchor.getBoundingClientRect();
      const left=Math.max(0,(bounds.left-viewport.left)*xScale);
      const top=Math.max(0,(bounds.top-viewport.top)*yScale);
      const right=Math.min(page.width,(bounds.right-viewport.left)*xScale);
      const bottom=Math.min(page.height,(bounds.bottom-viewport.top)*yScale);
      if(right>left && bottom>top)pdf.link(left,top,right-left,bottom-top,{url:url.href});
    }
  }

  window.downloadPosterPDF = async function downloadPosterPDF(svgElement) {
    if (!svgElement || svgElement.localName !== "svg") throw new TypeError("A poster SVG is required.");
    await loadRenderer();
    const svg = svgElement.cloneNode(true);
    const hidden = document.createElement("div");
    hidden.setAttribute("aria-hidden","true");
    Object.assign(hidden.style,{position:"fixed",left:"-50000px",top:"0",pointerEvents:"none",zIndex:"-1"});
    // Use the actual A0 dimensions as the SVG viewport in the document unit (pt).
    // The original viewBox preserves all illustration, text and plot coordinates.
    svg.setAttribute("width",String(page.width));
    svg.setAttribute("height",String(page.height));
    svg.style.width=page.width+"px";
    svg.style.height=page.height+"px";
    svg.style.maxWidth="none";
    svg.style.maxHeight="none";
    hidden.append(svg);
    document.body.append(hidden);
    try {
      const fonts = await prepareFonts(svg);
      normalizeVectorFeatures(svg);
      const pdf = new window.jspdf.jsPDF({
        orientation:"portrait", unit:"pt", format:[page.width,page.height],
        putOnlyUsedFonts:true, compress:true, precision:5
      });
      for (const {file,base64,family,style} of fonts) {
        pdf.addFileToVFS(file,base64);
        pdf.addFont(file,family,style);
      }
      pdf.setProperties({
        title:svg.querySelector("title")?.textContent || "SATB2 Neuroday poster",
        subject:"SATB2 CA1 atlas poster, A0 portrait",
        creator:"SATB2 CA1 atlas poster editor"
      });
      await pdf.svg(svg,{x:0,y:0,width:page.width,height:page.height,loadExternalStyleSheets:false});
      addLinkAnnotations(pdf,svg);
      const blob = pdf.output("blob");
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href=url;
      link.download="SATB2_Neuroday_A0_edited.pdf";
      document.body.append(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url),60000);
      return blob;
    } finally {
      hidden.remove();
    }
  };
})();
