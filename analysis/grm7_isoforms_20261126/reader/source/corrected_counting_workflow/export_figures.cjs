const fs=require('node:fs'),path=require('node:path');
const {chromium}=require('/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
(async()=>{const b=await chromium.launch({headless:true,executablePath:'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',args:['--disable-gpu']});
 const dir=process.argv[2];const page=await b.newPage({deviceScaleFactor:2});
 for(const f of fs.readdirSync(dir).filter(f=>f.endsWith('.svg'))){const svg=fs.readFileSync(path.join(dir,f),'utf8');const [_,w,h]=svg.match(/width="(\d+)" height="(\d+)"/);await page.setViewportSize({width:+w,height:+h});await page.setContent('<style>body{margin:0}svg{display:block}</style>'+svg);await page.screenshot({path:path.join(dir,f.replace('.svg','.png')),fullPage:true});await page.pdf({path:path.join(dir,f.replace('.svg','.pdf')),width:w+'px',height:h+'px',margin:{top:0,bottom:0,left:0,right:0},printBackground:true});}
 await b.close();})().catch(e=>{console.error(e);process.exit(1)});
