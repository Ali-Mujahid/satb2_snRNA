// Source-only JavaScript/unit and local-link checks. No browser or network use.
const fs = require('fs'), path = require('path'), vm = require('vm'), crypto = require('crypto');
const guide = path.resolve(process.argv[2]), out = path.resolve(process.argv[3]);
if (fs.existsSync(out)) throw Error('Review output already exists');
const html = fs.readFileSync(guide, 'utf8'), script = html.match(/<script>([\s\S]*)<\/script>/)[1];
new Function(script);
const data = JSON.parse(script.match(/const DATA=([\s\S]*?);\nconst el=/)[1]);
class Element {
 constructor(tag) { this.tag=tag;this.children=[];this.style={};this.attributes={};this.value=''; }
 append(...nodes){this.children.push(...nodes)}
 replaceChildren(...nodes){this.children=[...nodes]}
 setAttribute(key,value){this.attributes[key]=value}
 addEventListener(){}
}
const elements = new Map();
const document={createElement:tag=>new Element(tag),getElementById:id=>{if(!elements.has(id))elements.set(id,new Element(id));return elements.get(id)}};
for (const [id,value] of Object.entries({population:'all',branch:'original_full',tool:'all',scope:'original_frozen'}))document.getElementById(id).value=value;
const context={document,window:{open:()=>{}},console};vm.createContext(context);vm.runInContext(script,context);
function assert(ok,message){if(!ok)throw Error(message)}
assert(elements.get('books').children.length===8,'Default original-full filter must show both tools/all four populations');
assert(elements.get('program-cards').children.length===4,'Default program cards differ');
document.getElementById('scope').value='fresh_branch';context.render();
assert(document.getElementById('branch').value==='paired_original','Fresh scope did not choose a valid paired branch');
assert(elements.get('books').children.length===8,'Fresh paired filter has wrong book coverage');
assert(elements.get('program-cards').children.length===1,'Frozen programs must be replaced with a scope notice in fresh view');
for (const [id,value] of Object.entries({population:'CA1_ProS',branch:'paired_cellbender',tool:'scanpy',scope:'original_frozen'}))document.getElementById(id).value=value;
context.render();
assert(elements.get('books').children.length===1&&elements.get('program-cards').children.length===1,'Single CA1 tool/branch filter differs');
const book=elements.get('books').children[0], image=book.children.find(x=>x.tag==='img'), tabs=book.children[3];
assert(tabs.children.length>=5,'Relevant CA1 marker/program preview tabs absent');
tabs.children[tabs.children.length-1].onclick();assert(image.src.startsWith('data:image/png;base64,'),'Thumbnail tab did not assign an embedded preview');
const links=[];
function visit(value,key=''){
 if(Array.isArray(value)){value.forEach(x=>visit(x,key));return}
 if(value&&typeof value==='object'){Object.entries(value).forEach(([k,x])=>visit(x,k));return}
 if(typeof value!=='string'||!['png','pdf','svg','preview_png','cluster_png','mouse_png','annotation_tsv','nucleus_tsv','crosswalk_book','review','semantic_table','overlap_table'].includes(key))return;
 assert(!/^[a-z]+:/i.test(value),'Guide artifact links must be relative local paths');
 const destination=path.resolve(path.dirname(guide),decodeURIComponent(value));assert(fs.statSync(destination).isFile(),'Artifact link is not a file');links.push(destination);
}
visit(data);
assert(data.books.length===40&&data.programs.length===12&&data.crosswalks.length===5,'Guide supplementary coverage differs');
fs.mkdirSync(out,{recursive:true});
const receipt={status:'PASS',scope:'SOURCE_LOGIC_AND_RELATIVE_FILE_LINKS_ONLY',guide_sha256:crypto.createHash('sha256').update(fs.readFileSync(guide)).digest('hex'),source_sha256:crypto.createHash('sha256').update(fs.readFileSync(__filename)).digest('hex'),n_book_datasets:data.books.length,n_program_maps:data.programs.length,n_crosswalks:data.crosswalks.length,n_local_link_references:links.length,n_unique_local_artifacts:new Set(links).size,checks:['embedded_JS_syntax','default_filter_8_books','fresh_scope_separation','single_population_tool_branch_filter','marker_thumbnail_tab_logic','relative_artifact_file_links'],actual_browser_interaction_QA:'UNVERIFIED_BROWSER_URL_SECURITY_POLICY',browser_or_network_used:false,scientific_fitting_performed:false};
fs.writeFileSync(path.join(out,'COMPLETE.json'),JSON.stringify(receipt,null,2)+'\n');console.log(JSON.stringify(receipt));
