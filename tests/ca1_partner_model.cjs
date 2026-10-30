const assert = require('assert');
const M = require('../site/communication/ca1_partner/model.js');
const rows = [
  {passes_5pct:'True',passes_10pct:'False',majority_mouse_support_5pct:'True',majority_mouse_support_10pct:'False',direction:'CA1 to IN',ca1_population:'CA1_Deep',in_population:'IN_C1',in_status:'CA1 change only detected',ligand_complex:'Slit1',receptor_complex:'Robo1',ca1_significant_genes:'Slit1'},
  {passes_5pct:'True',passes_10pct:'True',majority_mouse_support_5pct:'False',majority_mouse_support_10pct:'True',direction:'IN to CA1',ca1_population:'CA1_Superficial',in_population:'IN_C6',in_status:'IN differential status unavailable',ligand_complex:'X',receptor_complex:'Y',ca1_significant_genes:'Y'},
  {passes_5pct:'False',passes_10pct:'False',majority_mouse_support_5pct:'False',majority_mouse_support_10pct:'False'}
];
const state={threshold:'5',direction:'',ca1:'',inpop:'',category:'',majority:false,query:''};
assert.equal(M.filter(rows,state).length,2);
assert.equal(M.filter(rows,{...state,threshold:'10'}).length,1);
assert.equal(M.filter(rows,{...state,majority:true}).length,1);
assert.equal(M.filter(rows,{...state,threshold:'10',majority:true}).length,1);
assert.equal(M.filter(rows,{...state,direction:'IN to CA1'}).length,1);
assert.equal(M.filter(rows,{...state,query:'robo1'}).length,1);
assert.equal(M.filter(rows,{...state,query:'none'}).length,0);
assert.equal(M.num(''),null);assert.equal(M.num('NA'),null);assert.equal(M.num('0'),0);
assert.equal(M.yes('False'),false);assert.equal(M.yes('True'),true);
console.log('model tests passed');
