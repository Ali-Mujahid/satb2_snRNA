'use strict';
// Display-only regression: trimming fractional zeros must never edit an exponent.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),crypto=require('crypto');
const site=process.env.SATB2_SITE_DIR||path.resolve(__dirname,'../site'),source=fs.readFileSync(path.join(site,'app.js'),'utf8');
const c={SATB2_NO_AUTO_INIT:true};c.window=c;c.globalThis=c;vm.createContext(c);vm.runInContext(source,c);const fmt=c.SATB2_APP.helpers.fmt;
const cases=[[1.234e-100,'1.234e-100'],[6.929e-110,'6.929e-110'],[1.828e-90,'1.828e-90'],[1e-100,'1e-100'],[1e-10,'1e-10'],[1e-20,'1e-20'],[1e-30,'1e-30'],[1e-300,'1e-300'],[0.01,'0.01'],[0.1,'0.1'],[0.0001,'0.0001'],[1,'1'],[10,'10'],[100,'100'],[0,'0'],[-1.234e-100,'-1.234e-100'],[-0.01,'-0.01'],[null,'NA'],[NaN,'NA'],[Infinity,'NA']];
for(const[value,expected]of cases){assert.equal(fmt(value),expected,'fmt('+String(value)+')');console.log('PASS '+String(value)+' -> '+expected);}
for(const exponent of[-10,-20,-30,-90,-100,-110,-200,-300]){const value=6.929*Math.pow(10,exponent),shown=fmt(value);assert.equal(shown.split('e')[1],value.toPrecision(4).split('e')[1],'Exponent must remain exact');}
assert.equal(fmt(1000),(1000).toLocaleString(undefined,{maximumFractionDigits:2}));
const report={status:'PASS',scope:'PURE_VM_FORMATTING_NO_DATA_MUTATION',example_cases:cases.length,exponent_cases:8,app_sha256:crypto.createHash('sha256').update(source).digest('hex')};fs.writeFileSync(path.join(__dirname,'number_format_acceptance.json'),JSON.stringify(report,null,2)+'\n');
