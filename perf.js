const BD=require('./beamdyn2d.js'), BS=require('./beamstat2d.js'), TL=require('./trainloads.js');
function bench(name, inp){
  const t0=Date.now(); const m=BD.buildModel(inp); const tm=Date.now()-t0;
  const t1=Date.now(); const s=BS.run(inp); const ts=Date.now()-t1;
  const t2=Date.now();
  let n=0; for(const tr of inp.F_train){ BD.runTrain(m,inp,tr,inp.v,inp.v); n++; }
  const td=Date.now()-t2;
  console.log(name.padEnd(34), 'model',String(tm).padStart(5),'ms  stat',String(ts).padStart(6),
     'ms  dyn',String(td).padStart(6),'ms  ('+n+' trains, '+inp.v.length+' speeds)');
}
const v=[]; for(let k=100;k<=200;k++) v.push(k);
const base={L:[20,20],m:15e3,EI:30e9,xi:1.5,dL:1,bc:[[0,1],[20,1],[40,1]],v,W:0,Nmod:3,fmax:30,statCorr:true,secForces:true};
bench('D2 20 cars, dL=1',       Object.assign({},base,{F_train:TL.buildTrains('D2',20)}));
bench('HSLM-A x10, dL=1',       Object.assign({},base,{F_train:TL.buildTrains('HSLM-A',20)}));
bench('HSLM-A x10, dL=0.5',     Object.assign({},base,{dL:0.5,F_train:TL.buildTrains('HSLM-A',20)}));
bench('D2 20 cars, dL=0.25',    Object.assign({},base,{dL:0.25,F_train:TL.buildTrains('D2',20)}));
bench('D2, secForces off',      Object.assign({},base,{secForces:false,F_train:TL.buildTrains('D2',20)}));
