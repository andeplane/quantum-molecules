(function(){const t=document.createElement("link").relList;if(t&&t.supports&&t.supports("modulepreload"))return;for(const o of document.querySelectorAll('link[rel="modulepreload"]'))r(o);new MutationObserver(o=>{for(const a of o)if(a.type==="childList")for(const i of a.addedNodes)i.tagName==="LINK"&&i.rel==="modulepreload"&&r(i)}).observe(document,{childList:!0,subtree:!0});function n(o){const a={};return o.integrity&&(a.integrity=o.integrity),o.referrerPolicy&&(a.referrerPolicy=o.referrerPolicy),o.crossOrigin==="use-credentials"?a.credentials="include":o.crossOrigin==="anonymous"?a.credentials="omit":a.credentials="same-origin",a}function r(o){if(o.ep)return;o.ep=!0;const a=n(o);fetch(o.href,a)}})();const $=1/.52917721092,Z=27.211386245988,it=2.541746473,lt=[["H",1.008,.31,[.95,.95,.95]],["He",4.0026,.28,[.85,1,1]],["Li",6.94,1.28,[.8,.5,1]],["Be",9.0122,.96,[.76,1,0]],["B",10.81,.84,[1,.71,.71]],["C",12.011,.76,[.38,.38,.38]],["N",14.007,.71,[.19,.31,.97]],["O",15.999,.66,[1,.13,.13]],["F",18.998,.57,[.56,.88,.31]],["Ne",20.18,.58,[.7,.89,.96]],["Na",22.99,1.66,[.67,.36,.95]],["Mg",24.305,1.41,[.54,1,0]],["Al",26.982,1.21,[.75,.65,.65]],["Si",28.085,1.11,[.94,.78,.63]],["P",30.974,1.07,[1,.5,0]],["S",32.06,1.05,[1,.9,.19]],["Cl",35.45,1.02,[.12,.94,.12]],["Ar",39.948,1.06,[.5,.82,.89]]],ct=lt.map(([e,t,n,r],o)=>({symbol:e,z:o+1,mass:t,covalentRadius:n,color:r})),ut=new Map(ct.map(e=>[e.symbol.toLowerCase(),e]));function H(e){const t=ut.get(e.toLowerCase());if(!t)throw new Error(`Unsupported element: ${e} (only H to Ar are available)`);return t}const ie=(e,t)=>[e[0]+t[0],e[1]+t[1],e[2]+t[2]],ze=(e,t)=>[e[0]-t[0],e[1]-t[1],e[2]-t[2]],V=(e,t)=>[e[0]*t,e[1]*t,e[2]*t],N=(e,t)=>Math.hypot(e[0]-t[0],e[1]-t[1],e[2]-t[2]);function le(e,t){return ie(ie(V(e[0],t[0]),V(e[1],t[1])),V(e[2],t[2]))}function Ee(e,t,n){const r=[];for(let o=t;o<=n;o++)for(let a=t;a<=n;a++)for(let i=t;i<=n;i++)for(const l of e.sites){const u=[l.frac[0]+o,l.frac[1]+a,l.frac[2]+i];r.push({symbol:l.symbol,frac:u,xyz:le(e.cell,u)})}return r}function de(e,t=e.repeats){const n=le(e.cell,e.center),r=e.bond*1.15,o=Math.ceil((e.radius+2*e.bond)/Math.min(...e.cell.map(f=>Math.hypot(...f))))+1,a=Ee(e,-o-1,o+1);let i=a.filter(f=>N(f.xyz,n)<e.radius);for(;;){const f=i.filter(y=>i.filter(x=>x!==y&&N(y.xyz,x.xyz)<r).length>=2);if(f.length===i.length)break;i=f}const l=[];if(e.capBond!==void 0){for(const f of i)for(const y of a)if(N(f.xyz,y.xyz)<r&&N(f.xyz,y.xyz)>.1&&!i.includes(y)){const x=V(ze(y.xyz,f.xyz),1/N(f.xyz,y.xyz));l.push({symbol:"H",xyz:ie(f.xyz,V(x,e.capBond))})}}const u=f=>V(ze(f,n),$),h=[...i.map(f=>({symbol:f.symbol,xyz:u(f.xyz)})),...l.map(f=>({symbol:"H",xyz:u(f.xyz)}))],d=[0,1,2].map(f=>e.fixedAxes?.includes(f)?1:t),g=e.center.map((f,y)=>Math.floor(f-d[y]/2+.5)),b=1e-6,p=Ee(e,Math.min(...g)-1,Math.max(...g)+t+1).filter(f=>f.frac.every((y,x)=>y>=g[x]-b&&y<=g[x]+d[x]+b)),C=f=>i.some(y=>N(y.xyz,f)<.001),M=p.filter(f=>!C(f.xyz)),B=M.map(f=>({symbol:f.symbol,xyz:u(f.xyz)})),w=[],v=(f,y)=>{const x=N(f,y);return x<r&&x>.1};M.forEach((f,y)=>{for(let x=y+1;x<M.length;x++)v(f.xyz,M[x].xyz)&&w.push({ghost:y,other:x,otherIsCluster:!1});i.forEach((x,k)=>{v(f.xyz,x.xyz)&&w.push({ghost:y,other:k,otherIsCluster:!0})})});const S=[],O=(f,y,x)=>u(le(e.cell,[g[0]+f,g[1]+y,g[2]+x]));for(let f=0;f<=d[0];f++)for(let y=0;y<=d[1];y++)for(let x=0;x<=d[2];x++)f<d[0]&&S.push([O(f,y,x),O(f+1,y,x)]),y<d[1]&&S.push([O(f,y,x),O(f,y+1,x)]),x<d[2]&&S.push([O(f,y,x),O(f,y,x+1)]);return{cluster:h,caps:l.length,ghosts:B,ghostBonds:w,cellEdges:S}}const Pe=e=>[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],[.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]].map(t=>({symbol:e,frac:t})),Ae=(e,t)=>[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]].flatMap(n=>[{symbol:e,frac:n},{symbol:t,frac:[n[0]+.5,n[1],n[2]]}]),K=e=>[[e,0,0],[0,e,0],[0,0,e]],_e={diamond:{cell:K(3.567),sites:Pe("C"),center:[.5,.5,.5],radius:1.9,bond:1.5445,capBond:1.09,repeats:2,maxRepeats:3},silicon:{cell:K(5.431),sites:Pe("Si"),center:[.5,.5,.5],radius:2.9,bond:2.3517,capBond:1.48,repeats:2,maxRepeats:3},nacl:{cell:K(5.64),sites:Ae("Na","Cl"),center:[.25,.25,.25],radius:2.5,bond:2.82,repeats:1,maxRepeats:3},mgo:{cell:K(4.212),sites:Ae("Mg","O"),center:[.25,.25,.25],radius:1.9,bond:2.106,repeats:1,maxRepeats:4},graphite:{cell:[[2.464,0,0],[-1.232,2.1339,0],[0,0,6.711]],sites:[{symbol:"C",frac:[0,0,.25]},{symbol:"C",frac:[1/3,2/3,.25]},{symbol:"C",frac:[0,0,.75]},{symbol:"C",frac:[2/3,1/3,.75]}],center:[1/6,1/3,.25],radius:3.8,bond:1.4226,capBond:1.08,repeats:4,maxRepeats:6,fixedAxes:[2]}};function ht(e){return H(e.symbol).z}function fe(e){return e.atoms.reduce((t,n)=>t+ht(n),0)-e.charge}function dt(e,t=1.25){const n=[];for(let r=0;r<e.length;r++)for(let o=0;o<r;o++){const[a,i,l]=e[r].xyz,[u,h,d]=e[o].xyz;Math.hypot(a-u,i-h,l-d)/$<t*(H(e[r].symbol).covalentRadius+H(e[o].symbol).covalentRadius)&&n.push([o,r])}return n}function Ve(e){const t=[];for(const n of e.split(/\r?\n/)){const r=n.trim().split(/\s+/);if(r.length<4)continue;const[o,...a]=r,i=a.slice(0,3).map(Number);if(!/^[A-Za-z]{1,2}$/.test(o)||i.some(u=>!Number.isFinite(u)))continue;const l=o[0].toUpperCase()+o.slice(1).toLowerCase();H(l),t.push({symbol:l,xyz:i.map(u=>u*$)})}return t}function ft(e,t=""){const n=e.map(r=>`${r.symbol.padEnd(2)} ${r.xyz.map(o=>(o/$).toFixed(6).padStart(12)).join(" ")}`);return`${e.length}
${t}
${n.join(`
`)}
`}const Le={small:"Small molecules",large:"Larger molecules",crystal:"Pieces of crystals"},E=["sto-3g","6-31g","6-31g*"];function pt(){const e=[];for(let t=0;t<6;t++){const n=t*Math.PI/3;e.push(`C ${(1.397*Math.cos(n)).toFixed(4)} ${(1.397*Math.sin(n)).toFixed(4)} 0`)}for(let t=0;t<6;t++){const n=t*Math.PI/3;e.push(`H ${(2.481*Math.cos(n)).toFixed(4)} ${(2.481*Math.sin(n)).toFixed(4)} 0`)}return e.join(`
`)}function mt(e){return de(_e[e]).cluster.map(t=>`${t.symbol} ${t.xyz.map(n=>(n/$).toFixed(4)).join(" ")}`).join(`
`)}const F=["sto-3g"],gt=`N 1.2823 -0.5147 -0.5204
C 0.5109 0.5487 0.1775
C -0.9994 0.3600 -0.0969
O -1.7660 1.1859 -0.5529
O -1.4161 -0.9179 0.2628
H 2.2759 -0.3804 -0.2705
H 1.0177 -1.4145 -0.0875
H 0.6645 0.5519 1.2627
H 0.8096 1.5223 -0.2117
H -2.3793 -0.9414 0.0368`,yt=`C 0.0000 -1.1889 -0.0000
C 1.1902 -0.4804 0.0000
C 1.1443 0.9069 -0.0000
N -0.0000 1.6294 -0.0000
C -1.1443 0.9069 0.0000
C -1.1902 -0.4804 0.0000
H 0.0000 -2.2723 -0.0000
H 2.1430 -0.9925 0.0000
H 2.0665 1.4821 0.0000
H -2.0666 1.4820 0.0000
H -2.1430 -0.9926 -0.0000`,bt=`C -3.4689 -0.5199 -0.0085
C -2.0025 -0.4696 -0.4637
O -1.5908 -0.3051 -1.5926
O -1.1518 -0.6764 0.6362
C 0.2190 -0.5834 0.3187
C 0.9646 -1.7585 0.2673
C 2.3291 -1.7078 0.0362
C 2.9623 -0.4854 -0.1413
C 2.2263 0.6857 -0.0824
C 0.8535 0.6546 0.1471
C 0.0563 1.9438 0.1941
O -1.1555 2.0453 0.2931
O 0.8687 3.0708 0.1013
H -3.6626 -1.4583 0.5025
H -4.1210 -0.4313 -0.8712
H -3.6613 0.2977 0.6811
H 0.4589 -2.7051 0.4067
H 2.9011 -2.6268 -0.0081
H 4.0277 -0.4466 -0.3281
H 2.7107 1.6450 -0.2148
H 0.2362 3.8312 0.1263`,xt=`O -3.3581 -0.4226 0.8698
C -2.4381 0.1337 -0.0805
C -1.0469 0.2517 0.6038
O -0.5215 -1.0443 0.9703
C 0.3704 -1.5787 -0.0096
O -0.2886 -1.6403 -1.2813
C 1.7078 -0.7394 -0.0736
O 2.6939 -1.2430 0.8367
C 1.4009 0.7770 0.1695
O 2.2495 1.6487 -0.5896
C -0.0595 1.0996 -0.2580
O -0.3711 2.4826 -0.0486
H -4.2037 -0.5061 0.3595
H -2.3469 -0.4941 -0.9735
H -2.7445 1.1386 -0.3973
H -1.1910 0.7830 1.5517
H 0.6257 -2.5902 0.3491
H 0.2821 -2.2417 -1.8210
H 2.1345 -0.8470 -1.0765
H 2.2483 -1.2037 1.7211
H 1.4974 0.9816 1.2444
H 3.1598 1.4722 -0.2436
H -0.1715 0.8305 -1.3158
H 0.3712 2.9520 -0.5070`,vt=[{id:"h2",name:"Hydrogen H₂",xyz:`H 0 0 -0.3705
H 0 0 0.3705`,charge:0,multiplicity:1,bases:E},{id:"h2o",name:"Water H₂O",xyz:`O 0 0 0.1173
H 0 0.7572 -0.4692
H 0 -0.7572 -0.4692`,charge:0,multiplicity:1,bases:E},{id:"nh3",name:"Ammonia NH₃",xyz:`N 0 0 0.1162
H 0 0.9377 -0.2711
H 0.8121 -0.4689 -0.2711
H -0.8121 -0.4689 -0.2711`,charge:0,multiplicity:1,bases:E},{id:"ch4",name:"Methane CH₄",xyz:`C 0 0 0
H 0.6291 0.6291 0.6291
H -0.6291 -0.6291 0.6291
H -0.6291 0.6291 -0.6291
H 0.6291 -0.6291 -0.6291`,charge:0,multiplicity:1,bases:E},{id:"hf",name:"Hydrogen fluoride HF",xyz:`F 0 0 0.0917
H 0 0 -0.8253`,charge:0,multiplicity:1,bases:E},{id:"n2",name:"Nitrogen N₂",xyz:`N 0 0 -0.5488
N 0 0 0.5488`,charge:0,multiplicity:1,bases:E},{id:"co2",name:"Carbon dioxide CO₂",xyz:`C 0 0 0
O 0 0 1.1621
O 0 0 -1.1621`,charge:0,multiplicity:1,bases:E},{id:"ch2o",name:"Formaldehyde CH₂O",xyz:`C 0 0 -0.5285
O 0 0 0.6745
H 0 0.9363 -1.1079
H 0 -0.9363 -1.1079`,charge:0,multiplicity:1,bases:E},{id:"c2h4",name:"Ethylene C₂H₄",xyz:`C 0 0 0.6695
C 0 0 -0.6695
H 0 0.9288 1.2321
H 0 -0.9288 1.2321
H 0 0.9288 -1.2321
H 0 -0.9288 -1.2321`,charge:0,multiplicity:1,bases:E},{id:"c6h6",name:"Benzene C₆H₆",xyz:pt(),charge:0,multiplicity:1,bases:["sto-3g"]},{id:"oh",name:"Hydroxyl radical OH",xyz:`O 0 0 0
H 0 0 0.9697`,charge:0,multiplicity:2,bases:E},{id:"o2",name:"Oxygen O₂ (triplet)",xyz:`O 0 0 -0.6038
O 0 0 0.6038`,charge:0,multiplicity:3,bases:E}],wt=[{id:"glycine",name:"Glycine C₂H₅NO₂",xyz:gt,charge:0,multiplicity:1,bases:E},{id:"pyridine",name:"Pyridine C₅H₅N",xyz:yt,charge:0,multiplicity:1,bases:E},{id:"aspirin",name:"Aspirin C₉H₈O₄",xyz:bt,charge:0,multiplicity:1,bases:F},{id:"glucose",name:"Glucose C₆H₁₂O₆",xyz:xt,charge:0,multiplicity:1,bases:F}],Ct=[{id:"diamond",name:"Diamond (C₁₀H₁₆ piece)",crystal:"diamond",charge:0,multiplicity:1,bases:F},{id:"silicon",name:"Silicon (Si₁₀H₁₆ piece)",crystal:"silicon",charge:0,multiplicity:1,bases:F},{id:"graphite",name:"Graphite (C₁₆H₁₀ flake)",crystal:"graphite",charge:0,multiplicity:1,bases:F},{id:"nacl",name:"Table salt NaCl (Na₄Cl₄ cube)",crystal:"nacl",charge:0,multiplicity:1,bases:F},{id:"mgo",name:"Magnesium oxide MgO (Mg₄O₄ cube)",crystal:"mgo",charge:0,multiplicity:1,bases:F}],pe=[...vt.map(e=>({...e,group:"small"})),...wt.map(e=>({...e,group:"large"})),...Ct.map(e=>({...e,group:"crystal",xyz:mt(e.crystal)}))];function Mt(e){const t=pe.find(n=>n.id===e);if(!t)throw new Error(`Unknown preset ${e}`);return t}const Bt={h2:"The simplest molecule: two hydrogen atoms sharing their two electrons in one bond.",h2o:"Oxygen pulls the shared electrons toward itself and leaves the hydrogens slightly positive. That imbalance is what makes water such a good solvent.",nh3:"Nitrogen bonds three hydrogens and keeps one pair of electrons to itself, which pushes the molecule into a pyramid.",ch4:"Carbon shares its four outer electrons with four hydrogens, which point to the corners of a tetrahedron.",hf:"Fluorine pulls on electrons harder than any other element, so this bond is very lopsided.",n2:"A triple bond, one of the strongest in chemistry. It is why the nitrogen that makes up most of the air hardly reacts.",co2:"Two double bonds in a straight line. Each bond is lopsided, but they point in opposite directions and cancel out.",ch2o:"A carbon–oxygen double bond, the group at the heart of sugars and many other molecules of life.",c2h4:"A carbon–carbon double bond. Show the HOMO orbital to see the extra bond lying above and below the molecule.",c6h6:"Six electrons shared around the whole ring instead of sitting in fixed bonds. It uses the quickest settings to stay smooth.",oh:"A radical: one electron has no partner, which makes it extremely reactive.",o2:"The oxygen we breathe has two unpaired electrons, which is why liquid oxygen sticks to a magnet.",glycine:"The simplest amino acid, one of the twenty building blocks of every protein: an amine group and an acid group on one carbon.",pyridine:"Benzene with one carbon swapped for a nitrogen. The nitrogen keeps a pair of electrons pointing out of the ring, which makes pyridine a base.",aspirin:"Twenty-one atoms: a benzene ring carrying an acid group and an ester. Its highest filled orbitals are the ring’s π electrons.",glucose:"The sugar our cells burn for energy, in the six-membered ring it forms in water. Twenty-four atoms, with five OH groups around the ring.",diamond:"Every carbon bonds to four others in a tetrahedron, without end. The ten carbons in the middle of one cubic cell are calculated; the bonds that would carry on into the crystal end in hydrogen atoms.",silicon:"The same structure as diamond with bigger atoms and longer, weaker bonds. Compare the gap between the highest filled and lowest empty levels with diamond’s: smaller, which is why silicon is a semiconductor.",graphite:"Carbon in flat sheets of hexagons, stacked loosely: pencil lead. The calculated flake has its π electrons spread over the whole sheet; show the HOMO to see them above and below the plane.",nacl:"Sodium hands an electron to chlorine and the ions stack in a checkerboard of cubes. Eight ions are calculated; the cloud piles up on the chlorines, which hold the extra electrons.",mgo:"The same checkerboard as table salt, with magnesium giving two electrons to each oxygen. The ions are smaller and more charged, so they pack tighter."};function De(e,t="alpha"){const n=t==="beta"&&e.unrestricted?e.epsB:e.epsA,r=t==="beta"?e.nbeta:e.nalpha;return n.map((o,a)=>{let i="";return a===r-1?i="HOMO":a===r?i="LUMO":a<r?i=`HOMO−${r-1-a}`:i=`LUMO+${a-r}`,{index:a,energy:o,occupied:a<r,label:i}})}function Ot(e){let t=e.epsA[e.nalpha-1]??-1/0,n=e.epsA[e.nalpha]??1/0;return e.unrestricted&&(e.nbeta>0&&(t=Math.max(t,e.epsB[e.nbeta-1])),n=Math.min(n,e.epsB[e.nbeta]??1/0)),{homo:t,lumo:Number.isFinite(n)?n:null}}class St{worker;nextId=1;pending=new Map;relaxes=new Map;busy=!1;queued=null;ready;constructor(t){this.worker=t;let n;this.ready=new Promise(r=>n=r),t.onmessage=r=>{const o=r.data;if(o.type==="ready")return n();if(o.type==="relax-step"||o.type==="relax-done"){const i=this.relaxes.get(o.id);if(!i)return;o.type==="relax-step"?i.onStep(o.output,o.step):(this.relaxes.delete(o.id),i.onDone(o.converged,o.steps));return}if(o.type==="error"&&this.relaxes.has(o.id)){const i=this.relaxes.get(o.id);this.relaxes.delete(o.id),i.onError(o.message);return}const a=this.pending.get(o.id);a&&(this.pending.delete(o.id),o.type==="result"?a.resolve(o.output):o.type==="error"&&a.reject(new Error(o.message)))}}send(t){this.worker.postMessage(t)}compute(t){const n=this.nextId++;return new Promise((r,o)=>{this.pending.set(n,{resolve:r,reject:o}),this.send({id:n,type:"compute",input:t})})}computeLatest(t){return this.busy?(this.queued?.resolve(null),new Promise(n=>this.queued={input:t,resolve:n})):(this.busy=!0,this.compute(t).finally(()=>{this.busy=!1;const n=this.queued;this.queued=null,n&&this.computeLatest(n.input).then(n.resolve,()=>n.resolve(null))}))}relax(t,n,r=60){const o=this.nextId++;return this.relaxes.set(o,n),this.send({id:o,type:"relax",input:t,maxSteps:r}),o}cancelRelax(t){this.relaxes.delete(t),this.send({id:t,type:"cancel"})}reset(){this.send({id:this.nextId++,type:"reset"})}configure(t){this.send({id:this.nextId++,type:"configure",threads:t})}}const zt=8;function Et(e){const t=Number(e.get("threads"))||navigator.hardwareConcurrency||4;return Math.max(1,Math.min(16,Math.round(t)))}const J=(e,t)=>[e[0]-t[0],e[1]-t[1],e[2]-t[2]],qe=(e,t)=>[e[0]+t[0],e[1]+t[1],e[2]+t[2]],me=(e,t)=>[e[0]*t,e[1]*t,e[2]*t],G=(e,t)=>e[0]*t[0]+e[1]*t[1]+e[2]*t[2],He=(e,t)=>[e[1]*t[2]-e[2]*t[1],e[2]*t[0]-e[0]*t[2],e[0]*t[1]-e[1]*t[0]],Pt=e=>Math.hypot(e[0],e[1],e[2]),ce=e=>me(e,1/(Pt(e)||1));function At(e,t){const n=new Float32Array(16);for(let r=0;r<4;r++)for(let o=0;o<4;o++){let a=0;for(let i=0;i<4;i++)a+=e[i*4+o]*t[r*4+i];n[r*4+o]=a}return n}function Lt(e,t,n,r){const o=1/Math.tan(e/2),a=new Float32Array(16);return a[0]=o/t,a[5]=o,a[10]=r/(n-r),a[11]=-1,a[14]=r*n/(n-r),a}function Ht(e,t,n){const r=ce(J(e,t)),o=ce(He(n,r)),a=He(r,o);return new Float32Array([o[0],a[0],r[0],0,o[1],a[1],r[1],0,o[2],a[2],r[2],0,-G(o,e),-G(a,e),-G(r,e),1])}function je(e){const t=new Float64Array(16),n=Array.from(e),r=Array.from({length:4},(o,a)=>[0,1,2,3].map(i=>n[i*4+a]).concat([0,1,2,3].map(i=>i===a?1:0)));for(let o=0;o<4;o++){let a=o;for(let l=o+1;l<4;l++)Math.abs(r[l][o])>Math.abs(r[a][o])&&(a=l);[r[o],r[a]]=[r[a],r[o]];const i=r[o][o];if(Math.abs(i)<1e-30)throw new Error("singular matrix");for(let l=0;l<8;l++)r[o][l]/=i;for(let l=0;l<4;l++){if(l===o)continue;const u=r[l][o];for(let h=0;h<8;h++)r[l][h]-=u*r[o][h]}}for(let o=0;o<4;o++)for(let a=0;a<4;a++)t[a*4+o]=r[o][4+a];return new Float32Array(t)}function Te(e,t){const n=e[0]*t[0]+e[4]*t[1]+e[8]*t[2]+e[12],r=e[1]*t[0]+e[5]*t[1]+e[9]*t[2]+e[13],o=e[2]*t[0]+e[6]*t[1]+e[10]*t[2]+e[14],a=e[3]*t[0]+e[7]*t[1]+e[11]*t[2]+e[15];return[n/a,r/a,o/a]}function Tt(e,t,n,r){const o=J(e,n),a=G(o,t),i=G(o,o)-r*r,l=a*a-i;if(l<0)return null;const u=-a-Math.sqrt(l);return u>0?u:null}class Rt{target=[0,0,0];distance=15;yaw=.6;pitch=.35;fovy=40*Math.PI/180;aspect=1;get eye(){const t=Math.cos(this.pitch);return qe(this.target,me([t*Math.sin(this.yaw),Math.sin(this.pitch),t*Math.cos(this.yaw)],this.distance))}get view(){return Ht(this.eye,this.target,[0,1,0])}get projection(){return Lt(this.fovy,this.aspect,Math.max(.05,this.distance*.02),this.distance*10+50)}get viewProjection(){return At(this.projection,this.view)}orbit(t,n){this.yaw-=t*.008,this.pitch=Math.max(-1.5,Math.min(1.5,this.pitch+n*.008))}zoom(t){this.distance=Math.max(2,Math.min(200,this.distance*t))}fit(t,n){this.target=t,this.distance=Math.max(6,(n+2.5)/Math.sin(this.fovy/2))}ray(t,n){const r=je(this.viewProjection),o=Te(r,[t,n,0]),a=Te(r,[t,n,1]);return{origin:o,dir:ce(J(a,o))}}basis(){const t=this.view;return{right:[t[0],t[4],t[8]],up:[t[1],t[5],t[9]],forward:[-t[2],-t[6],-t[10]]}}}const kt=`// Ball-and-stick atoms: spheres are ray-cast impostors, bonds are instanced cylinders.

struct Camera {
  view: mat4x4f,
  proj: mat4x4f,
}

@group(0) @binding(0) var<uniform> cam: Camera;

// Ghosts: the crystal around a computed cluster, drawn matte and see-through (premultiplied
// alpha, no depth writes) so the calculated atoms and their cloud stay visible behind them.
const GHOST_ALPHA = 0.5;

fn ghostShade(n: vec3f, base: vec3f) -> vec3f {
  let pale = mix(base, vec3f(0.8, 0.83, 0.92), 0.4);
  let diff = max(dot(n, normalize(vec3f(0.4, 0.6, 1.0))), 0.0);
  return pale * (0.3 + 0.6 * diff);
}

struct Sphere {
  @location(0) center: vec3f,
  @location(1) radius: f32,
  @location(2) color: vec3f,
  @location(3) highlight: f32,
}

struct SphereOut {
  @builtin(position) pos: vec4f,
  @location(0) viewPos: vec3f,
  @location(1) center: vec3f,
  @location(2) radius: f32,
  @location(3) color: vec3f,
  @location(4) highlight: f32,
}

@vertex
fn sphereVs(@builtin(vertex_index) vi: u32, s: Sphere) -> SphereOut {
  var corners = array<vec2f, 4>(vec2f(-1.0, -1.0), vec2f(1.0, -1.0), vec2f(-1.0, 1.0), vec2f(1.0, 1.0));
  let c = (cam.view * vec4f(s.center, 1.0)).xyz;
  // Grow the quad and push it towards the camera so the whole sphere is covered.
  let grow = 1.6;
  let vp = c + vec3f(corners[vi] * s.radius * grow, s.radius);
  var o: SphereOut;
  o.pos = cam.proj * vec4f(vp, 1.0);
  o.viewPos = vp;
  o.center = c;
  o.radius = s.radius;
  o.color = s.color;
  o.highlight = s.highlight;
  return o;
}

fn shade(n: vec3f, v: vec3f, base: vec3f) -> vec3f {
  let l = normalize(vec3f(0.4, 0.6, 1.0));
  let diff = max(dot(n, l), 0.0);
  let spec = pow(max(dot(n, normalize(l + v)), 0.0), 48.0);
  let rim = pow(1.0 - max(dot(n, v), 0.0), 3.0);
  return base * (0.22 + 0.78 * diff) + vec3f(0.35 * spec) + 0.15 * rim * base;
}

struct FragOut {
  @location(0) color: vec4f,
  @builtin(frag_depth) depth: f32,
}

@fragment
fn sphereFs(in: SphereOut) -> FragOut {
  let dir = normalize(in.viewPos);
  let b = dot(dir, in.center);
  let disc = b * b - dot(in.center, in.center) + in.radius * in.radius;
  if (disc < 0.0) {
    discard;
  }
  let t = b - sqrt(disc);
  let hit = dir * t;
  let n = normalize(hit - in.center);
  let clip = cam.proj * vec4f(hit, 1.0);
  var base = in.color;
  if (in.highlight > 0.5) {
    base = mix(base, vec3f(1.0, 0.85, 0.3), 0.45);
  }
  var o: FragOut;
  o.color = vec4f(shade(n, -dir, base), 1.0);
  o.depth = clip.z / clip.w;
  return o;
}

struct Bond {
  @location(0) a: vec3f,
  @location(1) radius: f32,
  @location(2) b: vec3f,
  @location(3) colorA: vec3f,
  @location(4) colorB: vec3f,
}

struct BondOut {
  @builtin(position) pos: vec4f,
  @location(0) normal: vec3f,
  @location(1) viewPos: vec3f,
  @location(2) along: f32,
  @location(3) colorA: vec3f,
  @location(4) colorB: vec3f,
}

// \`ring\` is (cos, sin, along) on a unit cylinder from along=0 to along=1.
@vertex
fn bondVs(@location(5) ring: vec3f, bond: Bond) -> BondOut {
  let axis = bond.b - bond.a;
  let w = normalize(axis);
  let helper = select(vec3f(1.0, 0.0, 0.0), vec3f(0.0, 1.0, 0.0), abs(w.x) > 0.9);
  let u = normalize(cross(w, helper));
  let v = cross(w, u);
  let n = u * ring.x + v * ring.y;
  let p = bond.a + axis * ring.z + n * bond.radius;
  var o: BondOut;
  let vp = cam.view * vec4f(p, 1.0);
  o.pos = cam.proj * vp;
  o.viewPos = vp.xyz;
  o.normal = (cam.view * vec4f(n, 0.0)).xyz;
  o.along = ring.z;
  o.colorA = bond.colorA;
  o.colorB = bond.colorB;
  return o;
}

@fragment
fn bondFs(in: BondOut) -> @location(0) vec4f {
  let base = select(in.colorB, in.colorA, in.along < 0.5);
  return vec4f(shade(normalize(in.normal), normalize(-in.viewPos), base * 0.9), 1.0);
}

@fragment
fn sphereGhostFs(in: SphereOut) -> FragOut {
  let dir = normalize(in.viewPos);
  let b = dot(dir, in.center);
  let disc = b * b - dot(in.center, in.center) + in.radius * in.radius;
  if (disc < 0.0) {
    discard;
  }
  let t = b - sqrt(disc);
  let hit = dir * t;
  let n = normalize(hit - in.center);
  let clip = cam.proj * vec4f(hit, 1.0);
  // Softer towards the rim, like a glass bead.
  let a = GHOST_ALPHA * (0.55 + 0.45 * max(dot(n, -dir), 0.0));
  var o: FragOut;
  o.color = vec4f(ghostShade(n, in.color) * a, a);
  o.depth = clip.z / clip.w;
  return o;
}

@fragment
fn bondGhostFs(in: BondOut) -> @location(0) vec4f {
  let base = select(in.colorB, in.colorA, in.along < 0.5);
  let a = GHOST_ALPHA * 0.8;
  return vec4f(ghostShade(normalize(in.normal), base) * a, a);
}
`,Nt=`// Ray-marches the density cloud and the orbital isosurfaces, stopping at the atoms' depth.

struct Uniforms {
  invViewProj: mat4x4f,
  eye: vec3f,
  cloudOpacity: f32,
  boxMin: vec3f,
  logMin: f32,
  boxMax: vec3f,
  logMax: f32,
  voxel: vec3f,
  iso: f32,
  showCloud: f32,
  showOrbital: f32,
  orbitalOpacity: f32,
  frame: f32,
}

@group(0) @binding(0) var<uniform> u: Uniforms;
@group(0) @binding(1) var volume: texture_3d<f32>;
@group(0) @binding(2) var samp: sampler;
@group(0) @binding(3) var depthTex: texture_depth_2d;

struct VsOut {
  @builtin(position) pos: vec4f,
  @location(0) ndc: vec2f,
}

@vertex
fn vs(@builtin(vertex_index) i: u32) -> VsOut {
  let xy = vec2f(f32((i << 1u) & 2u), f32(i & 2u)) * 2.0 - 1.0;
  var o: VsOut;
  o.pos = vec4f(xy, 0.0, 1.0);
  o.ndc = xy;
  return o;
}

fn unproject(ndc: vec2f, z: f32) -> vec3f {
  let p = u.invViewProj * vec4f(ndc, z, 1.0);
  return p.xyz / p.w;
}

fn toUvw(p: vec3f) -> vec3f {
  return (p - u.boxMin) / (u.boxMax - u.boxMin);
}

fn sampleVol(p: vec3f) -> vec4f {
  return textureSampleLevel(volume, samp, toUvw(p), 0.0);
}

// A muted violet-to-white ramp for the cloud, kept apart from the blue and orange orbital lobes.
// Keep in sync with CLOUD_RAMP in src/ui/legend.ts.
fn cloudColor(t: f32) -> vec3f {
  let c0 = vec3f(0.30, 0.28, 0.40);
  let c1 = vec3f(0.66, 0.64, 0.78);
  let c2 = vec3f(1.00, 0.98, 0.95);
  return select(mix(c1, c2, (t - 0.5) * 2.0), mix(c0, c1, t * 2.0), t < 0.5);
}

fn hash(p: vec2f) -> f32 {
  return fract(sin(dot(p, vec2f(12.9898, 78.233)) + u.frame * 0.618) * 43758.5453);
}

fn orbitalNormal(p: vec3f) -> vec3f {
  let h = u.voxel;
  let g = vec3f(
    sampleVol(p + vec3f(h.x, 0.0, 0.0)).g - sampleVol(p - vec3f(h.x, 0.0, 0.0)).g,
    sampleVol(p + vec3f(0.0, h.y, 0.0)).g - sampleVol(p - vec3f(0.0, h.y, 0.0)).g,
    sampleVol(p + vec3f(0.0, 0.0, h.z)).g - sampleVol(p - vec3f(0.0, 0.0, h.z)).g,
  ) / (2.0 * h);
  return normalize(g + vec3f(1e-12));
}

@fragment
fn fs(in: VsOut) -> @location(0) vec4f {
  let frag = in.pos;
  let far = unproject(in.ndc, 1.0);
  let dir = normalize(far - u.eye);
  var tmax = 1e9;
  let depth = textureLoad(depthTex, vec2i(frag.xy), 0);
  if (depth < 1.0) {
    tmax = distance(unproject(in.ndc, depth), u.eye);
  }
  let inv = 1.0 / dir;
  let ta = (u.boxMin - u.eye) * inv;
  let tb = (u.boxMax - u.eye) * inv;
  let t0 = max(max(max(min(ta.x, tb.x), min(ta.y, tb.y)), min(ta.z, tb.z)), 0.0);
  let t1 = min(min(min(max(ta.x, tb.x), max(ta.y, tb.y)), max(ta.z, tb.z)), tmax);
  if (t1 <= t0) {
    discard;
  }
  let dt = 0.6 * min(min(u.voxel.x, u.voxel.y), u.voxel.z);
  let steps = min(u32((t1 - t0) / dt) + 1u, 1024u);
  var t = t0 + dt * hash(frag.xy);
  var color = vec3f(0.0);
  var alpha = 0.0;
  var prevPsi = sampleVol(u.eye + dir * t).g;
  let light = normalize(-dir + vec3f(0.3, 0.5, 0.0));
  for (var i = 0u; i < steps; i++) {
    if (t > t1 || alpha > 0.985) {
      break;
    }
    let p = u.eye + dir * t;
    let s = sampleVol(p);
    if (u.showOrbital > 0.5) {
      let psi = s.g;
      // A crossing of +iso or -iso between the previous sample and this one.
      let crossPos = (prevPsi - u.iso) * (psi - u.iso) < 0.0;
      let crossNeg = (prevPsi + u.iso) * (psi + u.iso) < 0.0;
      if (crossPos || crossNeg) {
        let level = select(-u.iso, u.iso, crossPos);
        let f = clamp((level - prevPsi) / (psi - prevPsi), 0.0, 1.0);
        let hit = p - dir * dt * (1.0 - f);
        var n = orbitalNormal(hit);
        if (dot(n, dir) > 0.0) {
          n = -n;
        }
        // Keep in sync with LOBE_POSITIVE and LOBE_NEGATIVE in src/ui/legend.ts.
        let base = select(vec3f(0.95, 0.42, 0.18), vec3f(0.22, 0.52, 1.0), level > 0.0);
        let diff = max(dot(n, light), 0.0);
        let spec = pow(max(dot(n, normalize(light - dir)), 0.0), 40.0);
        let shaded = base * (0.25 + 0.75 * diff) + vec3f(0.35 * spec);
        let a = u.orbitalOpacity;
        color += (1.0 - alpha) * a * shaded;
        alpha += (1.0 - alpha) * a;
      }
      prevPsi = psi;
    }
    if (u.showCloud > 0.5) {
      let rho = max(s.r, 1e-12);
      let x = clamp((log2(rho) * 0.30103 - u.logMin) / (u.logMax - u.logMin), 0.0, 1.0);
      let a = 1.0 - exp(-u.cloudOpacity * x * x * x * dt);
      color += (1.0 - alpha) * a * cloudColor(x) * (0.8 + 0.8 * x);
      alpha += (1.0 - alpha) * a;
    }
    t += dt;
  }
  return vec4f(color, alpha);
}
`,Ut=`// Evaluates the electron density and one molecular orbital on a 3D grid.
// The basis layout matches \`gpu_basis\` in crates/qchem-wasm and \`evalAo\` in src/chem/orbitals.ts.

const MAX_BF: u32 = 128u;

struct Params {
  boxMin: vec3f,
  nbf: u32,
  boxMax: vec3f,
  nshell: u32,
  dims: vec3u,
  // Rows of \`coefs\`: \`nocc\` occupied orbitals weighted by \`weights\`, then the selected orbital.
  nocc: u32,
}

struct Shell {
  center: vec3f,
  l: u32,
  nprim: u32,
  prim: u32,
  fn0: u32,
  // Smallest exponent, for skipping shells that vanish at this point.
  minExp: f32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> shells: array<Shell>;
@group(0) @binding(2) var<storage, read> prims: array<vec2f>;
// Cartesian component norms: [s, px, py, pz, dxx, dxy, dxz, dyy, dyz, dzz].
@group(0) @binding(3) var<storage, read> norms: array<f32, 10>;
@group(0) @binding(4) var<storage, read> coefs: array<f32>;
@group(0) @binding(5) var<storage, read> weights: array<f32>;
@group(0) @binding(6) var volume: texture_storage_3d<rgba16float, write>;

var<private> ao: array<f32, MAX_BF>;

fn evalAo(p: vec3f) {
  for (var s = 0u; s < params.nshell; s++) {
    let sh = shells[s];
    let d = p - sh.center;
    let r2 = dot(d, d);
    let f = sh.fn0;
    var radial = 0.0;
    if (sh.minExp * r2 < 40.0) {
      for (var k = 0u; k < sh.nprim; k++) {
        let pr = prims[sh.prim + k];
        radial += pr.y * exp(-pr.x * r2);
      }
    }
    if (sh.l == 0u) {
      ao[f] = radial * norms[0];
    } else if (sh.l == 1u) {
      ao[f] = radial * d.x * norms[1];
      ao[f + 1u] = radial * d.y * norms[2];
      ao[f + 2u] = radial * d.z * norms[3];
    } else {
      ao[f] = radial * d.x * d.x * norms[4];
      ao[f + 1u] = radial * d.x * d.y * norms[5];
      ao[f + 2u] = radial * d.x * d.z * norms[6];
      ao[f + 3u] = radial * d.y * d.y * norms[7];
      ao[f + 4u] = radial * d.y * d.z * norms[8];
      ao[f + 5u] = radial * d.z * d.z * norms[9];
    }
  }
}

fn moValue(row: u32) -> f32 {
  var v = 0.0;
  let base = row * params.nbf;
  for (var i = 0u; i < params.nbf; i++) {
    v += ao[i] * coefs[base + i];
  }
  return v;
}

@compute @workgroup_size(4, 4, 4)
fn main(@builtin(global_invocation_id) id: vec3u) {
  if (any(id >= params.dims)) {
    return;
  }
  let t = vec3f(id) / vec3f(max(params.dims - 1u, vec3u(1u)));
  let p = mix(params.boxMin, params.boxMax, t);
  evalAo(p);
  var rho = 0.0;
  for (var k = 0u; k < params.nocc; k++) {
    let v = moValue(k);
    rho += weights[k] * v * v;
  }
  let psi = moValue(params.nocc);
  textureStore(volume, id, vec4f(rho, psi, 0.0, 1.0));
}
`;function Re(e){const t=new Float32Array(Math.max(1,e.length)*8);return e.forEach((n,r)=>t.set([...n.center,n.radius,...n.color,n.highlight?1:0],8*r)),t}function ke(e){const t=new Float32Array(Math.max(1,e.length)*13);return e.forEach((n,r)=>t.set([...n.a,n.radius,...n.b,...n.colorA,...n.colorB],13*r)),t}const re="depth24plus",se=18;function Ft(){const e=[];for(let t=0;t<se;t++){const n=t/se*2*Math.PI,r=(t+1)/se*2*Math.PI,o=(a,i)=>e.push(Math.cos(a),Math.sin(a),i);o(n,0),o(r,0),o(r,1),o(n,0),o(r,1),o(n,1)}return new Float32Array(e)}class ge{device;context;format;depth=null;cameraBuffer;rayBuffer;paramsBuffer;normsBuffer;ringBuffer;ringCount;sphereBuffer=null;sphereCount=0;bondBuffer=null;bondCount=0;ghostSphereBuffer=null;ghostSphereCount=0;ghostBondBuffer=null;ghostBondCount=0;spherePipeline;bondPipeline;ghostSpherePipeline;ghostBondPipeline;rayPipeline;volumePipeline;atomsBindGroup;sampler;volume=null;box=null;frame=0;adapter=null;background={r:.035,g:.045,b:.07,a:1};static async create(t){if(!navigator.gpu)throw new Error("WebGPU is not available in this browser.");const n=await navigator.gpu.requestAdapter({powerPreference:"high-performance"});if(!n)throw new Error("No WebGPU adapter was found.");const r=await n.requestDevice(),o=t?t.getContext("webgpu"):null;if(t&&!o)throw new Error("Could not create a WebGPU canvas context.");const a=new ge(r,o,navigator.gpu.getPreferredCanvasFormat());return a.adapter=n,a}constructor(t,n,r){this.device=t,this.context=n,this.format=r,n?.configure({device:t,format:r,alphaMode:"opaque"});const o=GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST;this.cameraBuffer=t.createBuffer({size:128,usage:o}),this.rayBuffer=t.createBuffer({size:144,usage:o}),this.paramsBuffer=t.createBuffer({size:48,usage:o}),this.normsBuffer=t.createBuffer({size:40,usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST});const a=Ft();this.ringCount=a.length/3,this.ringBuffer=t.createBuffer({size:a.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST}),t.queue.writeBuffer(this.ringBuffer,0,a),this.sampler=t.createSampler({magFilter:"linear",minFilter:"linear"});const i=t.createShaderModule({code:kt,label:"atoms"}),l=t.createPipelineLayout({bindGroupLayouts:[this.cameraLayout()]}),u={format:re,depthWriteEnabled:!0,depthCompare:"less"},h=(B,w)=>({shaderLocation:B,offset:w,format:"float32x3"}),d=(B,w)=>({shaderLocation:B,offset:w,format:"float32"}),g={format:r,blend:{color:{srcFactor:"one",dstFactor:"one-minus-src-alpha"},alpha:{srcFactor:"one",dstFactor:"one-minus-src-alpha"}}},b={format:re,depthWriteEnabled:!1,depthCompare:"less"},p=B=>t.createRenderPipeline({label:B?"ghost spheres":"spheres",layout:l,vertex:{module:i,entryPoint:"sphereVs",buffers:[{arrayStride:32,stepMode:"instance",attributes:[h(0,0),d(1,12),h(2,16),d(3,28)]}]},fragment:{module:i,entryPoint:B?"sphereGhostFs":"sphereFs",targets:[B?g:{format:r}]},primitive:{topology:"triangle-strip"},depthStencil:B?b:u}),C=B=>t.createRenderPipeline({label:B?"ghost bonds":"bonds",layout:l,vertex:{module:i,entryPoint:"bondVs",buffers:[{arrayStride:52,stepMode:"instance",attributes:[h(0,0),d(1,12),h(2,16),h(3,28),h(4,40)]},{arrayStride:12,stepMode:"vertex",attributes:[h(5,0)]}]},fragment:{module:i,entryPoint:B?"bondGhostFs":"bondFs",targets:[B?g:{format:r}]},primitive:{topology:"triangle-list",cullMode:"back"},depthStencil:B?b:u});this.spherePipeline=p(!1),this.bondPipeline=C(!1),this.ghostSpherePipeline=p(!0),this.ghostBondPipeline=C(!0),this.atomsBindGroup=t.createBindGroup({layout:this.cameraLayout(),entries:[{binding:0,resource:{buffer:this.cameraBuffer}}]});const M=t.createShaderModule({code:Nt,label:"raymarch"});this.rayPipeline=t.createRenderPipeline({label:"raymarch",layout:"auto",vertex:{module:M,entryPoint:"vs"},fragment:{module:M,entryPoint:"fs",targets:[{format:r,blend:{color:{srcFactor:"one",dstFactor:"one-minus-src-alpha"},alpha:{srcFactor:"one",dstFactor:"one-minus-src-alpha"}}}]}}),this.volumePipeline=t.createComputePipeline({label:"volume",layout:"auto",compute:{module:t.createShaderModule({code:Ut,label:"volume"}),entryPoint:"main"}})}cameraLayoutCache=null;cameraLayout(){return this.cameraLayoutCache??=this.device.createBindGroupLayout({entries:[{binding:0,visibility:GPUShaderStage.VERTEX|GPUShaderStage.FRAGMENT,buffer:{type:"uniform"}}]}),this.cameraLayoutCache}resize(t,n){this.depth&&this.depth.width===t&&this.depth.height===n||(this.depth?.destroy(),this.depth=this.device.createTexture({size:[t,n],format:re,usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.TEXTURE_BINDING}))}setAtoms(t,n){const r=t.filter(l=>!l.ghost),o=t.filter(l=>l.ghost);this.sphereBuffer=this.upload(this.sphereBuffer,Re(r)),this.sphereCount=r.length,this.ghostSphereBuffer=this.upload(this.ghostSphereBuffer,Re(o)),this.ghostSphereCount=o.length;const a=n.filter(l=>!l.ghost),i=n.filter(l=>l.ghost);this.bondBuffer=this.upload(this.bondBuffer,ke(a)),this.bondCount=a.length,this.ghostBondBuffer=this.upload(this.ghostBondBuffer,ke(i)),this.ghostBondCount=i.length}upload(t,n){return(!t||t.size<n.byteLength)&&(t?.destroy(),t=this.device.createBuffer({size:n.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST})),this.device.queue.writeBuffer(t,0,n),t}updateVolume(t){const{device:n}=this,[r,o,a]=t.box.dims;(!this.volume||this.volume.width!==r||this.volume.height!==o||this.volume.depthOrArrayLayers!==a)&&(this.volume?.destroy(),this.volume=n.createTexture({size:[r,o,a],dimension:"3d",format:"rgba16float",usage:GPUTextureUsage.STORAGE_BINDING|GPUTextureUsage.TEXTURE_BINDING|GPUTextureUsage.COPY_SRC})),this.box=t.box;const i=new ArrayBuffer(48),l=new Float32Array(i),u=new Uint32Array(i);l.set(t.box.min,0),u[3]=t.nbf,l.set(t.box.max,4),u[7]=t.nshell,u.set(t.box.dims,8),u[11]=t.nocc,n.queue.writeBuffer(this.paramsBuffer,0,i),n.queue.writeBuffer(this.normsBuffer,0,t.norms);const h=C=>{const M=n.createBuffer({size:Math.max(16,C.byteLength),usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST});return n.queue.writeBuffer(M,0,C),M},d=[h(t.shells),h(t.prims),h(t.coefs),h(t.weights.length?t.weights:new Float32Array(4))],g=n.createBindGroup({layout:this.volumePipeline.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.paramsBuffer}},{binding:1,resource:{buffer:d[0]}},{binding:2,resource:{buffer:d[1]}},{binding:3,resource:{buffer:this.normsBuffer}},{binding:4,resource:{buffer:d[2]}},{binding:5,resource:{buffer:d[3]}},{binding:6,resource:this.volume.createView()}]}),b=n.createCommandEncoder(),p=b.beginComputePass();p.setPipeline(this.volumePipeline),p.setBindGroup(0,g),p.dispatchWorkgroups(Math.ceil(r/4),Math.ceil(o/4),Math.ceil(a/4)),p.end(),n.queue.submit([b.finish()]);for(const C of d)C.destroy()}render(t,n,r){const{device:o}=this,a=r??this.context?.getCurrentTexture();if(!a)return;this.resize(a.width,a.height),t.aspect=a.width/a.height;const i=new Float32Array(32);i.set(t.view,0),i.set(t.projection,16),o.queue.writeBuffer(this.cameraBuffer,0,i);const l=o.createCommandEncoder(),u=l.beginRenderPass({colorAttachments:[{view:a.createView(),clearValue:this.background,loadOp:"clear",storeOp:"store"}],depthStencilAttachment:{view:this.depth.createView(),depthClearValue:1,depthLoadOp:"clear",depthStoreOp:"store"}});if(this.bondCount>0&&this.bondBuffer&&(u.setPipeline(this.bondPipeline),u.setBindGroup(0,this.atomsBindGroup),u.setVertexBuffer(0,this.bondBuffer),u.setVertexBuffer(1,this.ringBuffer),u.draw(this.ringCount,this.bondCount)),this.sphereCount>0&&this.sphereBuffer&&(u.setPipeline(this.spherePipeline),u.setBindGroup(0,this.atomsBindGroup),u.setVertexBuffer(0,this.sphereBuffer),u.draw(4,this.sphereCount)),this.ghostBondCount>0&&this.ghostBondBuffer&&(u.setPipeline(this.ghostBondPipeline),u.setBindGroup(0,this.atomsBindGroup),u.setVertexBuffer(0,this.ghostBondBuffer),u.setVertexBuffer(1,this.ringBuffer),u.draw(this.ringCount,this.ghostBondCount)),this.ghostSphereCount>0&&this.ghostSphereBuffer&&(u.setPipeline(this.ghostSpherePipeline),u.setBindGroup(0,this.atomsBindGroup),u.setVertexBuffer(0,this.ghostSphereBuffer),u.draw(4,this.ghostSphereCount)),u.end(),this.volume&&this.box&&(n.showCloud||n.showOrbital)){const h=this.box,d=new Float32Array(36);d.set(je(t.viewProjection),0),d.set(t.eye,16),d[19]=n.cloudOpacity,d.set(h.min,20),d[23]=n.logMin,d.set(h.max,24),d[27]=n.logMax,d.set(h.max.map((p,C)=>(p-h.min[C])/Math.max(1,h.dims[C]-1)),28),d[31]=n.iso,d[32]=n.showCloud?1:0,d[33]=n.showOrbital?1:0,d[34]=n.orbitalOpacity,d[35]=this.frame++%64,o.queue.writeBuffer(this.rayBuffer,0,d);const g=o.createBindGroup({layout:this.rayPipeline.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.rayBuffer}},{binding:1,resource:this.volume.createView()},{binding:2,resource:this.sampler},{binding:3,resource:this.depth.createView()}]}),b=l.beginRenderPass({colorAttachments:[{view:a.createView(),loadOp:"load",storeOp:"store"}]});b.setPipeline(this.rayPipeline),b.setBindGroup(0,g),b.draw(3),b.end()}o.queue.submit([l.finish()])}async snapshot(t,n,r,o){const a=this.device.createTexture({size:[r,o],format:this.format,usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.COPY_SRC});this.render(t,n,a);const i=Math.ceil(r*4/256)*256,l=this.device.createBuffer({size:i*o,usage:GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ}),u=this.device.createCommandEncoder();u.copyTextureToBuffer({texture:a},{buffer:l,bytesPerRow:i},[r,o]),this.device.queue.submit([u.finish()]),await l.mapAsync(GPUMapMode.READ);const h=new Uint8Array(l.getMappedRange()),d=new Uint8ClampedArray(r*o*4),g=this.format.startsWith("bgra");for(let b=0;b<o;b++)for(let p=0;p<r;p++){const C=b*i+4*p,M=4*(b*r+p);d[M]=h[C+(g?2:0)],d[M+1]=h[C+1],d[M+2]=h[C+(g?0:2)],d[M+3]=255}return l.unmap(),l.destroy(),a.destroy(),d}async readVolume(){if(!this.volume||!this.box)return null;const[t,n,r]=this.box.dims,o=Math.ceil(t*8/256)*256,a=this.device.createBuffer({size:o*n*r,usage:GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ}),i=this.device.createCommandEncoder();i.copyTextureToBuffer({texture:this.volume},{buffer:a,bytesPerRow:o,rowsPerImage:n},[t,n,r]),this.device.queue.submit([i.finish()]),await a.mapAsync(GPUMapMode.READ);const l=new Uint16Array(a.getMappedRange()),u=new Float32Array(t*n*r*2);for(let h=0;h<r;h++)for(let d=0;d<n;d++)for(let g=0;g<t;g++){const b=(h*n+d)*o/2+g*4,p=2*((h*n+d)*t+g);u[p]=Ne(l[b]),u[p+1]=Ne(l[b+1])}return a.unmap(),a.destroy(),{dims:[t,n,r],data:u}}}function Ne(e){const t=e&32768?-1:1,n=e>>10&31,r=e&1023;return n===0?t*r*2**-24:n===31?r?NaN:t*(1/0):t*(1+r/1024)*2**(n-15)}const Ue=128,Gt=5;function It(e,t=128,n=Gt){const r=[1/0,1/0,1/0],o=[-1/0,-1/0,-1/0];for(let u=0;u<e.length;u+=3)for(let h=0;h<3;h++)r[h]=Math.min(r[h],e[u+h]-n),o[h]=Math.max(o[h],e[u+h]+n);const a=o.map((u,h)=>u-r[h]),i=Math.max(...a)/(t-1),l=a.map(u=>Math.min(t,Math.ceil(u/i-1e-9)+1));for(let u=0;u<3;u++){const h=(r[u]+o[u])/2,d=(l[u]-1)*i/2;r[u]=h-d,o[u]=h+d}return{min:r,max:o,dims:l}}function $t(e,t,n=e,r=128){const{nbf:o,nmo:a}=e;if(o>Ue)throw new Error(`The renderer supports at most ${Ue} basis functions (this system has ${o}).`);const{shells:i,prims:l,componentNorms:u}=e.basis,h=i.length/8,d=new ArrayBuffer(32*h),g=new Float32Array(d),b=new Uint32Array(d);for(let w=0;w<h;w++){const v=8*w,S=i[v+4],O=i[v+5];let f=1/0;for(let y=0;y<S;y++)f=Math.min(f,l[2*(O+y)]);g[v]=i[v],g[v+1]=i[v+1],g[v+2]=i[v+2],b[v+3]=i[v+3],b[v+4]=S,b[v+5]=O,b[v+6]=i[v+6],g[v+7]=f}const p=e.unrestricted?n.cb:n.ca,C=[];for(let w=0;w<e.nalpha;w++)C.push({c:n.ca,k:w,w:e.unrestricted?1:2});if(e.unrestricted)for(let w=0;w<e.nbeta;w++)C.push({c:p,k:w,w:1});const M=C.length,B=new Float32Array((M+1)*o);if(C.forEach(({c:w,k:v},S)=>{for(let O=0;O<o;O++)B[S*o+O]=w[O*a+v]}),t){const w=t.spin==="beta"?p:n.ca;for(let v=0;v<o;v++)B[M*o+v]=w[v*a+t.index]}return{nbf:o,nshell:h,nocc:M,shells:d,prims:Float32Array.from(l),norms:Float32Array.from(u.flat()),coefs:B,weights:Float32Array.from(C.map(w=>w.w)),box:It(e.coordinates,r)}}const ae=e=>`${(e*Z).toFixed(2)} eV`;function _t(e){const{homo:t,lumo:n}=Ot(e),r=Math.hypot(...e.dipole)*it,o=[{label:"Total energy",value:`${e.energy.toFixed(6)} Ha`,title:`${(e.energy*Z).toFixed(3)} eV. The energy of the electrons plus the repulsion between the nuclei; lower is more stable (1 Ha = 27.211 eV).`},{label:"HOMO",value:ae(t),title:"Energy of the highest occupied orbital. Its negative is roughly the energy needed to remove an electron."},{label:"LUMO",value:n===null?"none":ae(n),title:"Energy of the lowest empty orbital, where an added electron would go."},{label:"HOMO–LUMO gap",value:n===null?"n/a":ae(n-t),title:"Roughly how much energy it takes to excite the molecule."},{label:"Dipole moment",value:`${r.toFixed(3)} D`,title:"How unevenly the charge is spread, in debye. Water is about 1.85 D in experiment."},{label:"SCF",value:`${e.converged?"converged":"NOT converged"} in ${e.iterations} iterations, ${Math.round(e.millis)} ms`,title:"The self-consistent field loop: guess the orbitals, compute the field they make, solve again, until nothing changes."},{label:"Basis functions",value:String(e.nbf),title:"How many atom-centred functions the orbitals are built from. Bigger basis sets are more accurate and slower."}];return e.maxForce!==void 0&&o.push({label:"Max force",value:`${e.maxForce.toExponential(2)} Ha/bohr`}),o}const Vt={H:"hydrogen",He:"helium",Li:"lithium",Be:"beryllium",B:"boron",C:"carbon",N:"nitrogen",O:"oxygen",F:"fluorine",Ne:"neon",Na:"sodium",Mg:"magnesium",Al:"aluminium",Si:"silicon",P:"phosphorus",S:"sulfur",Cl:"chlorine",Ar:"argon"},Dt=([e,t,n])=>`rgb(${Math.round(e*255)}, ${Math.round(t*255)}, ${Math.round(n*255)})`;function qt(e){const t=[...new Set(e.map(n=>H(n.symbol).symbol))];return t.sort((n,r)=>H(r).z-H(n).z),t.map(n=>({symbol:n,name:Vt[n]??n,color:Dt(H(n).color)}))}function jt(e){return e==="HOMO"?"the highest-energy orbital that holds electrons":e==="LUMO"?"the lowest-energy empty orbital":e.startsWith("HOMO")?"an occupied orbital below the HOMO":"an empty orbital above the LUMO"}const Wt=1.6;function Yt(e,t=8,n=6){const o=(e.unrestricted?["alpha","beta"]:["alpha"]).flatMap(h=>{const d=h==="alpha"?e.nalpha:e.nbeta,g=De(e,h),b=g[d-1]?.energy??g[0].energy;return g.filter(p=>p.index>=d-t&&p.index<=d+n-1&&p.energy>b-Wt).map(p=>({...p,spin:h}))}),a=o.map(h=>h.energy),i=Math.min(...a),l=Math.max(...a),u=l-i||1;return o.map(h=>({...h,y:.04+.92*((l-h.energy)/u)}))}const U="http://www.w3.org/2000/svg";function Xt(e,t,n,r){e.replaceChildren();const o=240,a=300;e.setAttribute("viewBox",`0 0 ${o} ${a}`);const i=Yt(t),l=t.unrestricted?2:1,u=70,h=g=>l===1?70:g==="alpha"?40:125,d=new Map;for(const g of i){const b=`${g.spin}:${Math.round(g.energy*Z*50)}`;d.set(b,[...d.get(b)??[],g])}for(const g of d.values()){const b=(u-6*(g.length-1))/g.length;g.forEach((p,C)=>{const M=h(p.spin)+C*(b+6),B=p.y*(a-20)+10,w=n?.spin===p.spin&&n.index===p.index,v=document.createElementNS(U,"g");v.classList.add("level"),p.occupied&&v.classList.add("occupied"),w&&v.classList.add("selected");const S=document.createElementNS(U,"rect");S.setAttribute("x",String(M-2)),S.setAttribute("y",String(B-7)),S.setAttribute("width",String(b+4)),S.setAttribute("height","14"),S.setAttribute("class","hit");const O=document.createElementNS(U,"line");O.setAttribute("x1",String(M)),O.setAttribute("x2",String(M+b)),O.setAttribute("y1",String(B)),O.setAttribute("y2",String(B));const f=document.createElementNS(U,"title");if(f.textContent=`${p.label}${l===2?` (${p.spin==="alpha"?"α":"β"})`:""}: ${(p.energy*Z).toFixed(2)} eV`,v.append(S,O,f),p.occupied){const y=t.unrestricted?1:2;for(let x=0;x<y;x++){const k=document.createElementNS(U,"text");k.setAttribute("x",String(M+b/2+(y===2?x===0?-5:5:0))),k.setAttribute("y",String(B+4)),k.setAttribute("class","electron"),k.textContent=p.spin==="beta"!=(x===1)?"↓":"↑",v.append(k)}}if(C===g.length-1&&(p.label==="HOMO"||p.label==="LUMO"||w)){const y=document.createElementNS(U,"text");y.setAttribute("x",String(M+b+6)),y.setAttribute("y",String(B+4)),y.setAttribute("class","label"),y.textContent=p.label,v.append(y)}v.addEventListener("click",()=>r({spin:p.spin,index:p.index})),e.append(v)})}if(l===2)for(const[g,b]of[["alpha","α"],["beta","β"]]){const p=document.createElementNS(U,"text");p.setAttribute("x",String(h(g)+u/2)),p.setAttribute("y",String(a-2)),p.setAttribute("class","spin-label"),p.textContent=b,e.append(p)}}const We=.17;function ye(e){return .35+.3*Math.min(H(e).covalentRadius,1.15)*$}function D(e){return H(e).color}function Kt(e,t,n){const r=e.map((i,l)=>({center:i.xyz,radius:ye(i.symbol),color:D(i.symbol),highlight:l===t})),o=([i,l])=>n===void 0||Math.hypot(...[0,1,2].map(u=>e[i].xyz[u]-e[l].xyz[u]))<n*$,a=dt(e).filter(o).map(([i,l])=>({a:e[i].xyz,b:e[l].xyz,radius:We,colorA:D(e[i].symbol),colorB:D(e[l].symbol)}));return{spheres:r,bonds:a}}function Fe(e){const t=[0,1,2].map(r=>e.reduce((o,a)=>o+a.xyz[r],0)/Math.max(1,e.length)),n=e.reduce((r,o)=>Math.max(r,Math.hypot(o.xyz[0]-t[0],o.xyz[1]-t[1],o.xyz[2]-t[2])),0);return{center:t,radius:n}}const Ge=[.79,.65,.35];function Zt(e,t){const n=e.ghosts.map(o=>({center:o.xyz,radius:ye(o.symbol)*.6,color:D(o.symbol),highlight:!1,ghost:!0})),r=e.ghostBonds.map(({ghost:o,other:a,otherIsCluster:i})=>{const l=e.ghosts[o],u=i?t[a]:e.ghosts[a];return{a:l.xyz,b:u.xyz,radius:We*.45,colorA:D(l.symbol),colorB:D(u.symbol),ghost:!0}});for(const[o,a]of e.cellEdges)r.push({a:o,b:a,radius:.045,colorA:Ge,colorB:Ge});return{spheres:n,bonds:r}}const Jt=[{until:"drag",html:"<b>Try it:</b> drag one of the atoms away from the others. Every frame is a fresh quantum calculation."},{until:"relax",html:"The energy on the right went up: the molecule is strained. Press <b>Relax geometry</b> and watch the atoms settle back."},{until:"orbital",html:"Now look inside. Click a line under <b>Energy levels</b> on the right to see where that pair of electrons lives."},{until:"preset",html:"That shape is an orbital: blue and orange are the two signs of its wave. Next, pick another molecule from the list on the left."}];function Qt(e){return Jt.find(t=>!e.has(t.until))??null}const m=e=>document.getElementById(e),s={canvas:m("view"),noGpu:m("no-gpu"),preset:m("preset"),xyzBox:m("xyz-box"),xyz:m("xyz"),loadXyz:m("load-xyz"),charge:m("charge"),multiplicity:m("multiplicity"),method:m("method"),basis:m("basis"),relax:m("relax"),resetGeometry:m("reset-geometry"),showCloud:m("show-cloud"),cloudOpacity:m("cloud-opacity"),cloudMin:m("cloud-min"),showOrbital:m("show-orbital"),iso:m("iso"),status:m("status"),readouts:m("readouts"),orbitalName:m("orbital-name"),diagram:document.getElementById("diagram"),saveImage:m("save-image"),playScf:m("play-scf"),scfFrame:m("scf-frame"),scfLabel:m("scf-label"),help:m("help"),intro:m("intro"),legendCloud:m("legend-cloud"),legendOrbital:m("legend-orbital"),legendOrbitalName:m("legend-orbital-name"),legendAtoms:m("legend-atoms"),legendHelp:m("legend-help"),legendTitle:m("legend-title"),presetBlurb:m("preset-blurb"),settingsSummary:m("settings-summary"),energy:m("energy"),labels:m("labels"),tip:m("tip"),tipText:m("tip-text"),tipClose:m("tip-close"),restartTips:m("restart-tips"),introMore:m("intro-more"),settings:m("settings"),crystalBox:m("crystal-box"),crystalNote:m("crystal-note"),showCrystal:m("show-crystal"),crystalRepeats:m("crystal-repeats"),crystalRepeatsLabel:m("crystal-repeats-label"),legendCrystal:m("legend-crystal")},c={atoms:[],initialAtoms:[],charge:0,multiplicity:1,basis:"6-31g",method:"hf",result:null,orbital:null,frame:null,relaxId:null,hover:null,crystal:null},A=new Rt,R=new St(new Worker(new URL(""+new URL("engine.worker-CB3EoyU4.js",import.meta.url).href,import.meta.url),{type:"module"}));let L=null;const Ie=new URLSearchParams(location.search).has("headless"),Ye=Math.max(16,Math.min(160,Number(new URLSearchParams(location.search).get("res"))||112)),ue=Et(new URLSearchParams(location.search));R.configure(ue);let T=!0,j=null;function Xe(){return{showCloud:s.showCloud.checked,showOrbital:s.showOrbital.checked&&c.orbital!==null,cloudOpacity:40*Number(s.cloudOpacity.value)**2,logMin:Number(s.cloudMin.value),logMax:1,iso:Number(s.iso.value),orbitalOpacity:.55}}function P(e,t="ok"){s.status.textContent=e,s.status.className=`status ${t==="ok"?"":t}`}function be(e={}){return{atoms:c.atoms.map(t=>({symbol:t.symbol,xyz:t.xyz})),charge:c.charge,multiplicity:c.multiplicity,basis:c.basis,method:c.method,...e}}function I(){if(!L)return;const e=c.crystal,{spheres:t,bonds:n}=Kt(c.atoms,c.hover,e?e.spec.bond*1.15:void 0);if(e&&s.showCrystal.checked){const r=Zt(e.scene,c.atoms);t.push(...r.spheres),n.push(...r.bonds)}L.setAtoms(t,n),T=!0}function xe(){const e=Fe(c.atoms),t=c.crystal;if(!t||!s.showCrystal.checked)return A.fit(e.center,e.radius);const n=Fe([...c.atoms,...t.scene.ghosts]);A.fit(e.center,Math.max(e.radius*1.15,Math.min(n.radius,e.radius*1.9)))}function Ke(e){const t=e?_e[e]:void 0;if(s.crystalBox.hidden=!t,s.legendCrystal.hidden=!t||!s.showCrystal.checked,!t){c.crystal=null;return}s.crystalRepeats.max=String(t.maxRepeats),s.crystalRepeats.value=String(t.repeats),c.crystal={spec:t,scene:de(t,t.repeats)},ve()}function ve(){const e=c.crystal;if(!e)return;const t=Number(s.crystalRepeats.value);e.scene=de(e.spec,t),s.crystalRepeatsLabel.textContent=[0,1,2].map(a=>e.spec.fixedAxes?.includes(a)?1:t).join(" × ");const{cluster:n,caps:r}=e.scene,o=fe({atoms:n,charge:0});s.crystalNote.textContent=`Calculated: the ${n.length} bright atoms (${o} electrons)`+(r?`, ${n.length-r} from the crystal plus ${r} hydrogens on the cut bonds. `:". ")+"The faded atoms are the rest of the crystal, drawn from its structure but not calculated: the program solves finite molecules, not endless crystals.",s.legendCrystal.hidden=!s.showCrystal.checked}function we(){const e=c.result;if(!L||!e)return;const t=c.frame===null?null:e.history[c.frame];try{L.updateVolume($t(e,c.orbital,t??e,Ye))}catch(n){P(n instanceof Error?n.message:String(n),"error")}T=!0}function Ze(){L&&T&&(T=!1,tn(),L.render(A,Xe()),en()),requestAnimationFrame(Ze)}function en(){const e=c.atoms;for(;s.labels.children.length<e.length;)s.labels.append(document.createElement("span"));for(;s.labels.children.length>e.length;)s.labels.lastElementChild.remove();const t=A.viewProjection,n=s.canvas.getBoundingClientRect(),r=e.length>14;e.forEach((o,a)=>{const i=s.labels.children[a],[l,u,h]=o.xyz,d=t[3]*l+t[7]*u+t[11]*h+t[15];if(i.textContent=o.symbol,i.hidden=d<=0||r&&o.symbol==="H",i.hidden)return;const g=(t[0]*l+t[4]*u+t[8]*h+t[12])/d,b=(t[1]*l+t[5]*u+t[9]*h+t[13])/d;i.style.transform=`translate(${n.left+(g+1)/2*n.width}px, ${n.top+(1-b)/2*n.height}px) translate(-50%, -50%)`})}function tn(){const e=Math.min(window.devicePixelRatio||1,2),t=Math.max(1,Math.floor(s.canvas.clientWidth*e)),n=Math.max(1,Math.floor(s.canvas.clientHeight*e));(s.canvas.width!==t||s.canvas.height!==n)&&(s.canvas.width=t,s.canvas.height=n)}function Je(e){return{spin:"alpha",index:e.nalpha-1}}function Ce(e,t=!1){const n=c.result&&c.result.nmo===e.nmo&&c.result.unrestricted===e.unrestricted;c.result=e,(!n||!c.orbital||c.orbital.index>=e.nmo)&&(c.orbital=Je(e)),t||W(),c.frame=null,s.scfFrame.max=String(Math.max(0,e.history.length-1)),s.scfFrame.value=s.scfFrame.max,s.scfFrame.disabled=e.history.length===0,s.scfLabel.textContent=e.history.length?`step ${e.history.length} of ${e.history.length}`:"calculation steps",s.energy.textContent=`${e.energy.toFixed(5)} Ha`,s.readouts.replaceChildren(..._t(e).flatMap(r=>{const o=document.createElement("dt");o.textContent=r.label;const a=document.createElement("dd");return a.textContent=r.value,r.title&&(a.title=r.title),[o,a]})),Me(),we()}function Me(){const e=c.result;if(!e)return;const t=s.showOrbital.checked&&c.orbital!==null;if(Xt(s.diagram,e,t?c.orbital:null,n=>{c.orbital=n,s.showOrbital.checked=!0,Me(),we(),X("orbital")}),s.orbitalName.textContent="",t&&c.orbital){const n=De(e,c.orbital.spin)[c.orbital.index],r=`${n.label}${e.unrestricted?c.orbital.spin==="alpha"?" α":" β":""}`;s.orbitalName.textContent=`· ${r}`,s.legendOrbitalName.textContent=`${r} orbital`,s.legendOrbital.title=`${r}: ${jt(n.label)}. Blue where its wave is positive, orange where it is negative; where they meet, these electrons are never found.`}Q()}function Q(){s.legendTitle.textContent=s.preset.value==="custom"?"Your molecule":pe.find(e=>e.id===s.preset.value)?.name??"",s.legendCloud.hidden=!s.showCloud.checked,s.legendOrbital.hidden=!(s.showOrbital.checked&&c.orbital),s.legendAtoms.replaceChildren(...qt(c.atoms).map(e=>{const t=document.createElement("span");t.className="legend-item",t.title=e.name;const n=document.createElement("span");return n.className="swatch",n.style.background=e.color,t.append(n,e.symbol),t}))}function Qe(e){const t=c.atoms.length>=zt&&ue>1?`, integrals shared over ${ue} CPU threads`:"";return`Solved in ${Math.round(e.millis)} ms (${e.iterations} steps${t})`}async function _(e={}){if(c.atoms.length!==0){P("Solving…","busy");try{const t=await R.computeLatest(be({history:e.history??!0}));if(!t)return;Ce(t),P(t.converged?Qe(t):`The calculation did not settle after ${t.iterations} steps`,t.converged?"ok":"error")}catch(t){P(t instanceof Error?t.message:String(t),"error")}}}function et(e,t,n){Y(),c.atoms=e,c.initialAtoms=e.map(r=>({symbol:r.symbol,xyz:[...r.xyz]})),c.charge=t,c.multiplicity=n,c.result=null,c.orbital=null,s.charge.value=String(t),s.multiplicity.value=String(n),s.xyz.value=ft(e).split(`
`).slice(2).join(`
`),Be(),R.reset(),xe(),oe(),s.presetBlurb.textContent=Bt[s.preset.value]??"Your own molecule. Drag its atoms, relax it, and look at its orbitals like any other.",I(),Q(),_()}function tt(e){const t=Mt(e);Ke(t.crystal);for(const n of Array.from(s.basis.options))n.disabled=!t.bases.includes(n.value);t.bases.includes(c.basis)||(c.basis=t.bases[t.bases.length-1],s.basis.value=c.basis),et(Ve(t.xyz),t.charge,t.multiplicity)}function Be(){const e=c.multiplicity!==1;for(const t of Array.from(s.method.options))t.disabled=e&&t.value!=="hf";e&&c.method!=="hf"&&(c.method="hf",s.method.value="hf")}function Y(){c.relaxId!==null&&R.cancelRelax(c.relaxId),c.relaxId=null,s.relax.textContent="Relax geometry"}function nn(){if(c.relaxId!==null)return Y();s.relax.textContent="Stop",P("Relaxing…","busy"),c.relaxId=R.relax(be({history:!1}),{onStep:(e,t)=>{c.atoms=c.atoms.map((n,r)=>({symbol:n.symbol,xyz:e.coordinates.slice(3*r,3*r+3)})),I(),Ce(e),P(`Relaxing: step ${t+1}, E = ${e.energy.toFixed(6)} Ha, max force ${e.maxForce?.toExponential(1)}`,"busy")},onDone:(e,t)=>{c.relaxId=null,X("relax"),s.relax.textContent="Relax geometry",P(e?`Relaxed in ${t} steps · ${c.result.energy.toFixed(6)} Ha`:`Stopped after ${t} steps without converging`,e?"ok":"error")},onError:e=>{c.relaxId=null,s.relax.textContent="Relax geometry",P(e,"error")}})}function W(){j!==null&&clearInterval(j),j=null,s.playScf.textContent="▶ Watch it solve"}function he(e){const t=c.result;if(!t||t.history.length===0)return;c.frame=e>=t.history.length-1?null:e,s.scfFrame.value=String(e);const n=t.history[e];s.scfLabel.textContent=`step ${e+1} of ${t.history.length} · E = ${n.energy.toFixed(6)} Ha · error ${n.error.toExponential(1)}`,we()}async function on(){if(j!==null)return W();Y(),P("Recording the calculation…","busy");const e=await R.compute(be({history:!0,warmStart:!1})).catch(n=>(P(n.message,"error"),null));if(!e)return;Ce(e,!0),P(Qe(e));let t=0;s.playScf.textContent="■ Stop",he(0),j=window.setInterval(()=>{if(t++,t>=e.history.length)return W();he(t)},450)}function nt(e){const t=s.canvas.getBoundingClientRect();return[(e.clientX-t.left)/t.width*2-1,1-(e.clientY-t.top)/t.height*2]}function ot(e){const{origin:t,dir:n}=A.ray(...nt(e));let r=null,o=1/0;return c.atoms.forEach((a,i)=>{const l=Tt(t,n,a.xyz,ye(a.symbol)*1.1);l!==null&&l<o&&(o=l,r=i)}),r}let z=null;s.canvas.addEventListener("pointerdown",e=>{s.canvas.setPointerCapture(e.pointerId);const t=ot(e);t!==null?(Y(),W(),z={kind:"atom",index:t,planePoint:[...c.atoms[t].xyz],normal:A.basis().forward}):z={kind:"orbit",x:e.clientX,y:e.clientY},s.canvas.classList.add("dragging")});s.canvas.addEventListener("pointermove",e=>{if(!z){const l=ot(e);s.canvas.classList.toggle("over-atom",l!==null),l!==c.hover&&(c.hover=l,I());return}if(z.kind==="orbit"){A.orbit(e.clientX-z.x,e.clientY-z.y),z.x=e.clientX,z.y=e.clientY,T=!0;return}const{origin:t,dir:n}=A.ray(...nt(e)),r=G(n,z.normal);if(Math.abs(r)<1e-6)return;const o=G(J(z.planePoint,t),z.normal)/r,a=qe(t,me(n,o)),i=z.index;c.atoms=c.atoms.map((l,u)=>u===i?{symbol:l.symbol,xyz:a}:l),I(),_({history:!1})});const rt=()=>{z?.kind==="atom"&&X("drag"),z=null,s.canvas.classList.remove("dragging")};s.canvas.addEventListener("pointerup",rt);s.canvas.addEventListener("pointercancel",rt);s.canvas.addEventListener("wheel",e=>{e.preventDefault(),A.zoom(Math.exp(e.deltaY*.001)),T=!0},{passive:!1});window.addEventListener("resize",()=>T=!0);for(const e of Object.keys(Le)){const t=document.createElement("optgroup");t.label=Le[e];for(const n of pe.filter(r=>r.group===e))t.append(new Option(n.name,n.id));s.preset.append(t)}s.preset.add(new Option("Custom (pasted XYZ)","custom"));s.preset.addEventListener("change",()=>{s.preset.value==="custom"?(s.settings.open=!0,s.xyzBox.open=!0):(tt(s.preset.value),X("preset"))});s.loadXyz.addEventListener("click",()=>{try{const e=Ve(s.xyz.value);if(e.length===0)throw new Error('No atoms found. Use lines like "O 0 0 0" in Angstrom.');s.preset.value="custom",Ke(void 0);for(const n of Array.from(s.basis.options))n.disabled=!1;const t=fe({atoms:e,charge:c.charge,multiplicity:1});et(e,c.charge,t%2===0?1:2)}catch(e){P(e instanceof Error?e.message:String(e),"error")}});s.charge.addEventListener("change",()=>{c.charge=Math.round(Number(s.charge.value)),(fe({atoms:c.atoms,charge:c.charge})+c.multiplicity-1)%2!==0&&(c.multiplicity=c.multiplicity===1?2:c.multiplicity-1),s.multiplicity.value=String(c.multiplicity),Be(),_()});s.multiplicity.addEventListener("change",()=>{c.multiplicity=Math.max(1,Math.round(Number(s.multiplicity.value))),Be(),_()});s.method.addEventListener("change",()=>{c.method=s.method.value,oe(),_()});s.basis.addEventListener("change",()=>{c.basis=s.basis.value,oe(),_()});s.relax.addEventListener("click",nn);s.showCrystal.addEventListener("input",()=>{ve(),I(),Q(),xe()});s.crystalRepeats.addEventListener("input",()=>{ve(),I(),xe()});s.resetGeometry.addEventListener("click",()=>{Y(),c.atoms=c.initialAtoms.map(e=>({symbol:e.symbol,xyz:[...e.xyz]})),I(),_()});for(const e of[s.cloudOpacity,s.cloudMin,s.iso])e.addEventListener("input",()=>T=!0);s.showCloud.addEventListener("input",()=>{Q(),T=!0});s.showOrbital.addEventListener("input",()=>{s.showOrbital.checked&&!c.orbital&&c.result&&(c.orbital=Je(c.result)),s.showOrbital.checked&&X("orbital"),Me(),T=!0});s.playScf.addEventListener("click",on);s.saveImage.addEventListener("click",async()=>{const e=Math.min(window.devicePixelRatio||1,2),t=document.createElement("a");t.href=await at(Math.round(s.canvas.clientWidth*e),Math.round(s.canvas.clientHeight*e)),t.download=`${s.preset.value}-${c.method}-${c.basis}.png`,t.click()});s.scfFrame.addEventListener("input",()=>{W(),he(Number(s.scfFrame.value))});const st="quantum-molecules:intro-seen";function Oe(e){s.intro.open||s.intro.showModal(),s.tip.hidden=!0,e?(s.introMore.open=!0,document.getElementById(e)?.scrollIntoView({block:"start"})):s.intro.scrollTop=0}s.help.addEventListener("click",()=>Oe());s.legendHelp.addEventListener("click",()=>Oe("intro-colors"));s.intro.addEventListener("close",()=>{ee(st),ne()});s.intro.addEventListener("click",e=>{e.target===s.intro&&s.intro.close()});function ee(e,t=!0){try{t?localStorage.setItem(e,"1"):localStorage.removeItem(e)}catch{}}function $e(e){try{return localStorage.getItem(e)!==null}catch{return!1}}const te="quantum-molecules:tips-done",Se=new Set;let q=!1;function ne(){const e=q?Qt(Se):null;e===null&&q&&(q=!1,ee(te)),s.tip.hidden=e===null||s.intro.open,e&&(s.tipText.innerHTML=e.html)}function X(e){Se.add(e),ne()}s.tipClose.addEventListener("click",()=>{q=!1,ee(te),ne()});s.restartTips.addEventListener("click",()=>{Se.clear(),q=!0,ee(te,!1),s.intro.close()});function oe(){s.settingsSummary.textContent=`(${s.method.selectedOptions[0]?.text}, ${s.basis.selectedOptions[0]?.text})`}async function rn(){const e=!Ie||new URLSearchParams(location.search).has("intro");q=e&&!$e(te),e&&!$e(st)?Oe():ne(),s.method.value=c.method,s.basis.value=c.basis,oe();try{L=await ge.create(Ie?null:s.canvas),L.device.lost.then(t=>P(`The GPU device was lost: ${t.message}`,"error"))}catch(t){s.noGpu.hidden=!1,s.noGpu.textContent=`${t instanceof Error?t.message:String(t)} The calculations still run and the numbers on the right are real, but the 3D view needs a browser with WebGPU (recent Chrome, Edge, or Safari).`}s.preset.value="h2o",await R.ready,tt("h2o"),requestAnimationFrame(Ze)}async function at(e=s.canvas.clientWidth||800,t=s.canvas.clientHeight||600){if(!L)throw new Error("WebGPU is not available.");A.aspect=e/t;const n=await L.snapshot(A,Xe(),e,t),r=document.createElement("canvas");return r.width=e,r.height=t,r.getContext("2d").putImageData(new ImageData(n,e,t),0,0),r.toDataURL("image/png")}Object.assign(window,{__qm:{state:c,engine:R,camera:A,snapshot:at,volumeResolution:Ye,get renderer(){return L}}});rn();
