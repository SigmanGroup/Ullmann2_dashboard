const FP_OPTIONS=JSON.stringify({radius:2,nBits:2048});
let DATA=null,RDKit=null;

const $=id=>document.getElementById(id);
function status(s,bad=false){$("status").textContent=s;$("status").style.color=bad?"#b42318":"#667085"}

function bitsFromString(fp){
  const s=new Set();
  for(let i=0;i<fp.length;i++) if(fp[i]==="1") s.add(i);
  return s;
}
function tanimoto(a,b){
  let inter=0;
  for(const x of a) if(b.has(x)) inter++;
  return inter/(a.size+b.size-inter || 1);
}
function molInfo(smiles){
  const mol=RDKit.get_mol(smiles);
  if(!mol) throw new Error("Invalid SMILES");
  const canonical=mol.get_smiles();
  const fp=bitsFromString(mol.get_morgan_fp(FP_OPTIONS));
  const svg=mol.get_svg(300,180);
  mol.delete();
  return {canonical,fp,svg};
}
function matchLibrary(q,library,canonKey){
  let exact=library.find(r=>r.canonical_smiles===q.canonical);
  if(exact) return {row:exact,score:1,exact:true};
  let best=null,bestScore=-1;
  for(const r of library){
    const b=new Set(r.fp_onbits||[]);
    const score=tanimoto(q.fp,b);
    if(score>bestScore){bestScore=score;best=r}
  }
  return {row:best,score:bestScore,exact:false};
}
function info(target,m,kind){
  const r=m.row;
  const id=kind==="amine"?r.Numerical_ID:r.id;
  const name=kind==="amine" && r.Compound_Name ? `<div class="info-line"><b>Name:</b> ${r.Compound_Name}</div>`:"";
  $(target).innerHTML=`<div class="info-line"><b>${m.exact?"Exact library match":"Closest library match"}</b></div>
  ${name}<div class="info-line"><b>Library ID:</b> ${id}</div>
  <div class="info-line"><b>Tanimoto similarity:</b> ${m.score.toFixed(3)}</div>
  <div class="info-line"><b>Matched SMILES:</b> ${r.canonical_smiles}</div>`;
}
function project(am,br){
  const av=DATA.pca.am_columns.map(c=>Number(am[c.replace(/^Am_/,"")]));
  const bv=DATA.pca.br_columns.map(c=>Number(br[c.replace(/^Br_/,"")]));
  if([...av,...bv].some(v=>!Number.isFinite(v))) throw new Error("Matched library entry is missing a PCA descriptor.");
  const z=[
    ...av.map((v,i)=>(v-DATA.pca.amine_mean[i])/DATA.pca.amine_scale[i]),
    ...bv.map((v,i)=>(v-DATA.pca.bromide_mean[i])/DATA.pca.bromide_scale[i])
  ];
  return DATA.pca.components.map(comp=>comp.reduce((s,w,i)=>s+w*z[i],0));
}
function thresholdPlot(user=null){
  const groups=["HTE Round 1","HTE Round 2","Validation"];
  const traces=groups.map(g=>{
    const p=DATA.threshold_points.filter(d=>d.group===g);
    return {
      x:p.map(d=>d.x),
      y:p.map(d=>d.y),
      text:p.map(d=>d.product),
      customdata:p.map(d=>({
        product:d.product,
        amine_smiles:d.amine_smiles,
        bromide_smiles:d.bromide_smiles,
        yield:d.yield
      })),
      mode:"markers",
      type:"scatter",
      name:g,
      marker:{
        size:g==="Validation" ? p.map(d=>Math.max(8,2*Math.sqrt((Number(d.yield)||0)*7/Math.PI))) : 10,
        symbol:g==="Validation"?"x":"circle",
        color:p.map(d=>d.yield),
        colorscale:"GnBu",
        cmin:0,
        cmax:100,
        showscale:g==="HTE Round 1",
        colorbar:g==="HTE Round 1"?{title:"Average yield"}:undefined,
        line:{color:"#777",width:.6}
      },
      hovertemplate:
        "<b>%{text}</b><br>Amine descriptor: %{x:.4f}<br>Bromide descriptor: %{y:.2f}<br>Average yield: %{customdata.yield:.1f}<extra></extra>"
    };
  });

  if(user) traces.push({
    x:[user.x],
    y:[user.y],
    mode:"markers+text",
    text:["User"],
    textposition:"top right",
    name:"User input",
    marker:{size:18,color:"red",line:{color:"black",width:2}},
    hovertemplate:"User input<extra></extra>"
  });

  Plotly.react(
    "thresholdPlot",
    traces,
    {
      margin:{l:70,r:30,t:55,b:65},
      hovermode:"closest",
      xaxis:{title:"Amine_LP_energy_min"},
      yaxis:{title:"Bromide_buried_vol_3.0Å_boltz"},
      shapes:[
        {type:"line",x0:DATA.thresholds.x,x1:DATA.thresholds.x,y0:0,y1:1,yref:"paper",line:{color:"black",dash:"dash"}},
        {type:"line",y0:DATA.thresholds.y,y1:DATA.thresholds.y,x0:0,x1:1,xref:"paper",line:{color:"firebrick",dash:"dash"}}
      ],
      legend:{
        orientation:"h",
        x:0.5,
        xanchor:"center",
        y:1.11,
        yanchor:"bottom"
      }
    },
    {responsive:true,displaylogo:false}
  ).then(() => attachStructureHover("thresholdPlot"));
}

function renderSmilesSVG(smiles, width=210, height=130){
  if(!smiles || !RDKit) return "";
  try{
    const mol = RDKit.get_mol(smiles);
    if(!mol) return "";
    const svg = mol.get_svg(width, height);
    mol.delete();
    return svg;
  }catch(e){
    console.warn("Could not render structure", e);
    return "";
  }
}

function showStructureHover(pointData, evt){
  const card = document.getElementById("structure-hover");
  if(!card || !pointData) return;

  const amine = pointData.amine_smiles || "";
  const bromide = pointData.bromide_smiles || "";
  const product = pointData.product || "";

  if(!amine && !bromide){
    card.classList.add("hidden-hover");
    return;
  }

  document.getElementById("hover-product").textContent = product;
  document.getElementById("hover-amine-svg").innerHTML = renderSmilesSVG(amine);
  document.getElementById("hover-bromide-svg").innerHTML = renderSmilesSVG(bromide);
  document.getElementById("hover-amine-smiles").textContent = amine;
  document.getElementById("hover-bromide-smiles").textContent = bromide;

  card.classList.remove("hidden-hover");

  // Place near pointer but keep inside viewport.
  const pad = 14;
  const w = 470;
  const h = 270;
  let x = (evt && evt.clientX ? evt.clientX : 40) + 16;
  let y = (evt && evt.clientY ? evt.clientY : 40) + 16;

  if(x + w > window.innerWidth - pad) x = Math.max(pad, x - w - 32);
  if(y + h > window.innerHeight - pad) y = Math.max(pad, window.innerHeight - h - pad);

  card.style.left = `${x}px`;
  card.style.top = `${y}px`;
}

function hideStructureHover(){
  const card = document.getElementById("structure-hover");
  if(card) card.classList.add("hidden-hover");
}

function attachStructureHover(plotId){
  const plot = document.getElementById(plotId);
  if(!plot || plot.__structureHoverAttached) return;
  plot.__structureHoverAttached = true;

  plot.on("plotly_hover", ev => {
    const p = ev.points && ev.points[0];
    if(!p || !p.customdata) return hideStructureHover();

    // customdata is an object for these traces
    showStructureHover(p.customdata, ev.event);
  });

  plot.on("plotly_unhover", hideStructureHover);
}

function pcaPlot(user=null){
  const pts = DATA.pca_points || [];

  const tested = pts.filter(d => d.group === "tested");
  const validation = pts.filter(d => d.group === "validation");

  const sizeFromYield = y => {
    const v = Number(y);
    if (!Number.isFinite(v) || v <= 0) return 5;
    return Math.max(5, 2 * Math.sqrt((v * 7) / Math.PI));
  };

  const traces = [
    {
      x: tested.map(d => d.pc1),
      y: tested.map(d => d.pc2),
      text: tested.map(d => d.product),
      customdata: tested.map(d => ({
        product: d.product,
        amine_smiles: d.amine_smiles,
        bromide_smiles: d.bromide_smiles,
        yield: d.yield
      })),
      mode: "markers",
      type: "scatter",
      name: "Tested products",
      marker: {
        size: tested.map(d => sizeFromYield(d.yield)),
        sizemode: "diameter",
        color: "sandybrown",
        opacity: 0.8,
        line: { color: "grey", width: 0.5 }
      },
      hovertemplate:
        "<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>Average yield: %{customdata.yield:.1f}<extra></extra>"
    },
    {
      x: validation.map(d => d.pc1),
      y: validation.map(d => d.pc2),
      text: validation.map(d => d.product),
      customdata: validation.map(d => ({
        product: d.product,
        amine_smiles: d.amine_smiles,
        bromide_smiles: d.bromide_smiles,
        yield: d.yield
      })),
      mode: "markers+text",
      type: "scatter",
      textposition: "bottom right",
      name: "Validation products",
      marker: {
        size: validation.map(d => sizeFromYield(d.yield)),
        sizemode: "diameter",
        color: "steelblue",
        opacity: 0.65,
        line: { color: "black", width: 1 }
      },
      hovertemplate:
        "<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>Average yield: %{customdata.yield:.1f}<extra></extra>"
    }
  ];

  if (user) {
    traces.push({
      x: [user[0]],
      y: [user[1]],
      mode: "markers+text",
      type: "scatter",
      text: ["User input"],
      textposition: "top right",
      name: "User input",
      marker: {
        size: 18,
        color: "red",
        line: { color: "black", width: 2 }
      },
      hovertemplate:
        "User input<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>"
    });
  }

  // Keep limits based on all PCA points so removing gray x marks
  // does not change the manuscript-space framing.
  const allX = pts.map(d => Number(d.pc1)).filter(Number.isFinite);
  const allY = pts.map(d => Number(d.pc2)).filter(Number.isFinite);
  if (user) {
    allX.push(Number(user[0]));
    allY.push(Number(user[1]));
  }

  const xmin = Math.min(...allX), xmax = Math.max(...allX);
  const ymin = Math.min(...allY), ymax = Math.max(...allY);
  const span = Math.max(xmax - xmin, ymax - ymin);
  const xmid = (xmin + xmax) / 2;
  const ymid = (ymin + ymax) / 2;
  const half = span / 2 + 0.35;

  Plotly.react(
    "pcaPlot",
    traces,
    {
      width: 560,
      height: 560,
      margin: { l: 60, r: 20, t: 70, b: 60 },
      plot_bgcolor: "white",
      paper_bgcolor: "white",
      xaxis: {
        title: "PC1",
        range: [xmid - half, xmid + half],
        showline: true,
        mirror: true,
        linewidth: 1,
        linecolor: "black",
        zeroline: false,
        showgrid: false,
        ticks: "outside"
      },
      yaxis: {
        title: "PC2",
        range: [ymid - half, ymid + half],
        showline: true,
        mirror: true,
        linewidth: 1,
        linecolor: "black",
        zeroline: false,
        showgrid: false,
        ticks: "outside",
        scaleanchor: "x",
        scaleratio: 1
      },
      legend: {
        orientation: "h",
        x: 0.5,
        xanchor: "center",
        y: 1.13,
        yanchor: "bottom",
        font: { size: 12 },
        bgcolor: "rgba(255,255,255,0)"
      },
      hovermode: "closest"
    },
    {
      responsive: true,
      displaylogo: false
    }
  ).then(() => attachStructureHover("pcaPlot"));
}
async function analyze(){
  try{
    const as=$("amine").value.trim(),bs=$("bromide").value.trim();
    if(!as||!bs) throw new Error("Please enter both SMILES.");
    status("Analyzing…");
    const aq=molInfo(as),bq=molInfo(bs);
    const am=matchLibrary(aq,DATA.amines,"canonical_smiles");
    const br=matchLibrary(bq,DATA.bromides,"canonical_smiles");
    $("matches").classList.remove("hidden");
    $("amineSvg").innerHTML=aq.svg;$("bromideSvg").innerHTML=bq.svg;
    info("amineInfo",am,"amine");info("bromideInfo",br,"bromide");
    const x=Number(am.row.NBO_LP_energy_N1_min),y=Number(br.row["%Vbur_Br_3.0Å_Boltz"]);
    thresholdPlot({x,y});
    const pc=project(am.row,br.row);pcaPlot(pc);
    const inside=x<=DATA.thresholds.x && y<=DATA.thresholds.y;
    $("prediction").textContent=inside?"Within threshold-defined region":"Outside threshold-defined region";
    $("prediction").className="badge "+(inside?"good":"bad");
    status(`Complete. Amine similarity ${am.score.toFixed(3)}; bromide similarity ${br.score.toFixed(3)}.`);
  }catch(e){status(e.message,true);console.error(e)}
}
Promise.all([
  fetch("data/dashboard_data.json").then(r=>r.json()),
  initRDKitModule()
]).then(([data,rdkit])=>{
  DATA=data;RDKit=rdkit;thresholdPlot();pcaPlot();
  $("go").disabled=false;$("go").textContent="Analyze";status("Ready.");
}).catch(e=>status("Could not load dashboard: "+e.message,true));
$("go").addEventListener("click",analyze);
