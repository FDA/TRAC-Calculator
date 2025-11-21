document.addEventListener('DOMContentLoaded', function() {
  try{ if(window.ChartAnnotation){ Chart.register(window.ChartAnnotation); } }catch(e){}

  let hillChart, modelChart;
  let Prob1=0, cmaxIsNM=false;

  // collapse toggles
  const cmaxHeader=document.getElementById('cmaxSectionHeader'),
        cmaxContent=document.getElementById('cmaxSectionContent'),
        hillHeader=document.getElementById('hillSectionHeader'),
        hillContent=document.getElementById('hillSectionContent');
  cmaxHeader.addEventListener('click',()=>{
    if(cmaxContent.style.display==='none'){cmaxContent.style.display='';cmaxHeader.textContent='▾ Cmax Interpolation';}
    else{cmaxContent.style.display='none';cmaxHeader.textContent='▸ Cmax Interpolation';}
  });
  hillHeader.addEventListener('click',()=>{
    if(hillContent.style.display==='none'){hillContent.style.display='';hillHeader.textContent='▾ Hill Fit Curve (only available when input is filled in Cmax Interpolation)';}
    else{hillContent.style.display='none';hillHeader.textContent='▸ Hill Fit Curve (only available when input is filled in Cmax Interpolation)';}
  });

  // QT collapse
  const qtHeader = document.getElementById('qtSectionHeader'),
        qtContent = document.getElementById('qtSectionContent');
  if (qtHeader && qtContent) {
    qtHeader.addEventListener('click', () => {
      if (qtContent.style.display === 'none') {
        qtContent.style.display = '';
        qtHeader.textContent = '▾ QT Prolongation Prediction (in progress)';
      } else {
        qtContent.style.display = 'none';
        qtHeader.textContent = '▸ QT Prolongation Prediction (in progress)';
      }
    });
  }


  // switch unit
  document.getElementById('switchCmaxUnit').addEventListener('click',()=>{const inp=document.getElementById('cmax'), lbl=document.getElementById('cmaxLabel');let v=parseFloat(inp.value); if(isNaN(v))return;if(cmaxIsNM){v/=1000;lbl.innerHTML='Cmax (<span class="plain-greek">µ</span>M)';const ch=document.getElementById('concHeader'); if(ch) ch.innerHTML='Concentration (<span class="plain-greek">µ</span>M)';}else{v*=1000;lbl.innerText='Cmax (nM)'; const ch=document.getElementById('concHeader'); if(ch) ch.innerText='Concentration (nM)';}inp.value=v.toFixed(4); cmaxIsNM=!cmaxIsNM;});

  window.addRow=()=>{const tb=document.getElementById('dataBody'),r=document.createElement('tr');r.innerHTML='<td><input name="concentration[]" type="number" step="any" required></td><td><input name="fpdc[]" type="number" step="any" required></td><td><button type="button" onclick="removeRow(this)">−</button></td>';tb.appendChild(r);};
  window.removeRow=btn=>btn.closest('tr').remove();
  document.getElementById('predictorCalcBtn').addEventListener('click',()=>calculate(true));
  document.getElementById('riskForm').addEventListener('submit',e=>{e.preventDefault();calculate(false);});

  function calculate(predOnly){
    let arr,p4,p7,cell=0,assay='30',Cmax,concs=[],fpdcs=[];
    if(predOnly){
      arr=parseInt(document.getElementById('predictor1').value);
      p4=parseFloat(document.getElementById('predictor4').value);
      p7=parseFloat(document.getElementById('predictor7').value);
      if(isNaN(arr)||isNaN(p4)||isNaN(p7)){return;}validatePredictorRanges();const c=document.getElementById('cmax'); if(c) c.value='';document.querySelectorAll('#dataBody input').forEach(el=>el.value='');
      if(p4===0&&p7===0){}
    } else {
      Cmax=parseFloat(document.getElementById('cmax').value);
      if(cmaxIsNM)Cmax/=1000;
      arr=parseInt(document.getElementById('arrhythmia').value);
      const cellEl=document.getElementById('celltype'); cell=cellEl?parseFloat(cellEl.value):0;
      
      document.querySelectorAll('#dataBody tr').forEach(r=>{const [ci,fi]=r.querySelectorAll('input'),c=parseFloat(ci.value),f=parseFloat(fi.value);if(!isNaN(c)&&!isNaN(f)){concs.push(c);fpdcs.push(f);}});
      // Check Cmax within input concentration range
      if (!isNaN(Cmax)) {
        const minConc = Math.min(...concs);
        const maxConc = Math.max(...concs);
        if (Cmax < minConc || Cmax > maxConc) {
          alert('STOP\n[Cmax Interpolation]\nPlease enter Concentration – Repolarization values with a range that covers the Cmax.');
          return;
        }
      }
      if(concs.length<4){alert('STOP\n[Cmax Interpolation]\nPlease enter at least 4 pairs of Concentration – Repolarization values for the Hill\'s Fit to converge.');return;}
      p4=Math.max(...fpdcs);
      // Determine trend and set Top/Bottom to allow negative or positive Hill fit
      const minY = Math.min(...fpdcs), maxY = Math.max(...fpdcs);
      // Simple slope estimate
      const meanX = concs.reduce((a,b)=>a+b,0)/concs.length;
      const meanY = fpdcs.reduce((a,b)=>a+b,0)/fpdcs.length;
      const slope = concs.map((x,i)=>(x-meanX)*(fpdcs[i]-meanY)).reduce((a,b)=>a+b,0) /
                    concs.map(x=>(x-meanX)*(x-meanX)).reduce((a,b)=>a+b,1e-9);
      const decreasing = slope < 0;
      const Bottom = decreasing ? maxY : minY;
      const Top    = decreasing ? minY : maxY;

      const med = a => { const s=[...a].sort((x,y)=>x-y), m=Math.floor(s.length/2); return a.length%2?s[m]:(s[m-1]+s[m])/2; };
      const guess = { Bottom, Top, EC50: med(concs), Hill: 1 }; // lock Hill coefficient at 1

      // Hill function with Hill fixed at 1
      const hillf = (x,p) => p.Bottom + (p.Top - p.Bottom) / (1 + (p.EC50 / x));

      // Fit EC50 on a grid while Hill fixed at 1
      let best = { ...guess }, minE = Infinity;
      const loss = p => fpdcs.reduce((s,y,i)=>s + Math.pow(hillf(concs[i],p) - y, 2), 0);
      for (let ec = Math.min(...concs)*0.1; ec <= Math.max(...concs)*10; ec += (Math.max(...concs)-Math.min(...concs))/100) {
        const t = { ...guess, EC50: ec };
        const e = loss(t);
        if (e < minE) { minE = e; best = t; }
      }

      const FPDc = hillf(Cmax||Math.min(...concs), best);p7=FPDc;document.getElementById('predictor1').value=String(arr);document.getElementById('predictor4').value=isFinite(p4)?Number(p4).toFixed(4):'';document.getElementById('predictor7').value=isFinite(p7)?Number(p7).toFixed(4):'';validatePredictorRanges();
      if(p4===0&&p7===0){}
      const Thr=assay==='30'?Bottom*1.103:Bottom*1.0794;
      const logM=assay==='30'?(Thr+0.35)/0.92:(Thr+0.17)/0.93;
      (()=>{const el=document.getElementById('estimatedQTc'); if(el){ el.innerHTML=`<strong>QTc (log M):</strong> ${logM.toFixed(4)}<br><strong>Conc >10ms QT:</strong> ${Math.pow(10,logM).toFixed(4)} µM`; }})();
      const fitX=Array.from({length:100},(_,i)=>Math.pow(10,Math.log10(Math.max(0.001,Math.min(...concs)))+i*(Math.log10(Math.max(...concs))-Math.log10(Math.max(0.001,Math.min(...concs))))/99));
      const fitY=fitX.map(x=>hillf(x,best));
      if(hillChart)hillChart.destroy();
      hillChart=new Chart(document.getElementById('hillPlot'),{type:'line',data:{labels:fitX,datasets:[{label:'Hill Fit',data:fitX.map((x,i)=>({x,y:fitY[i]})),borderWidth:3,fill:false},{label:'Data',type:'scatter',data:concs.map((x,i)=>({x,y:fpdcs[i]})),pointRadius:4},{label:'Cmax',type:'scatter',data:[{x:Cmax,y:FPDc}],pointRadius:6}]},options:{responsive:true,maintainAspectRatio:false,scales:{x:{type:'logarithmic', grid:{lineWidth:5}, ticks:{font:{size:20}}, title:{display:true, text:'Concentration (µM)', font:{size:18}}},y:{grid:{lineWidth:5}, ticks:{font:{size:20}}, title:{display:true, text:'ΔΔFPDc or ΔΔAPD90c (ms)', font:{size:18}}}}}});
    }
    // model probabilities (Model 1 only)
const map1=[0,0.6583,1.7944];
const logit1=-0.1311+map1[arr]+0.00687*p4+0.0232*p7;
Prob1=1/(1+Math.exp(-logit1));
if(Prob1<0){}
updateModelPanel();
}

  
  
function updateModelPanel(){
  const title=document.getElementById('modelTitle'),
        sub=document.getElementById('modelSubtitle'),
        res=document.getElementById('modelResults');
  const labels=['High or Intermediate','Low'];
  const data=[Prob1*100,(1-Prob1)*100];
  if(title) title.innerText='Model 1 TdP Risk';
  if(sub)   sub.innerHTML='The background calculation uses a logistic regression model. The model outputs are:';
  const colors = labels.map(l => (l.includes('High') ? 'rgb(230,75,53)' : 'rgb(3,160,135)'));
  if(res){ res.innerHTML = '<ul style="margin-left:20px;">' +
    labels.map((l,i)=>`<li><strong>${l} TdP Risk Probability:</strong> ${data[i].toFixed(1)}%</li>`).join('') +
    '</ul>'; }
  if(modelChart) modelChart.destroy();
  const datasets = labels.map((l,i)=>({ label: l + ' Risk', data: [data[i]], backgroundColor: colors[i], stack: 'risk', borderWidth: 0, borderRadius: 6, maxBarThickness: 72 }));
  modelChart = new Chart(document.getElementById('modelChart'), {
    type: 'bar',
    data: { labels: ['Predicted Risk'], datasets },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      layout: { padding: {top:8,right:8,bottom:4,left:8} },
      scales: {
        x: { stacked: true, grid: { display:false }, ticks: { display:false } },
        y:{ stacked: true, beginAtZero: true, max: 100, grid:{color:'rgba(0,0,0,0.08)', lineWidth:5}, ticks: { font: { size: 20 } }, title:{display:true, text:'Predicted Risk Probability (%)', font:{size:18}} }
      },
      plugins: {
        legend: { display: true, position: 'chartArea', align: 'end', labels:{ boxWidth:14, boxHeight:14, useBorderRadius:true, borderRadius:3, font:{size:14} } },
        tooltip: { callbacks: { label: ctx => `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(1)}%` } },
        annotation: { annotations: { riskThreshold: { type: 'line', yMin: 80, yMax: 80, borderColor: 'red', borderWidth: 2, borderDash:[6,6], label:{ display:true, content:'Risk Probability Threshold (80%)', position:'end', color:'red', font:{ size:14 }, yAdjust:-6 } } } }
      }
    }
  });
}


  const p4el = document.getElementById('predictor4');
  if (p4el) p4el.addEventListener('change', validatePredictorRanges);
  const p7el = document.getElementById('predictor7');
  if (p7el) p7el.addEventListener('change', validatePredictorRanges);
});

// V3.70: Range validation for Predictor4 and Predictor7
function validatePredictorRanges(){
  const p4 = parseFloat(document.getElementById('predictor4').value);
  if(!isNaN(p4) && (p4 < -372 || p4 > 1280)){
    alert('WARNING\n[Predictor Inputs]\nPredictor Inputs are outside the range for the prediction model to yield risk probability within acceptable confidence interval. Please enter Predictor Inputs within the designated ranges.');
  }
  const p7 = parseFloat(document.getElementById('predictor7').value);
  if(!isNaN(p7) && (p7 < -100 || p7 > 303)){
    alert('WARNING\n[Predictor Inputs]\nPredictor Inputs are outside the range for the prediction model to yield risk probability within acceptable confidence interval. Please enter Predictor Inputs within the designated ranges.');
  }
}