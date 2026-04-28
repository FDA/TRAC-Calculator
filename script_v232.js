document.addEventListener('DOMContentLoaded', function() {
  try{ if(window.ChartAnnotation){ Chart.register(window.ChartAnnotation); } }catch(e){}

  let hillChart, modelChart;
  let Prob1=0, cmaxIsNM=false;
  
  // Store CI values globally
  window.ciLower = null;
  window.ciUpper = null;

  // Covariance matrix from R code (TRAC_output_CI.R)
  const covBeta = [
    [0.024975,  0.000023, -0.00019, -0.0238,  -0.02417],
    [0.000023,  2.636e-6, -4.31e-6, -0.00014, -0.00015],
    [-0.00019, -4.31e-6,  0.000055, -0.00034, -0.00022],
    [-0.0238,  -0.00014, -0.00034,  0.188749,  0.037828],
    [-0.02417, -0.00015, -0.00022,  0.037828,  0.176284]
  ];

  // Function to calculate confidence intervals
  function calculateCI(arr, p4, p7, alpha = 0.05) {
    // Create design vector based on arrhythmia type
    let X_new;
    if (arr === 1) {
      X_new = [1, p4, p7, 1, 0];  // Type A arrhythmia
    } else if (arr === 2) {
      X_new = [1, p4, p7, 0, 1];  // Other arrhythmia type
    } else {
      X_new = [1, p4, p7, 0, 0];  // No arrhythmia (reference)
    }
    
    // Calculate variance of linear predictor: X^T * Cov * X
    let var_eta = 0;
    for (let i = 0; i < 5; i++) {
      for (let j = 0; j < 5; j++) {
        var_eta += X_new[i] * covBeta[i][j] * X_new[j];
      }
    }
    
    const se_eta = Math.sqrt(var_eta);
    
    // Linear predictor (logit scale)
    const map1 = [0, 0.6583, 1.7944];
    const eta = -0.1311 + map1[arr] + 0.00687 * p4 + 0.0232 * p7;
    
    // Critical value for 95% CI (z-score for 97.5th percentile)
    const z_crit = 1.96;
    
    // CI on logit scale
    const eta_lower = eta - z_crit * se_eta;
    const eta_upper = eta + z_crit * se_eta;
    
    // Transform to probability scale
    const prob_lower = 1 / (1 + Math.exp(-eta_lower));
    const prob_upper = 1 / (1 + Math.exp(-eta_upper));
    
    return { 
      prob_lower, 
      prob_upper,
      se_eta,
      eta
    };
  }

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
      
      // Get user-specified Hill coefficient (with default of 1)
      const hillCoeffInput = document.getElementById('hillCoefficient');
      let userHillCoeff = hillCoeffInput ? parseFloat(hillCoeffInput.value) : 1;
      
      // Validate Hill coefficient
      if (isNaN(userHillCoeff) || userHillCoeff <= 0) {
        alert('WARNING\n[Hill Coefficient]\nHill coefficient must be a positive number. Using default value of 1.');
        userHillCoeff = 1;
        if (hillCoeffInput) hillCoeffInput.value = '1';
      }
      
      // Reasonable range check
      if (userHillCoeff < 0.1 || userHillCoeff > 10) {
        if (!confirm(`WARNING\n[Hill Coefficient]\nHill coefficient of ${userHillCoeff} is outside typical range (0.1-10).\n\nDo you want to continue?`)) {
          return;
        }
      }
      
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
      const guess = { Bottom, Top, EC50: med(concs), Hill: userHillCoeff };

      // Hill function with variable Hill coefficient
      const hillf = (x, p) => {
        return p.Bottom + (p.Top - p.Bottom) / (1 + Math.pow(p.EC50 / x, p.Hill));
      };

      // Fit EC50 using grid search with user-specified Hill coefficient
      let best = { ...guess }, minE = Infinity;
      const loss = p => fpdcs.reduce((s,y,i)=>s + Math.pow(hillf(concs[i],p) - y, 2), 0);
      
      // Grid search over EC50 values
      const minConc = Math.min(...concs);
      const maxConc = Math.max(...concs);
      const searchMin = minConc * 0.01;  // Wider search range
      const searchMax = maxConc * 100;
      const numSteps = 200;  // More steps for better precision
      
      for (let i = 0; i <= numSteps; i++) {
        // Logarithmic spacing for EC50 search
        const logMin = Math.log10(searchMin);
        const logMax = Math.log10(searchMax);
        const ec = Math.pow(10, logMin + (i / numSteps) * (logMax - logMin));
        
        const t = { ...guess, EC50: ec };
        const e = loss(t);
        if (e < minE) { 
          minE = e; 
          best = t; 
        }
      }
      
      // Refine search around best EC50 found
      const refinedSearchMin = best.EC50 * 0.5;
      const refinedSearchMax = best.EC50 * 2.0;
      for (let i = 0; i <= 100; i++) {
        const logMin = Math.log10(refinedSearchMin);
        const logMax = Math.log10(refinedSearchMax);
        const ec = Math.pow(10, logMin + (i / 100) * (logMax - logMin));
        
        const t = { ...guess, EC50: ec };
        const e = loss(t);
        if (e < minE) { 
          minE = e; 
          best = t; 
        }
      }

      const FPDc = hillf(Cmax||Math.min(...concs), best);
      p7=FPDc;
      document.getElementById('predictor1').value=String(arr);
      document.getElementById('predictor4').value=isFinite(p4)?Number(p4).toFixed(4):'';
      document.getElementById('predictor7').value=isFinite(p7)?Number(p7).toFixed(4):'';
      validatePredictorRanges();
      
      if(p4===0&&p7===0){}
      const Thr=assay==='30'?Bottom*1.103:Bottom*1.0794;
      const logM=assay==='30'?(Thr+0.35)/0.92:(Thr+0.17)/0.93;
      (()=>{
        const el=document.getElementById('estimatedQTc'); 
        if(el){ 
          el.innerHTML=`<strong>QTc (log M):</strong> ${logM.toFixed(4)}<br><strong>Conc >10ms QT:</strong> ${Math.pow(10,logM).toFixed(4)} µM<br><strong>Hill Coefficient:</strong> ${best.Hill.toFixed(2)}<br><strong>EC50:</strong> ${best.EC50.toFixed(4)} µM`; 
        }
      })();
      
      const fitX=Array.from({length:100},(_,i)=>Math.pow(10,Math.log10(Math.max(0.001,Math.min(...concs)))+i*(Math.log10(Math.max(...concs))-Math.log10(Math.max(0.001,Math.min(...concs))))/99));
      const fitY=fitX.map(x=>hillf(x,best));
      if(hillChart)hillChart.destroy();
      hillChart=new Chart(document.getElementById('hillPlot'),{
        type:'line',
        data:{
          labels:fitX,
          datasets:[
            {label:`Hill Fit (n=${best.Hill.toFixed(2)})`,data:fitX.map((x,i)=>({x,y:fitY[i]})),borderWidth:3,fill:false,borderColor:'rgb(75, 192, 192)'},
            {label:'Data',type:'scatter',data:concs.map((x,i)=>({x,y:fpdcs[i]})),pointRadius:4,backgroundColor:'rgb(54, 162, 235)'},
            {label:'Cmax',type:'scatter',data:[{x:Cmax,y:FPDc}],pointRadius:6,backgroundColor:'rgb(255, 99, 132)'}
          ]
        },
        options:{
          responsive:true,
          maintainAspectRatio:false,
          scales:{
            x:{type:'logarithmic', grid:{lineWidth:5}, ticks:{font:{size:20}}, title:{display:true, text:'Concentration (µM)', font:{size:18}}},
            y:{grid:{lineWidth:5}, ticks:{font:{size:20}}, title:{display:true, text:'ΔΔFPDc or ΔΔAPD90c (ms)', font:{size:18}}}
          },
          plugins: {
            legend: {
              display: true,
              position: 'top'
            },
            title: {
              display: true,
              text: `EC50 = ${best.EC50.toFixed(4)} µM, Hill = ${best.Hill.toFixed(2)}`,
              font: { size: 14 }
            }
          }
        }
      });
    }
    
    // model probabilities (Model 1 only)
    const map1=[0,0.6583,1.7944];
    const logit1=-0.1311+map1[arr]+0.00687*p4+0.0232*p7;
    Prob1=1/(1+Math.exp(-logit1));
    if(Prob1<0){}
    
    // Calculate confidence intervals
    const ciResults = calculateCI(arr, p4, p7);
    window.ciLower = ciResults.prob_lower;
    window.ciUpper = ciResults.prob_upper;
    
    updateModelPanel();
  }

  
  
function updateModelPanel(){
  const title=document.getElementById('modelTitle'),
        sub=document.getElementById('modelSubtitle'),
        res=document.getElementById('modelResults');
  const labels=['High or Intermediate','Low'];
  const data=[Prob1*100,(1-Prob1)*100];
  
  // Get CI values
  const ciLower = window.ciLower !== null ? window.ciLower : Prob1;
  const ciUpper = window.ciUpper !== null ? window.ciUpper : Prob1;
  
  if(title) title.innerText='Model 1 TdP Risk';
  if(sub)   sub.innerHTML='The background calculation uses a logistic regression model. The model outputs are:';
  const colors = labels.map(l => (l.includes('High') ? 'rgb(230,75,53)' : 'rgb(3,160,135)'));
  
  if(res){ 
    res.innerHTML = '<ul style="margin-left:20px;">' +
      `<li><strong>High or Intermediate TdP Risk Probability:</strong> ${data[0].toFixed(1)}% <span style="color:#666;">(95% CI: ${(ciLower*100).toFixed(1)}% – ${(ciUpper*100).toFixed(1)}%)</span></li>` +
      `<li><strong>Low TdP Risk Probability:</strong> ${data[1].toFixed(1)}%</li>` +
      '</ul>'; 
  }
  
  if(modelChart) modelChart.destroy();
  const datasets = labels.map((l,i)=>({ label: l + ' Risk', data: [data[i]], backgroundColor: colors[i], stack: 'risk', borderWidth: 0, borderRadius: 6, maxBarThickness: 72 }));
  
  // Custom plugin to draw error bars
  const errorBarPlugin = {
    id: 'errorBars',
    afterDatasetsDraw(chart) {
      if (window.ciLower === null || window.ciUpper === null) return;
      
      const ctx = chart.ctx;
      const meta = chart.getDatasetMeta(0); // High/Intermediate risk dataset
      const bar = meta.data[0];
      
      if (!bar) return;
      
      const x = bar.x;
      const yScale = chart.scales.y;
      
      // Calculate error bar positions
      const probPercent = Prob1 * 100;
      const lowerPercent = ciLower * 100;
      const upperPercent = ciUpper * 100;
      
      const yCenter = yScale.getPixelForValue(probPercent);
      const yLower = yScale.getPixelForValue(lowerPercent);
      const yUpper = yScale.getPixelForValue(upperPercent);
      
      // Draw error bar
      ctx.save();
      ctx.strokeStyle = 'rgba(0, 0, 0, 0.8)';
      ctx.lineWidth = 2.5;
      
      // Vertical line
      ctx.beginPath();
      ctx.moveTo(x, yLower);
      ctx.lineTo(x, yUpper);
      ctx.stroke();
      
      // Lower cap
      ctx.beginPath();
      ctx.moveTo(x - 8, yLower);
      ctx.lineTo(x + 8, yLower);
      ctx.stroke();
      
      // Upper cap
      ctx.beginPath();
      ctx.moveTo(x - 8, yUpper);
      ctx.lineTo(x + 8, yUpper);
      ctx.stroke();
      
      ctx.restore();
    }
  };
  
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
        tooltip: { 
          callbacks: { 
            label: ctx => {
              if (ctx.datasetIndex === 0 && window.ciLower !== null) {
                return `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(1)}% (95% CI: ${(ciLower*100).toFixed(1)}% – ${(ciUpper*100).toFixed(1)}%)`;
              }
              return `${ctx.dataset.label}: ${ctx.parsed.y.toFixed(1)}%`;
            }
          } 
        },
        annotation: { annotations: { riskThreshold: { type: 'line', yMin: 80, yMax: 80, borderColor: 'red', borderWidth: 2, borderDash:[6,6], label:{ display:true, content:'Risk Probability Threshold (80%)', position:'end', color:'red', font:{ size:14 }, yAdjust:-6 } } } }
      }
    },
    plugins: [errorBarPlugin]
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
