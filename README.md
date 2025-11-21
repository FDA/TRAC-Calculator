# TdP Risk Assessment Categorization (TRAC) Calculator

The TRAC Calculator is an open-source, web-based tool for assessing **torsades de pointes (TdP) proarrhythmic risk** of small molecule compounds. It uses electrophysiological response predictors measured in human induced pluripotent stem cell-derived cardiomyocytes (hiPSC-CMs) as part of the **Comprehensive In Vitro Proarrhythmia Assay (CiPA)** paradigm.

This tool implements a logistic regression–based model trained on large, multisite datasets, providing researchers and regulators with reproducible, human-relevant risk predictions.

---

## Running the Program

You can directly run the program in this package by visiting [the live TRAC Calculator webpage](https://tkfeaster.github.io/TRAC-Calculator/) and follow the instructions on the page.

## Features

- **Three validated predictors**:
  1. Arrhythmia-like event occurrence (categorical)
  2. Maximum repolarization change (ms) observed at any concentration
  3. Repolarization change (ms) at expected *Cmax*
- **Two input modes**:
  - Direct predictor input  
  - *Cmax* interpolation with concentration–response curve fitting (Hill curve with Hill coefficient fixed at 1)
- **Program outputs**:
  - Predicted probability of **Low risk** vs **High or Intermediate risk**
  - Bar chart of risk probabilities
  - Hill-fit concentration–response plot with *Cmax* overlay
- **Context of Use (COU)** and **Limitations** documented with links to key references.

---

## Authors

**Developer:** Marshall Ma (marshma@bu.edu, shendu.ma@hhs.fda.gov)  
**Advisor(s):** TK Feaster (tromondae.feaster@hhs.fda.gov), Ksenia Blinova (ksenia.blinova@hhs.fda.gov)

## Getting Started

The TRAC Calculator is a browser-based app (no installation required). Clone or download this repository and open `index.html` in a web browser.

### Input Instructions
1. Review the [Context of Use](cou.html#coU) and [Limitations of Use](cou.html#limitations).
2. Enter predictor values in the **Predictor Inputs** panel.  
   - Predictor ranges are shown in the input fields upon startup.
3. Optionally, use the **Cmax Interpolation** panel to fit a concentration–FPDc curve and automatically generate predictors.
4. Click **Calculate** to view the outputs.

### Outputs
- Logistic regression model results are displayed under **Program Outputs**.  
- Graphical plots update dynamically:
  - Risk probability bar chart
  - Hill-fit concentration–response curve (if Cmax interpolation is used)

---

## References

- Blinova K. *et al.* (2018). *International Multisite Study of Human iPSC-CMs for Drug Proarrhythmic Potential Assessment.* **Cell Reports**, 24(13): 3582–3592. [doi:10.1016/j.celrep.2018.08.079](https://doi.org/10.1016/j.celrep.2018.08.079)  
- Patel D. *et al.* (2019). *Assessment of Proarrhythmic Potential of Drugs in Optogenetically Paced hiPSC-CMs.* **Toxicological Sciences**, 170(1): 167–179. [doi:10.1093/toxsci/kfz076](https://doi.org/10.1093/toxsci/kfz076)  
- Gintant G. *et al.* (2020). *Repolarization studies using human stem cell-derived cardiomyocytes: Validation studies and best practice recommendations.* **Regul Toxicol Pharmacol**, 117: 104756. [doi:10.1016/j.yrtph.2020.104756](https://doi.org/10.1016/j.yrtph.2020.104756)  

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## DISCLAIMER

*This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Or by authors in their capacities as Cognizance Technologies contract scientist at the FDA. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties. The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.*
