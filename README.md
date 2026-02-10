# Causal Machine Learning–Guided Personalized Immunochemotherapy Strategies in Intrahepatic Cholangiocarcinoma (CMICC)

## Abstract
### Background & Aims

Immunochemotherapy (IO-chemo) has become the first-line standard care for patients with advanced intrahepatic cholangiocarcinoma (iCCA), but the incremental benefit of adding IO varies greatly among individuals, making it urgently to develop systems for identifying high-benefit populations
 based on individualized treatment effect (ITE) estimation.
 
### Methods
This multi-cohort study included patients with iCCA who underwent either IO-chemo or chemo. The discovery cohort (2018–2022) included patients from three centers, while an independent cohort (2017–2023) from seven other centers served for external validation. Target trial emulation was employed to obtain unbiased average treatment effect estimation. We developed a causal machine-learning (ML) model (CMICC) to estimate heterogeneous treatment effects for IO-chemo. Based on predicted ITE, patients were stratified into high-benefit, no-to-moderate-benefit, and negative-benefit groups for each treatment. Counterfactual analyses compared overall survival between model-guided treatment selection and treatment actually received. The Qini and TOC curves were used to evaluate model performance. Model interpretability was assessed using SHAP.

### Results 
A total of 1485 patients, 696 in the IO-chemo group and 789 in the chemo group, were included in the discovery cohort, and 562 patients, 246 in the IO-chemo group and 316 in the chemo group were included in the external validation cohort. Using a recursive feature addition approach, we selected 17 of 55 multidimensional variables to build the CMICC model in the discovery cohort. In the high-benefit group, compared with chemo, the hazard ratio (HR) for IO-chemo was 0.39 (95% CI, 0.30, 0.52), P<0.001, with mortality rates reduced by 24.1%, 31.2%, and 28.5% at 12, 24, and 36 months, respectively. In the no-to-moderate-benefit groups, IO–chemo did not significantly differ from chemotherapy (HR, 0.91; 95% CI, 0.70, 1.18; P=0.488). In the negative-benefit groups, IO–chemo was associated with worse outcomes than chemo (HR, 1.91; 95% CI, 1.47, 2.48; P<0.001), with 12-, 24-, and 36-month mortality increased by 29.9%, 11.0%, and 5.1%, respectively. Counterfactual analysis indicated that CMICC-guided treatment decisions could improve survival compared with actually received treatment (RecomIC-RecIC vs RecomIC-RecC: HR: 0.53 [95%CI, 0.43, 0.65], P<0.001; RecomC-RecC vs RecomC-RecIC: HR: 0.62 [95%CI, 0.50, 0.77], P<0.001). These results were independently confirmed in external validation cohort.

### Conclusion
This study demonstrated that the causal ML-based CMICC model estimated ITE in iCCA patients receiving IO-chemo, performing effective benefit stratification, and may improve survival through individualized, CMICC-guided treatment decisions. 
