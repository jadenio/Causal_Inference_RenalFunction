# Causal Inference - Effect of Perfluoroalkyl Substances on Renal Function
## Bayesian Additive Regression Trees
Bayesian Additive Regression Trees (BART) is a sum-of-trees regression model that shows great performance for causal inference tasks.

An example of a 3-tree BART model is shown below

![Screenshot](https://i.imgur.com/A6L4CeK.png)

## Simulation Study
To compare the peformance of BART to other methods when the treatment is continuous, a simulation study is performed. BART is compared to Generalized Additive Models (GAM), Inverse Probability of Treatment Weighting (IPTW), Targeted Maximum Likelihood Estimation (TMLE), and Generalized Propensity Scores (GPS) methods. 

Three simulation scenarios are created, each with three levels of model misspecification: no misspecification, moderate misspecification, strong misspecification.

The results of the third simulation scenario are shown in the figure below. BART clearly performs best in the majority of situations.
![Screenshot](https://i.imgur.com/zlpa7lU.png)

## Causal Effect of Perfluoroalkyl acids (PFAS) on Renal Function
### Data
The data for this causal analysis comes from the National Health and Nutrition Examination Survey (NHANES) conducted by the Centers for Disease Control and Prevention (CDC) in the United States.

Bayesian Additive Regression Trees (BART) is used to derive the Average Dose-Response Function (ADRF) between various Perfluoroalykl substances (PFAS) and kidney function, as measured by estimated glomerular filtration rate (eGFR). 

The PFAS of interest include Total PFAS (the sum of all PFAS), PFOS, PFOA, PFHxS, and PFNA. Each PFAS is measured by serum blood concentration in ng/mL.

The Average Dose-Response Functions for the overall sample (not adjusted by gender)

![Screenshot](https://i.imgur.com/ByXWNKE.png)

The Average Dose-Response Functions adjusted by gender

![Screenshot](https://i.imgur.com/zRXuPhF.png)
