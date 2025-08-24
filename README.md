# Caveats on using Firth’s penalization in the model-based regression standardization for rare diseases

This repository contains simulation code and input data accompanying the manuscript:  
**"Caveats on using Firth’s penalization in the model-based regression standardization for rare diseases"**  
submitted to *Statistics in Medicine*.

## Abstract

Model-based regression standardization, also known as the parametric g-formula, is widely used to estimate
marginal effect measures. However, in rare disease settings, the small number of observed events relative to
the number of covariates can lead to (quasi-)complete separation, resulting in non-convergent estimates in the
regression models. Firth’s penalized likelihood is a common solution to this issue, ensuring finite parameter
estimates even for separated data. While effective for estimating regression coefficients, Firth’s method introduces
bias into model-based regression standardization because of its tendency to shrink the predicted probabilities
to 0.5, leading to discrepancies between the predicted and observed event rates.  

We examined the implications
of applying Firth’s method to model-based regression standardization, illustrating its potential bias through an
empirical study on surgical site infections in orthopedic surgeries. We also proposed two ad hoc corrections (i.e.,
Firth’s logistic regression with intercept correction and added covariate) to mitigate this bias and evaluated these
methods via simulation studies, comparing them with propensity score-based approaches. Finally, we applied the
proposed method to assess the association between SSI, a rare disease, and smoking status in a clinical database of
orthopedic surgeries.

## Repository structure

- `Main.R`: R script to run the simulation.
- `data/`: Contains 18 CSV files with true values for different scenarios.
- `results/`: Stores simulation outputs (not pushed to GitHub, see `.gitignore`).

## Requirements

- **R version:** 4.2.2 (2022-10-31 ucrt) – *Innocent and Trusting*
- **Required R packages:**
  - MASS  
  - logistf  
  - geepack  
  - detectseparation  
  - parallel  
  - pbapply  

Install them with:
```r
install.packages(c("MASS", "logistf", "geepack", "detectseparation", "parallel", "pbapply"))

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/shashibe/statmed.git
   cd statmed
source("Main.R")
