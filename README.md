# Caveats on using Firth’s penalization in the model-based regression standardization for rare diseases

This repository contains simulation code and input data for the manuscript:  
**"Caveats on using Firth’s penalization in the model-based regression standardization for rare diseases"**  
(*submitted to Statistics in Medicine*).

The project investigates statistical methods for handling separation and rare events in logistic regression, focusing on the use of Firth’s penalized likelihood in model-based regression standardization (parametric g-formula).

---

## Repository structure

- `Main.R` : Main simulation script  
- `data/` : Contains 18 CSV files with true values for different scenarios  
- `results/` : Output directory for simulation results (excluded from GitHub via `.gitignore`)  

---

## Abstract (Summary)

Model-based regression standardization, also known as the parametric g-formula, is widely used to estimate marginal effect measures.  
However, in rare disease settings, the small number of observed events relative to the number of covariates can lead to (quasi-)complete separation, resulting in non-convergent estimates.  
Firth’s penalized likelihood is a common solution, ensuring finite estimates even for separated data.  

While effective for regression coefficients, Firth’s method introduces bias into regression standardization by shrinking predicted probabilities toward 0.5, leading to discrepancies between predicted and observed event rates.  
We examined this bias through an empirical study on surgical site infections in orthopedic surgeries, and proposed two ad hoc corrections (Firth with intercept correction and added covariate).  
These were evaluated via simulation studies and compared with propensity score–based approaches.  
Finally, we applied the proposed method to a clinical database of orthopedic surgeries.

---

## Requirements

This project was developed in **R version 4.2.2 (2022-10-31 ucrt)**.  
The following R packages are required:

- `MASS`  
- `logistf`  
- `geepack`  
- `detectseparation`  
- `parallel`  
- `pbapply`  

You can install them with:

```r
install.packages(c("MASS", "logistf", "geepack", "detectseparation", "parallel", "pbapply"))


## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/shashibe/statmed.git
   cd statmed
source("Main.R")
