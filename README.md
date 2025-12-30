# Simulation-Causal-Inference-With-Missing-Data

Simulation code for my **bachelor project** on **Missing Data in Causal Inference**.

---

## How to Execute the Code

1. Load all packages in `main.R`.  
2. Run all functions (all R scripts **except** `main.R`).
4. Run the examples provided in `main.R`.

---

## File Overview

### `data.generation.R` :
- `generate_data()` generate treatment, outcome and missingness in a CRE.
- `generate_data()` generate treatment, outcome and missingness in a SRE.
- All missingness mechanisms are in `generate_data()`.

### `Li.R` :
- Use `worst_case_randomization_test()` for worst-case bounding in a CRE.
- `worst_case_randomization_test_two_step()` doesn't work, as `wang_upper_M` doesn't correctly estimate $M$ for some reason I didn't have time to fully investigate.
- Use `worst_case_stratified_randomization_test()` for worst-case bounding in SRE.

### `Heussen.R` :
- Use `conditional_randomization_test()` for the conditioning method in a CRE.

### `Heng.R` :
- Identical to the `iArt` package, run it if you can't download the `iArt` package for some reason.

### `evaluation.R` :
- `compute_pvalue()` computes the p-value for Heussen's, Heng's or Li's methods in a CRE.
- `run_sim()` runs the full simulation grid in a CRE.
- `compute_rejection_rates()` computes the rejection rates of $H_0$, given an output of the simulation grid.
- `run_sim_SRE()` and `compute_rejection_rates_SRE` are similar to above, but in a SRE.
- `compute_pvalue_CRE_on_SRE()` and `run_sim_SRE_as_CRE()` is used to compare SRE results and CRE results by computing the p-values using the CRE methods, but on stratified data.

### `visualization.R` :
- `plot_rejection_rates()` is used to plot the rejection rats of $H_0$ versus the true $\delta$ used to generate the data.
DISCLAIMER : I got help from AI models for better plotting using `ggplot2`.

### `main.R` :
- contains all necessary packages.
- The examples provided in the script are mainly the ones used in my report (but not in chronological order).

---

## Plots Overview (by folder name)

- Heng : Rejection rate using Heng's method.
- Heussen : Rejection rate using Heussen's method.
- Li : Rejection rate using Li's method.
- SRE_aligned : Rejection rate using the SRE aligned method
- SRE_combined : Rejection rate using the SRE combined method.
- Heng_v_Heussen : Compares Heng's and Heussen's methods.
- Li_v_Li : Compares Li's general and monotone methods.
- CRE_vSRE : Compares Bounding methods under a SRE vs ignoring stratification in a CRE.
- randomization_distribution : Plot showing the randomization distribution of the Wilcoxon Rank Sum statistic.
- old_plots : Contains outdated plots, not relevant anymore, but contains similar results.

---
 
### Notes
- Heng's simulation are **computationally heavy**.
- The **two-step procedure for Li assumption 2** in `Li.R` does not work.

