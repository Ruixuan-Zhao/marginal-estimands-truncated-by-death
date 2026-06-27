# Replication Files for “Causal Inference for All: Marginal Estimands for Outcomes Truncated by Death”

This supplementary file contains R code for the simulation studies and real data analysis for the manuscript. 

## Software Requirements

The code uses R and the following R packages:

- `ggplot2`
- `VIM` (not required for the simulated example; used only for the real-data preprocessing code)
- `patchwork`

## R Files

### Method Files

- `shared_estimand_helpers.R`: shared helper functions.
- `method_SubA_all.R`: estimators for the while guaranteed-survival and while extended-survival estimands.
- `method_SE_simple.R`: estimators for the marginal separable effect estimands.
- `method_WhileAlive.R`: estimators for while-alive estimands.
- `method_SACE.R`: estimator for SACE(t) using the method of Wang (2017).

### Simulation Files

- `Simul_SubA.R`: generates simulated data for the while guaranteed-survival and while extended-survival estimands.
- `Simul_trueValue_SubA.R`: computes truth values for `mu(1)-mu(0)` and `mu.ext`.
- `run_simul_SubA.R`: runs the simulation for `mu(1)-mu(0)` and `mu.ext`.
- `Simul_SE.R`: generates simulated data for marginal separable effect estimands.
- `Simul_trueValue_SE.R`: computes truth values for the marginal separable effect contrasts.
- `run_Simul_SE.R`: runs the simulation for the marginal separable effect contrasts.

### Real Data Files

- `Real Data/RealData_1_readSynthetic.R`: reads the simulated example dataset as `data.RealData`.
- `Real Data/RealData_2_EDA.R`: defines the survival probability plot function for Figure 7.
- `Real Data/RealData_3_estimators.R`: defines estimator functions for real data analysis, including point estimation and bootstrap.
- `Real Data/RealData_4_naive_SACE_plot.R`: defines the SACE(t) plot function.
- `Real Data/RealData_5_run_all.R`: runs the real data analysis and outputs the real data tables and figures.

### Running the Synthetic Example for Figures 7 - 8 and Tables 4 - 5

Since our real data is not publicly available, we provide the replication code based on a simulated dataset (`example_realdata.rds`). Because the simulated dataset differs from our real data, the resulting estimates and figures are not expected to match those reported in the manuscript.

The file for running the synthetic example corresponding to Figures 7-8 and Tables 4-5 is:

- `Real Data/RealData_5_run_all.R`

Run the script from inside the `code` directory:

```r
source("Real Data/RealData_5_run_all.R")
```

Tables 4 and 5 can be read from the `realdata_table` object created by the script. The script also saves:

- `Real Data/realdata_survival_probability.png`: Figure 7.
- `Real Data/realdata_combined_plotSACE.png`: Figure 8.

### Running the Simulations for Supplementary Table S4

The file for reproducing Supplementary Table S4 is:

- `run_results_TableS4.R`

Run the script from inside the `code` directory:

```r
source("run_results_TableS4.R")
```

Tables S4.1 - S4.4 can be constructed from the `sim_res_SubA` and `sim_res_SE` objects created by the script.



#### Use of AI tools

The code was developed and reviewed by the authors. OpenAI ChatGPT and Codex were used to assist with code organization and documentation.
