################################################################################
############################## Simulation results ##############################
################################### Table S4 ###################################
################################################################################

source("shared_estimand_helpers.R")
source("Simul_SubA.R")
source("Simul_trueValue_SubA.R")
source("method_SubA_all.R")
source("run_simul_SubA.R")
source("Simul_SE.R")
source("Simul_trueValue_SE.R")
source("method_SE_simple.R")
source("run_Simul_SE.R")

## =============================================================================
##    Output: Tables S4.1 and S4.2
## =============================================================================
sim_res_SubA <- run_sim_SubA_grid(
  n.sample.vec = c(500, 2000),
  MC = 500,
  tau_time = c(0, 1/4, 1/2, 1),
  seed = 20260515,
  n.truth = 1000000,
  verbose = TRUE
)
print(sim_res_SubA)

## =============================================================================
##    Output: Tables S4.3 and S4.4
## =============================================================================
sim_res_SE <- run_sim_SE_grid(
  n.sample.vec = c(500, 2000),
  MC = 500,
  tau_time = c(0, 1/4, 1/2, 1),
  seed = 20260515,
  n.truth = 1000000,
  verbose = TRUE
)
print(sim_res_SE)

