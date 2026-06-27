################################################################################
############################ Run real data analysis ############################
################################################################################

# library(VIM)
library(ggplot2)
library(patchwork)

source("shared_estimand_helpers.R")
source("method_SubA_all.R")
source("method_SACE.R")
source("method_WhileAlive.R")
source("method_SE_simple.R")

# source("Real Data/RealData_1_process.R")
source("Real Data/RealData_1_readSynthetic.R")
source("Real Data/RealData_2_EDA.R")
source("Real Data/RealData_3_estimators.R")
source("Real Data/RealData_4_naive_SACE_plot.R")

tau_time <- c(0, 1/4, 1/2, 1)
use_L_history <- FALSE
boot_reps <- getOption("realdata_boot_reps", 500)

realdata_results <- run_realdata_analysis(
  data = data.RealData,
  tau_time = tau_time,
  boot_reps = boot_reps,
  use_L_history = use_L_history
)

est_points <- realdata_results$est_points
re_Boot <- realdata_results$re_Boot
result_list <- realdata_results$result_list

realdata_table <- do.call(
  rbind,
  lapply(names(result_list), function(nm) {
    data.frame(
      Estimand = nm,
      result_list[[nm]],
      row.names = NULL
    )
  })
)
realdata_table[sapply(realdata_table, is.numeric)] <-
  lapply(realdata_table[sapply(realdata_table, is.numeric)], round, 2)

## =============================================================================
##    Output: Tables 4 and 5
## =============================================================================

print(realdata_table, row.names = FALSE)

write.csv(
  realdata_table,
  file = "Real Data/realdata_table.csv",
  row.names = FALSE
)


## =============================================================================
##    Output: Figure 7
## =============================================================================

p_survival <- plot_survival_probability(data.RealData, tau_time)
p_sace <- plot_sace_results(data.RealData, result_list, tau_time)

ggsave(
  "Real Data/realdata_survival_probability.png",
  plot = p_survival,
  width = 6,
  height = 4,
  dpi = 300
)

## =============================================================================
##    Output: Figure 8
## =============================================================================
ggsave(
  "Real Data/realdata_combined_plotSACE.png",
  plot = p_sace,
  width = 12,
  height = 5,
  dpi = 300
)

