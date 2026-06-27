################################################################################
################################ Figure 8 ######################################
################################################################################

plot_sace_results <- function(data, result_list, tau_time) {
  t_max <- length(tau_time) - 1
  non_covariate_cols <- c(
    "Z", "SURVT",
    paste0("Y", 0:t_max),
    paste0("S", seq_len(t_max)),
    paste0("L", seq_len(t_max))
  )
  covariate_cols <- setdiff(names(data), non_covariate_cols)

  Prob_St_Z0 <- compute_St_Zz(
    data = data,
    t_max = t_max,
    covariate_cols = covariate_cols
  )$Prob_St_Z0

  df_subgroup <- data.frame(
    Time = tau_time[-1],
    Value = Prob_St_Z0[-1]
  )

  p_subgroup <- ggplot(df_subgroup, aes(x = Time, y = Value)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(size = 2, color = "steelblue") +
    labs(
      x = expression(tau[t] ~ "(year)"),
      y = "Subgroup probability"
    ) +
    ylim(0.5, 1) +
    scale_x_continuous(breaks = tau_time[-1]) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    )

  df_sace <- data.frame(
    Time = tau_time[-1],
    Value = result_list$results.SACE$Point_Estimate,
    lower = result_list$results.SACE$CI_Lower,
    upper = result_list$results.SACE$CI_Upper
  )

  p_sace <- ggplot(df_sace, aes(x = Time, y = Value)) +
    geom_point(size = 2, color = "steelblue") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, color = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      x = expression(tau[t] ~ "(year)"),
      y = "SACE(t) with 95% CI"
    ) +
    scale_x_continuous(breaks = tau_time[-1]) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    )

  p_sace <- p_sace + theme(plot.margin = margin(5.5, 30, 5.5, 5.5))
  p_subgroup <- p_subgroup + theme(plot.margin = margin(5.5, 5.5, 5.5, 30))

  p_sace + p_subgroup + plot_layout(ncol = 2)
}
