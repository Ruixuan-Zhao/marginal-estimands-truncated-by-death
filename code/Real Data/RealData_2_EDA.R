################################################################################
############################ Figure 7 function ################################
################################################################################

plot_survival_probability <- function(data, tau_time) {
  t_max <- length(tau_time) - 1

  df_surv <- data.frame(
    Time = rep(tau_time, 2),
    Value = c(
      1,
      sapply(seq_len(t_max), function(t) {
        mean(data[[paste0("S", t)]][data$Z == 0])
      }),
      1,
      sapply(seq_len(t_max), function(t) {
        mean(data[[paste0("S", t)]][data$Z == 1])
      })
    ),
    Treatment = rep(c("MP", "DE"), each = length(tau_time))
  )

  ggplot(df_surv, aes(x = Time, y = Value, color = Treatment)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = tau_time) +
    scale_y_continuous(limits = c(0.6, 1)) +
    labs(
      x = expression(tau[t] ~ "(year)"),
      y = "Survival probability",
      color = "Treatment"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
}
