# benchmark botplots
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

loglik_data = read.csv('../opt_results_glike_true_rep/benchmarks.csv')
loglik_data

library(ggplot2)
library(tidyr)
library(dplyr)


plot_data <- loglik_data |>
  pivot_longer(
    cols = c(logp, true),
    names_to = "type",
    values_to = "value"
  )

ggplot(plot_data, aes(x = type, y = value, fill = type)) +
  geom_boxplot() +
  geom_jitter(
    color='black',
    width = 0.0,
    size = 2,
    alpha = 0.8
  ) +  
  facet_wrap(
    ~ NUM_TREES,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c("logp" = "#1f77b4", "true" = "#ff7f0e"),
    labels = c("logp" = "logp_x", "true" = "logp_true")
  ) +
  labs(
    x = NULL,
    y = "Log Probability",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )

param_data = read.csv('../opt_results_glike_true_rep/params.csv')
param_data

plot_param <- function(param, hline = NULL) {

  plot_data <- param_data |>
    select(NUM_TREES, all_of(param)) |>
    rename(value = all_of(param))

  p <- ggplot(plot_data, aes(x = "", y = value)) +
    geom_boxplot() +
    geom_jitter(
      color = "black",
      width = 0.0,
      size = 2,
      alpha = 0.8
    ) +
    facet_wrap(
      ~ NUM_TREES,
      nrow = 1
    ) +
    labs(
      x = NULL,
      y = param
    ) +
    theme_classic()

  if (!is.null(hline)) {
    p <- p +
      geom_hline(
        yintercept = hline,
        linetype = "dashed",
        color = "red"
      )
  }

  return(p)
}
#         x_true = {'t1':19, 't2':411, 't3':1040, 't4':2004, 'r1':0.0, 
#                   'r2':0.198, 'r3':0.334, 'N_admix':35682, 'N_afr':10000, 
#                   'N_eur':13388, 'N_asia':25234, 'N_pol':15695, 'N_aa':2702, 
#                  'N_ooa':2470, 'N_anc':2665, 'gr':0.078}

plot_param("t1", hline=19) + ggtitle('Parameter = t1')
plot_param("t2", hline=411) + ggtitle('Parameter = t2')
plot_param("t3", hline=1040) + ggtitle('Parameter = t3')
plot_param("t4", hline=2004) + ggtitle('Parameter = t4')
plot_param("r1", hline=0.0) + ggtitle('Parameter = r1')
plot_param("r2", hline=0.198) + ggtitle('Parameter = r2')
plot_param("r3", hline=0.334) + ggtitle('Parameter = r3')
plot_param("N_admix", hline=35682) + ggtitle('Parameter = n_admix')
plot_param("N_afr", hline=10000) + ggtitle('Parameter = N_afr')
plot_param("N_eur", hline=13388) + ggtitle('Parameter = N_eur')
plot_param("N_asia", hline=25234) + ggtitle('Parameter = N_asia')
plot_param("N_pol", hline=15695) + ggtitle('Parameter = N_pol')
plot_param("N_aa", hline=2702) + ggtitle('Parameter = N_aa')
plot_param("N_ooa", hline=2470) + ggtitle('Parameter = N_ooa')
plot_param("N_anc", hline=2665) + ggtitle('Parameter = N_anc')
plot_param("gr", hline=0.078) + ggtitle('Parameter = gr')

