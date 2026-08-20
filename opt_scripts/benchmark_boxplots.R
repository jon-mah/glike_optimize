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

init_true = read.csv('../opt_results_glike_true_rep/params.csv')
init_arbitrary = read.csv('../opt_results_glike_arbitrary_rep/params.csv')

plot_param <- function(param, true_param = NULL) {

  # Combine the two datasets
  plot_data <- bind_rows(
    init_true |>
      select(NUM_TREES, all_of(param)) |>
      rename(value = all_of(param)) |>
      mutate(initialization = "init_true"),

    init_arbitrary |>
      select(NUM_TREES, all_of(param)) |>
      rename(value = all_of(param)) |>
      mutate(initialization = "init_arbitrary")
  ) |>
    mutate(
      percent_error = ((value - true_param) / true_param) * 100
    )

  ggplot(
    plot_data,
    aes(
      x = initialization,
      y = percent_error,
      fill = initialization
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.8)
    ) +
    facet_wrap(
      ~ NUM_TREES,
      nrow = 1
    ) +
    scale_fill_manual(
      values = c(
        "init_true" = "#1f77b4",
        "init_arbitrary" = "#ff7f0e"
      ),
      labels = c(
        "init_true" = "True",
        "init_arbitrary" = "Arbitrary"
      )
    ) +
    labs(
      x = NULL,
      y = "Percent Error (%)",
      fill = NULL,
      title = 'Initial values'
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "red"
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) + 
    scale_y_continuous(limits = NULL)
}
#         x_true = {'t1':19, 't2':411, 't3':1040, 't4':2004, 'r1':0.0, 
#                   'r2':0.198, 'r3':0.334, 'N_admix':35682, 'N_afr':10000, 
#                   'N_eur':13388, 'N_asia':25234, 'N_pol':15695, 'N_aa':2702, 
#                  'N_ooa':2470, 'N_anc':2665, 'gr':0.078}

plot_param("t1", true_param=19) + ggtitle('Parameter = t1')
plot_param("t2", true_param=411) + ggtitle('Parameter = t2')
plot_param("t3", true_param=1040) + ggtitle('Parameter = t3')
plot_param("t4", true_param=2004) + ggtitle('Parameter = t4')
plot_param("r1", true_param=0.01) + ggtitle('Parameter = r1, NOTE: TRUE VALUE = 0, set to 0.01 for plot')
plot_param("r2", true_param=0.198) + ggtitle('Parameter = r2')
plot_param("r3", true_param=0.334) + ggtitle('Parameter = r3')
plot_param("N_admix", true_param=35682) + ggtitle('Parameter = N_admix')
plot_param("N_afr", true_param=10000) + ggtitle('Parameter = N_afr')
plot_param("N_eur", true_param=13388) + ggtitle('Parameter = N_eur')
plot_param("N_asia", true_param=25234) + ggtitle('Parameter = N_asia')
plot_param("N_pol", true_param=15695) + ggtitle('Parameter = N_pol')
plot_param("N_aa", true_param=2702) + ggtitle('Parameter = N_aa')
plot_param("N_ooa", true_param=2470) + ggtitle('Parameter = N_ooa')
plot_param("N_anc", true_param=2665) + ggtitle('Parameter = N_anc')
plot_param("gr", true_param=0.078) + ggtitle('Parameter = gr')

