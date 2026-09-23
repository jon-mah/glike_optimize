# benchmark botplots
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# loglik_data = read.csv('../opt_results_glike_true_rep/benchmarks.csv')
# loglik_data
# 
# plot_data <- loglik_data |>
#   pivot_longer(
#     cols = c(logp, true),
#     names_to = "type",
#     values_to = "value"
#   )
# 
# ggplot(plot_data, aes(x = type, y = value, fill = type)) +
#   geom_boxplot() +
#   geom_jitter(
#     color='black',
#     width = 0.0,
#     size = 2,
#     alpha = 0.8
#   ) +  
#   facet_wrap(
#     ~ NUM_TREES,
#     nrow = 1,
#     scales = "free_y"
#   ) +
#   scale_fill_manual(
#     values = c("logp" = "#1f77b4", "true" = "#ff7f0e"),
#     labels = c("logp" = "logp_x", "true" = "logp_true")
#   ) +
#   labs(
#     x = NULL,
#     y = "Log Probability",
#     fill = NULL
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "right"
#   )
# 
# init_true = read.csv('../opt_results_glike_true_rep/params.csv')
# init_arbitrary = read.csv('../opt_results_glike_arbitrary_rep/params.csv')
# 
# plot_param <- function(param, true_param = NULL) {
# 
#   # Combine the two datasets
#   plot_data <- bind_rows(
#     init_true |>
#       select(NUM_TREES, all_of(param)) |>
#       rename(value = all_of(param)) |>
#       mutate(initialization = "init_true"),
# 
#     init_arbitrary |>
#       select(NUM_TREES, all_of(param)) |>
#       rename(value = all_of(param)) |>
#       mutate(initialization = "init_arbitrary")
#   ) |>
#     mutate(
#       percent_error = ((value - true_param) / true_param) * 100
#     )
# 
#   ggplot(
#     plot_data,
#     aes(
#       x = initialization,
#       y = percent_error,
#       fill = initialization
#     )
#   ) +
#     geom_boxplot(
#       position = position_dodge(width = 0.8)
#     ) +
#     facet_wrap(
#       ~ NUM_TREES,
#       nrow = 1
#     ) +
#     scale_fill_manual(
#       values = c(
#         "init_true" = "#1f77b4",
#         "init_arbitrary" = "#ff7f0e"
#       ),
#       labels = c(
#         "init_true" = "True",
#         "init_arbitrary" = "Arbitrary"
#       )
#     ) +
#     labs(
#       x = NULL,
#       y = "Percent Error (%)",
#       fill = NULL,
#       title = 'Initial values'
#     ) +
#     geom_hline(
#       yintercept = 0,
#       linetype = "dashed",
#       color = "red"
#     ) +
#     theme_classic() +
#     theme(
#       legend.position = "right",
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank()
#     ) + 
#     scale_y_continuous(limits = NULL)
# }
#         x_true = {'t1':19, 't2':411, 't3':1040, 't4':2004, 'r1':0.0, 
#                   'r2':0.198, 'r3':0.334, 'N_admix':35682, 'N_afr':10000, 
#                   'N_eur':13388, 'N_asia':25234, 'N_pol':15695, 'N_aa':2702, 
#                  'N_ooa':2470, 'N_anc':2665, 'gr':0.078}

# plot_param("t1", true_param=19) + ggtitle('Parameter = t1')
# plot_param("t2", true_param=411) + ggtitle('Parameter = t2')
# plot_param("t3", true_param=1040) + ggtitle('Parameter = t3')
# plot_param("t4", true_param=2004) + ggtitle('Parameter = t4')
# plot_param("r1", true_param=0.01) + ggtitle('Parameter = r1, NOTE: TRUE VALUE = 0, set to 0.01 for plot')
# plot_param("r2", true_param=0.198) + ggtitle('Parameter = r2')
# plot_param("r3", true_param=0.334) + ggtitle('Parameter = r3')
# plot_param("N_admix", true_param=35682) + ggtitle('Parameter = N_admix')
# plot_param("N_afr", true_param=10000) + ggtitle('Parameter = N_afr')
# plot_param("N_eur", true_param=13388) + ggtitle('Parameter = N_eur')
# plot_param("N_asia", true_param=25234) + ggtitle('Parameter = N_asia')
# plot_param("N_pol", true_param=15695) + ggtitle('Parameter = N_pol')
# plot_param("N_aa", true_param=2702) + ggtitle('Parameter = N_aa')
# plot_param("N_ooa", true_param=2470) + ggtitle('Parameter = N_ooa')
# plot_param("N_anc", true_param=2665) + ggtitle('Parameter = N_anc')
# plot_param("gr", true_param=0.078) + ggtitle('Parameter = gr')

## 3G09
# gLike_3G09_benchmarks = read.csv('../opt_results_glike_3G09/benchmarks.csv')
# gLike_3G09_params = read.csv('../opt_results_glike_3G09/params.csv')
# CMA_ES_3G09_benchmarks = read.csv('../opt_results_CMA_ES_3G09/benchmarks.csv')
# CMA_ES_3G09_params = read.csv('../opt_results_CMA_ES_3G09/params.csv')

gLike_3G09_benchmarks = read.csv('../opt_results_glike_3G09_no_m1/benchmarks.csv')
gLike_3G09_params = read.csv('../opt_results_glike_3G09_no_m1/params.csv')
CMA_ES_3G09_benchmarks = read.csv('../opt_results_CMA_ES_3G09_no_m1/benchmarks.csv')
CMA_ES_3G09_params = read.csv('../opt_results_CMA_ES_3G09_no_m1/params.csv')

gLike_NH_benchmarks = read.csv('../opt_results_glike_NH/benchmarks.csv')
gLike_NH_params = read.csv('../opt_results_glike_NH/params.csv')
CMA_ES_NH_benchmarks = read.csv('../opt_results_CMA_ES_NH/benchmarks.csv')
CMA_ES_NH_params = read.csv('../opt_results_CMA_ES_NH/params.csv')


true_logp_3G09 = gLike_3G09_benchmarks$true[1]
true_logp_NH = CMA_ES_NH_benchmarks$true[1]

# Convert runtime strings such as "4h 10m 44.3s" to seconds
runtime_to_seconds <- function(x) {
  h <- as.numeric(str_extract(x, "\\d+(?:\\.\\d+)?(?=h)"))
  m <- as.numeric(str_extract(x, "\\d+(?:\\.\\d+)?(?=m)"))
  s <- as.numeric(str_extract(x, "\\d+(?:\\.\\d+)?(?=s)"))
  
  h[is.na(h)] <- 0
  m[is.na(m)] <- 0
  s[is.na(s)] <- 0
  
  h * 3600 + m * 60 + s
}

benchmark_data <- bind_rows(
  gLike_3G09_benchmarks |>
    mutate(
      Method = "gLike (3G09)",
      runtime = runtime_to_seconds(runtime)
    ),
  
  CMA_ES_3G09_benchmarks |>
    mutate(
      Method = "CMA-ES (3G09)",
      runtime = runtime_to_seconds(runtime)
    ),
  gLike_NH_benchmarks |>
    mutate(
      Method = "gLike (NH)",
      runtime = runtime_to_seconds(runtime)
    ),
  CMA_ES_NH_benchmarks |>
    mutate(
      Method = "CMA-ES (NH)",
      runtime = runtime_to_seconds(runtime)
    )
) |>
  pivot_longer(
    cols = c(runtime, logp),
    names_to = "metric",
    values_to = "value"
  )

# ggplot(
#   benchmark_data,
#   aes(x = Method, y = value, fill = Method)
# ) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(
#     width = 0.1,
#     size = 2,
#     alpha = 0.8
#   ) +
#   facet_wrap(
#     ~ metric,
#     scales = "free_y",
#     labeller = as_labeller(c(
#       runtime = "Runtime in seconds",
#       logp = "Log likelihood"
#     ))
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "right"
#   ) +
#   ggtitle('Runtime') +
#   labs(y = NULL)
  

ggplot(
  benchmark_data,
  aes(x = Method, y = value, fill = Method)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.1,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    labeller = as_labeller(c(
      runtime = "Runtime in seconds",
      logp = "Log likelihood"
    ))
  ) +
  geom_hline(
    data = data.frame(
      metric = "logp",
      true_logp_3G09 = true_logp_3G09
    ),
    aes(yintercept = true_logp_3G09),
    color = "red",
    linetype = "dashed",
  ) +
  geom_text(
    data = data.frame(
      metric = "logp",
      true_logp_3G09 = true_logp_3G09
    ),
    aes(
      x = Inf,
      y = true_logp_3G09,
      label = "True log likelihood (3G09)"
    ),
    color = "red",
    hjust = 1.05,
    vjust = -0.5,
    inherit.aes = FALSE
  ) +
  geom_hline(
    data = data.frame(
      metric = "logp",
      true_logp_NH = true_logp_NH
    ),
    aes(yintercept = true_logp_NH),
    color = "blue",
    linetype = "dashed",
  ) +
  geom_text(
    data = data.frame(
      metric = "logp",
      true_logp_NH = true_logp_NH
    ),
    aes(
      x = Inf,
      y = true_logp_NH - 1200,
      label = "True log likelihood (NH)"
    ),
    color = "blue",
    hjust = 1.05,
    vjust = -0.5,
    inherit.aes = FALSE
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  ) +
  ggtitle('Runtime and log likelihood benchmarking') +
  labs(y = NULL)

## 3G09

# Parameter fit
# Specify the true values for each parameter
true_values_3G09 <- c(
  N_ceu = 10000,
  N_nea = 10000,
  N_yri = 10000,
  m1    = 0.029,
  t1    = 30,
  t2    = 50,
  t3    = 73.95,
  t4    = 290
)

# Combine the two datasets
param_data <- bind_rows(
  gLike_3G09_params |>
    mutate(method = "gLike (3G09)"),
  
  CMA_ES_3G09_params |>
    mutate(method = "CMA-ES (3G09)")
)

# Reshape parameters into long format and calculate percent error
plot_data <- param_data |>
  pivot_longer(
    cols = c(N_ceu, N_nea, N_yri, m1, t1, t2, t3, t4),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    true_value = true_values_3G09[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )

# Plot
ggplot(
  plot_data,
  aes(x = method, y = percent_error, fill = method)
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    ncol = 4,
    scales = "free_y"
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = "Percent Error (%)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank()
  ) +
  ggtitle('3G09 Demographic scenario, no migration')

## NH

# Parameter fit
# Specify the true values for each parameter
true_values_NH <- c(
  t1 = 19,
  t2 = 411, 
  t3 = 1040,
  t4 = 2004,
  r1 = 0.0,
  r2 = 0.198,
  r3 = 0.334,
  N_admix = 35682,
  N_afr = 10000,
  N_eur = 13388,
  N_asia = 25234,
  N_pol = 15695,
  N_aa = 2702,
  N_ooa = 2470,
  N_anc = 2665,
  gr = 0.078
)

# Combine the two datasets
param_data <- bind_rows(
  gLike_NH_params |>
    mutate(method = "gLike (NH)"),
  
  CMA_ES_NH_params |>
    mutate(method = "CMA-ES (NH)")
)

# Reshape parameters into long format and calculate percent error
plot_data <- param_data |>
  pivot_longer(
    cols = c(t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, N_pol, N_aa, N_ooa, N_anc, gr),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    true_value = true_values_NH[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )

# Plot
ggplot(
  plot_data,
  aes(x = method, y = percent_error, fill = method)
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    ncol = 4,
    scales = "free_y"
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = "Percent Error (%)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank()
  ) +
  ggtitle('NH demographic scenario')

## gLike (3G09)

param_data_gLike <- bind_rows(
  gLike_3G09_params |>
    mutate(method = "gLike (3G09)"),
)

### No migration

plot_data_gLike <- param_data_gLike |>
  pivot_longer(
    cols = c(N_ceu, N_nea, N_yri, t1, t2, t3, t4),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    true_value = true_values_3G09[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )

# Plot
ggplot(
  plot_data_gLike,
  aes(x = method, y = percent_error, fill = method)
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    nrow = 1,
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = "Percent Error (%)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank()
  ) +
  ggtitle('gLike (3G09)')

## CMA-ES (3G09)

param_data_CMA_ES <- bind_rows(
  CMA_ES_3G09_params |>
    mutate(method = "CMA_ES (3G09)"),
)

plot_data_CMA_ES <- param_data_CMA_ES |>
  pivot_longer(
    cols = c(N_ceu, N_nea, N_yri, t1, t2, t3, t4),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    true_value = true_values_3G09[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )

# Plot
ggplot(
  plot_data_CMA_ES,
  aes(x = method, y = percent_error, fill = method)
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    nrow = 1,
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = "Percent Error (%)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank()
  ) +
  ggtitle('CMA_ES (3G09)')

## gLike (NH)

param_data_gLike <- bind_rows(
  gLike_NH_params |>
    mutate(method = "gLike (NH)"),
)

plot_data_gLike <- param_data_gLike |>
  pivot_longer(
    cols = c(t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, N_pol, N_aa, N_ooa, N_anc, gr),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    true_value = true_values_NH[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )

# Plot
ggplot(
  plot_data_gLike,
  aes(x = method, y = percent_error, fill = method)
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    nrow = 1,
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = "Percent Error (%)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank()
  ) +
  ggtitle('gLike (NH)')

## CMA-ES (NH)

param_data_CMA_ES <- bind_rows(
  CMA_ES_NH_params |>
    mutate(method = "CMA_ES (NH)"),
)

plot_data_CMA_ES <- param_data_CMA_ES |>
  pivot_longer(
    cols = c(t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, N_pol, N_aa, N_ooa, N_anc, gr),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    true_value = true_values_NH[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )

# Plot
ggplot(
  plot_data_CMA_ES,
  aes(x = method, y = percent_error, fill = method)
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    nrow = 1,
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = "Percent Error (%)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank()
  ) +
  ggtitle('CMA_ES (NH)')
