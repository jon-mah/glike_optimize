setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

## gLike (3G09)

true_values_3G09 <- c(
  t1 = 848,
  t2 = 5600,
  t3 = 8800,
  N_anc = 7300,
  N_yri = 12300,
  N_ooa = 2100,
  N_ceu = 1000,
  N_chb = 510,
  gr_ceu = 0.004,
  gr_chb = 0.0055
)

CMA_ES_3G09_benchmarks <- read.csv(
  "../opt_results_CMA_ES_3G09_no_m1/benchmarks.csv"
)

CMA_ES_3G09_params <- read.csv(
  "../opt_results_CMA_ES_3G09_no_m1/params.csv"
)


# ============================================================
# Parameter error plot
# ============================================================

# ============================================================
# Parameter error plot
# ============================================================

# Initial parameter values
x0 <- c(
  t1 = 10,
  t2 = 50,
  t3 = 100,
  N_anc = 10000,
  N_yri = 10000,
  N_ooa = 10000,
  N_ceu = 10000,
  N_chb = 10000,
  gr_ceu = 0.1,
  gr_chb = 0.1
)


# CMA-ES parameter data
param_data_CMA_ES <- CMA_ES_3G09_params |>
  mutate(
    method = "CMA-ES"
  )


# CMA-ES results
plot_data_CMA_ES <- param_data_CMA_ES |>
  pivot_longer(
    cols = c(
      t1, t2, t3,
      N_anc, N_yri, N_ooa, N_ceu, N_chb,
      gr_ceu, gr_chb
    ),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    true_value = true_values_3G09[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )


# Initial values
plot_data_x0 <- data.frame(
  parameter = names(x0),
  value = as.numeric(x0),
  method = "CMA-ES"
) |>
  mutate(
    true_value = true_values_3G09[parameter],
    percent_error = ((value - true_value) / true_value) * 100
  )


# Plot
ggplot(
  plot_data_CMA_ES,
  aes(
    x = method,
    y = percent_error,
    fill = method
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +
  geom_point(
    data = plot_data_x0,
    aes(
      x = method,
      y = percent_error
    ),
    inherit.aes = FALSE,
    shape = 12,
    size = 4,
    color = "blue"
  ) +
  facet_wrap(
    ~ parameter,
    nrow = 1
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
  ggtitle("CMA-ES (3G09)")


# ============================================================
# Runtime / likelihood benchmark plot
# ============================================================

true_logp_3G09 <- CMA_ES_3G09_benchmarks$true[1]


# Convert runtime strings such as
# "4h 10m 44.3s" to seconds
runtime_to_seconds <- function(x) {

  h <- as.numeric(
    str_extract(x, "\\d+(?:\\.\\d+)?(?=h)")
  )

  m <- as.numeric(
    str_extract(x, "\\d+(?:\\.\\d+)?(?=m)")
  )

  s <- as.numeric(
    str_extract(x, "\\d+(?:\\.\\d+)?(?=s)")
  )

  h[is.na(h)] <- 0
  m[is.na(m)] <- 0
  s[is.na(s)] <- 0

  h * 3600 + m * 60 + s
}


# First modify the benchmark data, then pivot it
benchmark_data <- CMA_ES_3G09_benchmarks |>
  mutate(
    Method = "CMA-ES (3G09)",
    runtime = runtime_to_seconds(runtime)
  ) |>
  pivot_longer(
    cols = c(runtime, logp),
    names_to = "metric",
    values_to = "value"
  )


# Plot runtime and likelihood
ggplot(
  benchmark_data,
  aes(
    x = Method,
    y = value,
    fill = Method
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    labeller = as_labeller(
      c(
        runtime = "Runtime in seconds",
        logp = "Log likelihood"
      )
    )
  ) +
  geom_hline(
    data = data.frame(
      metric = "logp",
      true_logp_3G09 = true_logp_3G09
    ),
    aes(
      yintercept = true_logp_3G09
    ),
    color = "red",
    linetype = "dashed"
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
  theme_classic() +
  theme(
    legend.position = "right"
  ) +
  ggtitle("Runtime and log likelihood benchmarking") +
  labs(
    y = NULL
  )

### Chronological
gLike_3G09_1 = read.csv('../Analysis/gLike_3G09_chronological_1.csv')
gLike_3G09_2 = read.csv('../Analysis/gLike_3G09_chronological_2.csv')
gLike_3G09_3 = read.csv('../Analysis/gLike_3G09_chronological_3.csv')
gLike_3G09_4 = read.csv('../Analysis/gLike_3G09_chronological_4.csv')
gLike_3G09_5 = read.csv('../Analysis/gLike_3G09_chronological_5.csv')
gLike_3G09_6 = read.csv('../Analysis/gLike_3G09_chronological_6.csv')
gLike_3G09_7 = read.csv('../Analysis/gLike_3G09_chronological_7.csv')
gLike_3G09_8 = read.csv('../Analysis/gLike_3G09_chronological_8.csv')
gLike_3G09_9 = read.csv('../Analysis/gLike_3G09_chronological_9.csv')
gLike_3G09_10 = read.csv('../Analysis/gLike_3G09_chronological_10.csv')

gLike_3G09_all <- bind_rows(
  mget(
    paste0("gLike_3G09_", 1:10)
  ),
  .id = "replicate"
) |>
  mutate(
    replicate = as.integer(replicate),
    iteration = row_number(),
    .by = replicate
  )

# ------------------------------------------------------------
# Parameters to plot
# ------------------------------------------------------------

parameters_3G09 <- c(
  "t1", "t2", "t3",
  "N_anc", "N_yri", "N_ooa", "N_ceu", "N_chb",
  "gr_ceu", "gr_chb"
)

# ------------------------------------------------------------
# Convert parameters to long format
# ------------------------------------------------------------

plot_data_3G09 <- gLike_3G09_all |>
  pivot_longer(
    cols = any_of(parameters_3G09),
    names_to = "parameter",
    values_to = "value"
  )

ggplot(
  gLike_3G09_all,
  aes(
    x = iteration,
    y = log_likelihood,
    group = replicate
  )
) +
  geom_line() +
  geom_point() +
  scale_x_continuous(
    breaks = 1:10
  ) +
  labs(
    x = "Iteration",
    y = "Likelihood"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  ggtitle("gLike (3G09): Likelihood Trajectories, kappa = 50000")

true_values_3G09 <- c(
  t1 = 848,
  t2 = 5600,
  t3 = 8800,
  N_anc = 7300,
  N_yri = 12300,
  N_ooa = 2100,
  N_ceu = 1000,
  N_chb = 510,
  gr_ceu = 0.004,
  gr_chb = 0.0055
)

parameter_plots <- list()

for (param in parameters_3G09) {

  plot_data <- gLike_3G09_all |>
    select(replicate, iteration, all_of(param)) |>
    rename(value = all_of(param))

  parameter_plots[[param]] <- ggplot(
    plot_data,
    aes(
      x = iteration,
      y = value,
      group = replicate,
      color = replicate
    )
  ) +
    geom_hline(
      yintercept = true_values_3G09[[param]],
      linetype = "dashed",
      color = "black"
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = 1:10
    ) +
    labs(
      title = paste("gLike (3G09), kappa = 50000:", param),
      x = "Iteration",
      y = param,
      color = "Replicate"
    ) +
    theme_classic()
}

for (param in parameters_3G09) {
  print(parameter_plots[[param]])
}
