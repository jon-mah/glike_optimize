setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(tidyverse)

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

CMA_kappa_10000_benchmarks = read.csv(
  '../opt_results_CMA_ES_3G09_no_m1_kappa_10000/benchmarks.csv'
)

CMA_kappa_10000_params = read.csv(
  '../opt_results_CMA_ES_3G09_no_m1_kappa_10000/params.csv'
)

CMA_kappa_25000_benchmarks = read.csv(
  '../opt_results_CMA_ES_3G09_no_m1_kappa_25000/benchmarks.csv'
)

CMA_kappa_25000_params = read.csv(
  '../opt_results_CMA_ES_3G09_no_m1_kappa_25000/params.csv'
)

CMA_kappa_50000_benchmarks = read.csv(
  '../opt_results_CMA_ES_3G09_no_m1_kappa_50000/benchmarks.csv'
)

CMA_kappa_50000_params = read.csv(
  '../opt_results_CMA_ES_3G09_no_m1_kappa_50000/params.csv'
)

# ============================================================
# Parameter error plot
# ============================================================

# ============================================================
# Parameter value plot
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


# Parameters to plot
parameters <- c(
  "t1", "t2", "t3",
  "N_anc", "N_yri", "N_ooa", "N_ceu", "N_chb",
  "gr_ceu", "gr_chb"
)


# ------------------------------------------------------------
# CMA-ES results
# ------------------------------------------------------------

plot_data_CMA_ES <- CMA_ES_3G09_params |>
  pivot_longer(
    cols = all_of(parameters),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(
    method = "CMA-ES",
    true_value = true_values_3G09[parameter]
  )


# ------------------------------------------------------------
# Initial values
# ------------------------------------------------------------

plot_data_x0 <- data.frame(
  parameter = names(x0),
  x0 = as.numeric(x0)
) |>
  mutate(
    true_value = true_values_3G09[parameter]
  )


# ------------------------------------------------------------
# True values
# ------------------------------------------------------------

true_values_plot <- data.frame(
  parameter = parameters,
  true_value = true_values_3G09[parameters]
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

ggplot(
  plot_data_CMA_ES,
  aes(
    x = method,
    y = value,
    fill = method
  )
) +

  # CMA-ES boxplots
  geom_boxplot(
    outlier.shape = NA
  ) +

  # Individual CMA-ES replicates
  geom_jitter(
    width = 0.1,
    size = 1.8,
    alpha = 0.8
  ) +

  # Initial x0 value
  geom_point(
    data = plot_data_x0,
    aes(
      x = "CMA-ES",
      y = x0
    ),
    inherit.aes = FALSE,
    shape = 18,
    size = 4
  ) +

  # True parameter value
  geom_hline(
    data = true_values_plot,
    aes(
      yintercept = true_value
    ),
    color = "red",
    linetype = "dashed"
  ) +

  facet_wrap(
    ~ parameter,
    nrow = 1,
    scales = "free_y"
  ) +

  labs(
    x = NULL,
    y = "Parameter Value",
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

### Kappa tests

convert_runtime <- function(x) {
    h <- as.numeric(sub(".*([0-9]+)h.*", "\\1", x))
    m <- as.numeric(sub(".*h ([0-9]+)m.*", "\\1", x))
    s <- as.numeric(sub(".*m ([0-9.]+)s.*", "\\1", x))

    h * 3600 + m * 60 + s
}


# Combine into one data frame
# Read benchmark files
benchmarks_10000 <- read_csv(
    '../opt_results_CMA_ES_3G09_no_m1_kappa_10000/benchmarks.csv',
    show_col_types = FALSE
) %>%
    mutate(kappa = 10000)

benchmarks_25000 <- read_csv(
    '../opt_results_CMA_ES_3G09_no_m1_kappa_25000/benchmarks.csv',
    show_col_types = FALSE
) %>%
    mutate(kappa = 25000)

benchmarks_50000 <- read_csv(
    '../opt_results_CMA_ES_3G09_no_m1_kappa_50000/benchmarks.csv',
    show_col_types = FALSE
) %>%
    mutate(kappa = 50000)

# Combine into one data frame
benchmarks <- bind_rows(
    benchmarks_10000,
    benchmarks_25000,
    benchmarks_50000
)

benchmarks <- benchmarks %>%
    mutate(runtime_seconds = convert_runtime(runtime))

# Treat kappa as a categorical variable
benchmarks$kappa <- factor(
    benchmarks$kappa,
    levels = c(10000, 25000, 50000),
    labels = c("κ = 10,000", "κ = 25,000", "κ = 50,000")
)

ggplot(benchmarks, aes(x = kappa, y = runtime_seconds)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        width = 0.1,
        alpha = 0.5
    ) +
    labs(
        x = "κ",
        y = "Runtime (seconds)",
        title = "CMA-ES Runtime across κ Values"
    ) +
    theme_classic()

true_logp <- benchmarks$true[1]

ggplot(benchmarks, aes(x = kappa, y = logp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        width = 0.1,
        alpha = 0.5
    ) +
    geom_hline(
      yintercept = true_logp,
      linetype = "dashed",
      color = 'red'
    ) +
    labs(
        x = "κ",
        y = "Log-likelihood",
        title = "CMA-ES Log-likelihood across κ Values"
    ) +
    theme_classic()

# Read parameter files
params_10000 <- read_csv(
    '../opt_results_CMA_ES_3G09_no_m1_kappa_10000/params.csv',
    show_col_types = FALSE
) %>%
    mutate(kappa = 10000)

params_25000 <- read_csv(
    '../opt_results_CMA_ES_3G09_no_m1_kappa_25000/params.csv',
    show_col_types = FALSE
) %>%
    mutate(kappa = 25000)

params_50000 <- read_csv(
    '../opt_results_CMA_ES_3G09_no_m1_kappa_50000/params.csv',
    show_col_types = FALSE
) %>%
    mutate(kappa = 50000)

# Combine files
params <- bind_rows(
    params_10000,
    params_25000,
    params_50000
)

# Make kappa categorical
params$kappa <- factor(
    params$kappa,
    levels = c(10000, 25000, 50000),
    labels = c("κ = 10,000", "κ = 25,000", "κ = 50,000")
)

parameter_names <- c(
    "N_anc",
    "N_ceu",
    "N_chb",
    "N_ooa",
    "N_yri",
    "gr_ceu",
    "gr_chb",
    "t1",
    "t2",
    "t3"
)

plots <- list()

for (parameter in parameter_names) {

    plots[[parameter]] <- ggplot(
        params,
        aes(
            x = kappa,
            y = .data[[parameter]]
        )
    ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(
            width = 0.1,
            alpha = 0.5
        ) +
        geom_hline(
          yintercept = true_values_3G09[parameter],
          linetype = "dashed",
          color='red'
        ) +
        labs(
            x = "κ",
            y = parameter,
            title = paste("CMA-ES estimates of", parameter)
        ) +
        theme_classic()
}
plots[1]
plots[2]
plots[3]
plots[4]
plots[5]
plots[6]
plots[7]
plots[8]
plots[9]
plots[10]
