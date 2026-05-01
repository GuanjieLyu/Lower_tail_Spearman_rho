## Degree selection for the Bernstein lower-tail Spearman estimator
##
## This script studies the sensitivity of the Bernstein-smoothed
## lower-tail Spearman estimator to the polynomial degree m. For each
## copula family, dependence setting, and sample size, it estimates the
## MISE of the Bernstein estimator over a grid of m values and marks the
## rule-of-thumb degree floor(n^(2/3)) in the resulting figures.
##
## The objective is to assess whether floor(n^(2/3)) is a reasonable
## practical choice for the Bernstein degree.
##
## Dependence settings:
##   - FGM:      theta in {-1, -0.5, 0, 0.5, 1}
##   - Gaussian: Kendall's tau in {-0.5, -0.25, 0, 0.25, 0.5}
##   - Frank:    Kendall's tau in {-0.5, -0.25, 0, 0.25, 0.5}
##
## Parallel note:
##   The simulation uses parallel::mclapply(), which is available on
##   macOS, Linux, and WSL/Unix-like environments. It is not suitable for
##   native Windows R sessions, where the script automatically falls back
##   to sequential execution. For parallel execution on Windows, run this
##   script through WSL.
##
## Run:
##   Rscript MISE_degree_selection.R

options(stringsAsFactors = FALSE, warn = 1)

## -------------------------------------------------------------------------
## Simulation settings
## -------------------------------------------------------------------------

family_values <- c("FGM", "Gaussian", "Frank")
theta_values_fgm <- c(-1, -0.5, 0, 0.5, 1)
tau_values <- c(-0.5, -0.25, 0, 0.25, 0.5)
n_values <- c(50L, 200L)
p_grid <- seq(0.05, 1.00, length.out = 20L)
m_grid <- 1:80

K <- 500L
base_seed <- 20260501L
truth_halton_size <- 300000L
ncores <- if (.Platform$OS.type == "windows") {
  1L
} else {
  max(1L, parallel::detectCores() - 2L)
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- sprintf("mise_degree_selection_results_%s", timestamp)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

curve_summary_file <- file.path(results_dir, "mise_degree_curve_summary.csv")
rule_summary_file <- file.path(results_dir, "mise_degree_rule_summary.csv")
truth_file <- file.path(results_dir, "mise_degree_truth_curves.csv")
raw_file <- file.path(results_dir, "mise_degree_raw_ise.rds")
combined_figure_file <- file.path(results_dir, "mise_degree_curves.pdf")

parallel_apply <- function(x, fun) {
  if (ncores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(x, fun, mc.cores = ncores, mc.set.seed = FALSE)
  } else {
    lapply(x, fun)
  }
}

## -------------------------------------------------------------------------
## Copula parameter conversions
## -------------------------------------------------------------------------

debye1 <- function(theta) {
  if (abs(theta) < 1e-8) return(1)

  integrand <- function(t) {
    ifelse(abs(t) < 1e-7, 1 - t / 2 + t^2 / 12, t / expm1(t))
  }

  integral_value <- if (theta > 0) {
    integrate(integrand, lower = 0, upper = theta,
              rel.tol = 1e-10, subdivisions = 200L)$value
  } else {
    -integrate(integrand, lower = theta, upper = 0,
               rel.tol = 1e-10, subdivisions = 200L)$value
  }

  integral_value / theta
}

frank_tau <- function(theta) {
  if (abs(theta) < 1e-8) return(0)
  1 - 4 / theta + 4 * debye1(theta) / theta
}

frank_theta_from_tau <- function(tau) {
  if (abs(tau) < 1e-12) return(0)

  if (tau > 0) {
    uniroot(function(x) frank_tau(x) - tau,
            lower = 1e-6, upper = 60, tol = 1e-10)$root
  } else {
    uniroot(function(x) frank_tau(x) - tau,
            lower = -60, upper = -1e-6, tol = 1e-10)$root
  }
}

gaussian_rho_from_tau <- function(tau) sin(pi * tau / 2)

## -------------------------------------------------------------------------
## Copula simulation
## -------------------------------------------------------------------------

rfgm <- function(n, theta) {
  u <- runif(n)
  z <- runif(n)
  a <- theta * (1 - 2 * u)
  v <- numeric(n)

  near_independence <- abs(a) < 1e-12
  v[near_independence] <- z[near_independence]

  idx <- !near_independence
  disc <- (1 + a[idx])^2 - 4 * a[idx] * z[idx]
  v[idx] <- ((1 + a[idx]) - sqrt(pmax(disc, 0))) / (2 * a[idx])

  cbind(u = u, v = pmin(pmax(v, 0), 1))
}

rgaussian_copula <- function(n, rho) {
  z1 <- rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  cbind(u = pnorm(z1), v = pnorm(z2))
}

rfrank_copula <- function(n, theta) {
  if (abs(theta) < 1e-12) {
    return(cbind(u = runif(n), v = runif(n)))
  }

  u <- runif(n)
  w <- runif(n)
  a <- exp(-theta * u)
  d <- expm1(-theta)
  x <- w * d / (a - w * (a - 1))
  exp_minus_theta_v <- pmax(1 + x, .Machine$double.eps)
  v <- -log(exp_minus_theta_v) / theta

  cbind(u = u, v = pmin(pmax(v, 0), 1))
}

r_copula_sample <- function(family, n, copula_parameter) {
  if (family == "FGM") {
    rfgm(n, copula_parameter)
  } else if (family == "Gaussian") {
    rgaussian_copula(n, copula_parameter)
  } else if (family == "Frank") {
    rfrank_copula(n, copula_parameter)
  } else {
    stop("Unknown copula family: ", family)
  }
}

## -------------------------------------------------------------------------
## True lower-tail Spearman curves
## -------------------------------------------------------------------------

true_rho_fgm <- function(theta, p) {
  d <- p^3 / 3 - p^4 / 4
  a <- p^2 / 2 - p^3 / 3
  theta * a^2 / d
}

van_der_corput <- function(n, base) {
  index <- seq_len(n)
  result <- numeric(n)
  factor <- 1 / base

  while (any(index > 0)) {
    result <- result + factor * (index %% base)
    index <- floor(index / base)
    factor <- factor / base
  }

  result
}

halton_2d <- function(n) {
  cbind(van_der_corput(n, 2L), van_der_corput(n, 3L))
}

truth_rho_halton <- function(family, copula_parameter, p_grid,
                             n_halton = truth_halton_size) {
  if (abs(copula_parameter) < 1e-12 && family != "FGM") {
    return(rep(0, length(p_grid)))
  }

  h <- halton_2d(n_halton)

  if (family == "Gaussian") {
    rho <- copula_parameter
    z1 <- qnorm(h[, 1])
    z2 <- rho * z1 + sqrt(1 - rho^2) * qnorm(h[, 2])
    uv <- cbind(u = pnorm(z1), v = pnorm(z2))
  } else if (family == "Frank") {
    theta <- copula_parameter
    u <- h[, 1]
    w <- h[, 2]
    a <- exp(-theta * u)
    d <- expm1(-theta)
    x <- w * d / (a - w * (a - 1))
    exp_minus_theta_v <- pmax(1 + x, .Machine$double.eps)
    v <- -log(exp_minus_theta_v) / theta
    uv <- cbind(u = u, v = pmin(pmax(v, 0), 1))
  } else {
    stop("Halton truth is only used for Gaussian and Frank copulas.")
  }

  integral_values <- vapply(p_grid, function(p) {
    mean(pmax(p - uv[, 1], 0) * pmax(p - uv[, 2], 0))
  }, numeric(1))

  d <- p_grid^3 / 3 - p_grid^4 / 4
  (integral_values - p_grid^4 / 4) / d
}

## -------------------------------------------------------------------------
## Estimators
## -------------------------------------------------------------------------

empirical_rho_from_ranks <- function(rx, ry, p_grid) {
  n <- length(rx)
  ux <- rx / n
  uy <- ry / n

  vapply(p_grid, function(p) {
    integral_value <- mean(pmax(p - ux, 0) * pmax(p - uy, 0))
    d <- p^3 / 3 - p^4 / 4
    (integral_value - p^4 / 4) / d
  }, numeric(1))
}

rank_copula_cpad <- function(rx, ry) {
  n <- length(rx)
  tab <- matrix(0, n, n)
  tab[cbind(rx, ry)] <- 1

  cum <- apply(apply(tab, 2, cumsum), 1, cumsum)
  cpad <- matrix(0, n + 1, n + 1)
  cpad[-1, -1] <- cum / n
  cpad
}

bernstein_weights <- function(p_grid, m) {
  k <- 0:m
  sapply(p_grid, function(p) pbeta(p, k + 1, m - k + 1) / (m + 1))
}

bernstein_rho_from_cpad <- function(cpad, n, p_grid, m, weights, grid_index) {
  cmat <- cpad[grid_index, grid_index, drop = FALSE]
  integral_values <- colSums(weights * (cmat %*% weights))
  d <- p_grid^3 / 3 - p_grid^4 / 4
  (integral_values - p_grid^4 / 4) / d
}

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

weights_cache <- setNames(
  lapply(m_grid, function(m) bernstein_weights(p_grid, m)),
  as.character(m_grid)
)

grid_index_cache <- setNames(
  lapply(n_values, function(n) {
    setNames(
      lapply(m_grid, function(m) floor(n * (0:m) / m + 1e-12) + 1L),
      as.character(m_grid)
    )
  }),
  as.character(n_values)
)

## -------------------------------------------------------------------------
## Simulation design
## -------------------------------------------------------------------------

fgm_cases <- data.frame(
  family = "FGM",
  dependence_type = "theta",
  dependence = theta_values_fgm,
  copula_parameter = theta_values_fgm
)

gaussian_cases <- data.frame(
  family = "Gaussian",
  dependence_type = "kendall_tau",
  dependence = tau_values,
  copula_parameter = vapply(tau_values, gaussian_rho_from_tau, numeric(1))
)

frank_cases <- data.frame(
  family = "Frank",
  dependence_type = "kendall_tau",
  dependence = tau_values,
  copula_parameter = vapply(tau_values, frank_theta_from_tau, numeric(1))
)

dependence_cases <- rbind(fgm_cases, gaussian_cases, frank_cases)
dependence_cases$truth_id <- seq_len(nrow(dependence_cases))

truth_rows <- vector("list", nrow(dependence_cases))
truth_cache <- vector("list", nrow(dependence_cases))

for (i in seq_len(nrow(dependence_cases))) {
  fam <- dependence_cases$family[i]
  par <- dependence_cases$copula_parameter[i]

  truth <- if (fam == "FGM") {
    true_rho_fgm(par, p_grid)
  } else {
    truth_rho_halton(fam, par, p_grid)
  }

  truth_cache[[i]] <- truth
  truth_rows[[i]] <- data.frame(
    family = fam,
    dependence_type = dependence_cases$dependence_type[i],
    dependence = dependence_cases$dependence[i],
    copula_parameter = par,
    p = p_grid,
    truth = truth
  )
}

truth_curves <- do.call(rbind, truth_rows)

case_grid <- do.call(rbind, lapply(seq_len(nrow(dependence_cases)), function(i) {
  data.frame(
    dependence_cases[rep(i, length(n_values)), ],
    n = n_values,
    row.names = NULL
  )
}))
case_grid$case_id <- seq_len(nrow(case_grid))

## -------------------------------------------------------------------------
## Monte Carlo simulation
## -------------------------------------------------------------------------

one_replication <- function(job) {
  set.seed(job$seed)

  uv <- r_copula_sample(job$family, job$n, job$copula_parameter)
  rx <- rank(uv[, 1], ties.method = "first")
  ry <- rank(uv[, 2], ties.method = "first")
  cpad <- rank_copula_cpad(rx, ry)

  rho_empirical <- empirical_rho_from_ranks(rx, ry, job$p_grid)
  ise_empirical <- trapz(job$p_grid, (rho_empirical - job$truth)^2)

  bernstein_ise <- vapply(job$m_grid, function(m) {
    rho_bernstein <- bernstein_rho_from_cpad(
      cpad = cpad,
      n = job$n,
      p_grid = job$p_grid,
      m = m,
      weights = job$weights_cache[[as.character(m)]],
      grid_index = job$grid_index_cache[[as.character(m)]]
    )
    trapz(job$p_grid, (rho_bernstein - job$truth)^2)
  }, numeric(1))

  data.frame(
    family = job$family,
    dependence_type = job$dependence_type,
    dependence = job$dependence,
    copula_parameter = job$copula_parameter,
    n = job$n,
    rep = job$rep,
    m = job$m_grid,
    ise_bernstein = bernstein_ise,
    ise_empirical = ise_empirical
  )
}

run_case <- function(case_row) {
  jobs <- lapply(seq_len(K), function(rep_id) {
    list(
      family = case_row$family,
      dependence_type = case_row$dependence_type,
      dependence = case_row$dependence,
      copula_parameter = case_row$copula_parameter,
      n = case_row$n,
      rep = rep_id,
      p_grid = p_grid,
      m_grid = m_grid,
      truth = truth_cache[[case_row$truth_id]],
      weights_cache = weights_cache,
      grid_index_cache = grid_index_cache[[as.character(case_row$n)]],
      seed = base_seed + 100000L * case_row$case_id + rep_id
    )
  })

  parallel_apply(jobs, one_replication)
}

case_results <- vector("list", nrow(case_grid))
for (i in seq_len(nrow(case_grid))) {
  case_results[[i]] <- run_case(case_grid[i, ])
}

ise_results <- do.call(rbind, unlist(case_results, recursive = FALSE))

## -------------------------------------------------------------------------
## Summaries
## -------------------------------------------------------------------------

summarize_degree <- function(dat) {
  data.frame(
    family = dat$family[1],
    dependence_type = dat$dependence_type[1],
    dependence = dat$dependence[1],
    copula_parameter = dat$copula_parameter[1],
    n = dat$n[1],
    m = dat$m[1],
    K = length(unique(dat$rep)),
    mise_bernstein = mean(dat$ise_bernstein),
    sd_ise_bernstein = sd(dat$ise_bernstein),
    mise_empirical = mean(dat$ise_empirical),
    sd_ise_empirical = sd(dat$ise_empirical)
  )
}

curve_summary <- do.call(
  rbind,
  lapply(
    split(
      ise_results,
      list(ise_results$family, ise_results$dependence, ise_results$n, ise_results$m),
      drop = TRUE
    ),
    summarize_degree
  )
)

curve_summary$mise_reduction_percent <- with(
  curve_summary,
  100 * (mise_empirical - mise_bernstein) / mise_empirical
)

curve_summary <- curve_summary[
  order(match(curve_summary$family, family_values),
        curve_summary$n,
        curve_summary$dependence,
        curve_summary$m),
]

summarize_rule <- function(dat) {
  dat <- dat[order(dat$m), ]
  m_nominal <- floor(dat$n[1]^(2 / 3))
  rule_idx <- which.min(abs(dat$m - m_nominal))
  oracle_idx <- which.min(dat$mise_bernstein)

  data.frame(
    family = dat$family[1],
    dependence_type = dat$dependence_type[1],
    dependence = dat$dependence[1],
    copula_parameter = dat$copula_parameter[1],
    n = dat$n[1],
    m_rule_nominal = m_nominal,
    m_rule_used = dat$m[rule_idx],
    m_oracle = dat$m[oracle_idx],
    mise_rule = dat$mise_bernstein[rule_idx],
    mise_oracle = dat$mise_bernstein[oracle_idx],
    mise_empirical = dat$mise_empirical[rule_idx],
    rule_excess_percent =
      100 * (dat$mise_bernstein[rule_idx] / dat$mise_bernstein[oracle_idx] - 1),
    rule_reduction_vs_empirical_percent =
      100 * (dat$mise_empirical[rule_idx] - dat$mise_bernstein[rule_idx]) /
      dat$mise_empirical[rule_idx]
  )
}

rule_summary <- do.call(
  rbind,
  lapply(
    split(curve_summary,
          list(curve_summary$family, curve_summary$dependence, curve_summary$n),
          drop = TRUE),
    summarize_rule
  )
)

rule_summary <- rule_summary[
  order(match(rule_summary$family, family_values),
        rule_summary$n,
        rule_summary$dependence),
]

## -------------------------------------------------------------------------
## Figures
## -------------------------------------------------------------------------

safe_name <- function(x) tolower(gsub("[^A-Za-z0-9]+", "_", x))

dependence_label <- function(family_name, dep_values) {
  symbol <- if (family_name == "FGM") quote(theta) else quote(tau)
  as.expression(lapply(dep_values, function(dep) {
    bquote(.(symbol) == .(dep))
  }))
}

draw_degree_panel <- function(curve_summary, family_name, n_value) {
  dat <- curve_summary[curve_summary$family == family_name &
                         curve_summary$n == n_value, ]
  dep_values <- sort(unique(dat$dependence))
  cols <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")
  cols <- cols[seq_along(dep_values)]
  y_max <- max(dat$mise_bernstein, na.rm = TRUE)
  m_rule <- floor(n_value^(2 / 3))

  plot(NA,
       xlim = range(m_grid),
       ylim = c(0, 1.08 * y_max),
       xlab = "m",
       ylab = "MISE",
       main = sprintf("%s copula, n = %d", family_name, n_value))

  for (i in seq_along(dep_values)) {
    dat_i <- dat[dat$dependence == dep_values[i], ]
    lines(dat_i$m, dat_i$mise_bernstein, col = cols[i], lwd = 3)
  }

  abline(v = m_rule, col = "red", lwd = 3, lty = 4)

  legend("topright",
         legend = c(dependence_label(family_name, dep_values),
                    expression(group(lfloor, n^{2/3}, rfloor))),
         col = c(cols, "red"),
         lty = c(rep(1, length(dep_values)), 4),
         lwd = c(rep(3, length(dep_values)), 3),
         bty = "n",
         cex = 1.35)
}

plot_degree_panel <- function(curve_summary, family_name, n_value, figure_file) {
  pdf(figure_file, width = 9, height = 6.2, pointsize = 12, bg = "white")
  op <- par(mar = c(5.2, 5.8, 4.2, 1.4),
            mgp = c(3.3, 1.0, 0),
            cex.lab = 1.7, cex.axis = 1.4, cex.main = 1.6,
            font.main = 1,
            bg = "white")
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  draw_degree_panel(curve_summary, family_name, n_value)
}

plot_combined_degree_panels <- function(curve_summary, figure_file) {
  pdf(figure_file, width = 14, height = 14, pointsize = 12, bg = "white")
  op <- par(mfrow = c(length(family_values), length(n_values)),
            mar = c(5.2, 5.8, 4.2, 1.4),
            mgp = c(3.3, 1.0, 0),
            cex.lab = 1.45, cex.axis = 1.15, cex.main = 1.35,
            font.main = 1,
            bg = "white")
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (fam in family_values) {
    for (n_i in n_values) {
      draw_degree_panel(curve_summary, fam, n_i)
    }
  }
}

plot_combined_degree_panels(curve_summary, combined_figure_file)

panel_figure_files <- character(0)
for (fam in family_values) {
  for (n_i in n_values) {
    figure_file <- file.path(
      results_dir,
      sprintf("%s_mise_degree_curve_n%d.pdf", safe_name(fam), n_i)
    )
    plot_degree_panel(curve_summary, fam, n_i, figure_file)
    panel_figure_files <- c(panel_figure_files, figure_file)
  }
}
