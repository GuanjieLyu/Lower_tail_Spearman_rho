## MISE comparison for lower-tail Spearman's rho estimators
##
## This script compares the empirical-copula estimator and the
## Bernstein-smoothed estimator of the lower-tail Spearman curve
## p -> rho(p). The simulation covers FGM, Gaussian, and Frank copulas.
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
##   Rscript MISE_comparison.R

options(stringsAsFactors = FALSE, warn = 1)

## -------------------------------------------------------------------------
## Simulation settings
## -------------------------------------------------------------------------

family_values <- c("FGM", "Gaussian", "Frank")
theta_values_fgm <- c(-1, -0.5, 0, 0.5, 1)
tau_values <- c(-0.5, -0.25, 0, 0.25, 0.5)
n_values <- c(50L, 200L)
p_grid <- seq(0.05, 1.00, length.out = 20L)

K <- 500L
base_seed <- 20260430L
truth_halton_size <- 300000L
ncores <- if (.Platform$OS.type == "windows") {
  1L
} else {
  max(1L, parallel::detectCores() - 2L)
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- sprintf("mise_comparison_results_%s", timestamp)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

ise_file <- file.path(results_dir, "mise_ise_by_replication.csv")
summary_file <- file.path(results_dir, "mise_summary.csv")
pointwise_file <- file.path(results_dir, "mise_pointwise_mse.csv")
decomposition_file <- file.path(results_dir, "mise_bias_variance_decomposition.csv")
truth_file <- file.path(results_dir, "mise_truth_curves.csv")
raw_file <- file.path(results_dir, "mise_raw_results.rds")
figure_file <- file.path(results_dir, "mise_summary.pdf")

bernstein_degree <- function(n) floor(n^(2 / 3))

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

empirical_rho <- function(u, v, p_grid) {
  n <- length(u)
  rx <- rank(u, ties.method = "first") / n
  ry <- rank(v, ties.method = "first") / n

  vapply(p_grid, function(p) {
    integral_value <- mean(pmax(p - rx, 0) * pmax(p - ry, 0))
    d <- p^3 / 3 - p^4 / 4
    (integral_value - p^4 / 4) / d
  }, numeric(1))
}

rank_copula_grid <- function(u, v, m) {
  n <- length(u)
  rx <- rank(u, ties.method = "first")
  ry <- rank(v, ties.method = "first")

  tab <- matrix(0, n, n)
  tab[cbind(rx, ry)] <- 1

  cum <- apply(apply(tab, 2, cumsum), 1, cumsum)
  cpad <- matrix(0, n + 1, n + 1)
  cpad[-1, -1] <- cum / n

  grid_index <- floor(n * (0:m) / m + 1e-12) + 1L
  cpad[grid_index, grid_index, drop = FALSE]
}

bernstein_weights <- function(p_grid, m) {
  k <- 0:m
  sapply(p_grid, function(p) pbeta(p, k + 1, m - k + 1) / (m + 1))
}

bernstein_rho <- function(u, v, p_grid, m, weights = NULL) {
  cmat <- rank_copula_grid(u, v, m)
  if (is.null(weights)) weights <- bernstein_weights(p_grid, m)

  integral_values <- colSums(weights * (cmat %*% weights))
  d <- p_grid^3 / 3 - p_grid^4 / 4
  (integral_values - p_grid^4 / 4) / d
}

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

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
  cbind(dependence_cases[i, ], n = n_values)
}))
case_grid$case_id <- seq_len(nrow(case_grid))

## -------------------------------------------------------------------------
## Monte Carlo simulation
## -------------------------------------------------------------------------

one_replication <- function(job) {
  set.seed(job$seed)

  uv <- r_copula_sample(job$family, job$n, job$copula_parameter)
  m <- bernstein_degree(job$n)
  weights <- bernstein_weights(job$p_grid, m)

  truth <- job$truth
  rho_empirical <- empirical_rho(uv[, 1], uv[, 2], job$p_grid)
  rho_bernstein <- bernstein_rho(uv[, 1], uv[, 2], job$p_grid, m, weights)

  se_empirical <- (rho_empirical - truth)^2
  se_bernstein <- (rho_bernstein - truth)^2

  common <- data.frame(
    family = job$family,
    dependence_type = job$dependence_type,
    dependence = job$dependence,
    copula_parameter = job$copula_parameter,
    n = job$n,
    m = m,
    rep = job$rep
  )

  ise <- cbind(
    common,
    data.frame(
      ise_empirical = trapz(job$p_grid, se_empirical),
      ise_bernstein = trapz(job$p_grid, se_bernstein)
    )
  )
  ise$aise_empirical <- ise$ise_empirical / diff(range(job$p_grid))
  ise$aise_bernstein <- ise$ise_bernstein / diff(range(job$p_grid))

  pointwise <- cbind(
    data.frame(common[rep(1, length(job$p_grid)), ], row.names = NULL),
    data.frame(
      p = job$p_grid,
      truth = truth,
      rho_empirical = rho_empirical,
      rho_bernstein = rho_bernstein,
      se_empirical = se_empirical,
      se_bernstein = se_bernstein
    )
  )

  list(ise = ise, pointwise = pointwise)
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
      truth = truth_cache[[case_row$truth_id]],
      seed = base_seed + 100000L * case_row$case_id + rep_id
    )
  })

  parallel_apply(jobs, one_replication)
}

case_results <- vector("list", nrow(case_grid))
for (i in seq_len(nrow(case_grid))) {
  case_results[[i]] <- run_case(case_grid[i, ])
}

rep_results <- unlist(case_results, recursive = FALSE)
ise_results <- do.call(rbind, lapply(rep_results, `[[`, "ise"))
pointwise_results <- do.call(rbind, lapply(rep_results, `[[`, "pointwise"))

## -------------------------------------------------------------------------
## Summaries
## -------------------------------------------------------------------------

summarize_case <- function(dat) {
  data.frame(
    family = dat$family[1],
    dependence_type = dat$dependence_type[1],
    dependence = dat$dependence[1],
    copula_parameter = dat$copula_parameter[1],
    n = dat$n[1],
    m = dat$m[1],
    K = nrow(dat),
    mise_empirical = mean(dat$ise_empirical),
    mise_bernstein = mean(dat$ise_bernstein),
    sd_ise_empirical = sd(dat$ise_empirical),
    sd_ise_bernstein = sd(dat$ise_bernstein),
    mean_aise_empirical = mean(dat$aise_empirical),
    mean_aise_bernstein = mean(dat$aise_bernstein)
  )
}

summary_results <- do.call(
  rbind,
  lapply(
    split(ise_results,
          list(ise_results$family, ise_results$dependence, ise_results$n),
          drop = TRUE),
    summarize_case
  )
)

summary_results$mise_reduction_percent <- with(
  summary_results,
  100 * (mise_empirical - mise_bernstein) / mise_empirical
)
summary_results$mise_ratio_bernstein_to_empirical <- with(
  summary_results,
  mise_bernstein / mise_empirical
)

summary_results <- summary_results[
  order(match(summary_results$family, family_values),
        summary_results$n,
        summary_results$dependence),
]

pointwise_mse <- aggregate(
  cbind(se_empirical, se_bernstein) ~
    family + dependence_type + dependence + copula_parameter + n + m + p + truth,
  data = pointwise_results,
  FUN = mean
)

names(pointwise_mse)[names(pointwise_mse) == "se_empirical"] <- "mse_empirical"
names(pointwise_mse)[names(pointwise_mse) == "se_bernstein"] <- "mse_bernstein"
pointwise_mse$mse_reduction_percent <- with(
  pointwise_mse,
  100 * (mse_empirical - mse_bernstein) / mse_empirical
)

pointwise_mse <- pointwise_mse[
  order(match(pointwise_mse$family, family_values),
        pointwise_mse$n,
        pointwise_mse$dependence,
        pointwise_mse$p),
]

decompose_case <- function(dat, estimator) {
  est_col <- if (estimator == "Empirical") "rho_empirical" else "rho_bernstein"

  pointwise_parts <- lapply(split(dat, dat$p), function(dp) {
    est <- dp[[est_col]]
    truth <- dp$truth[1]
    mean_est <- mean(est)

    data.frame(
      p = dp$p[1],
      bias_sq = (mean_est - truth)^2,
      variance = mean((est - mean_est)^2),
      mse = mean((est - truth)^2)
    )
  })

  pointwise_parts <- do.call(rbind, pointwise_parts)
  pointwise_parts <- pointwise_parts[order(pointwise_parts$p), ]

  data.frame(
    family = dat$family[1],
    dependence_type = dat$dependence_type[1],
    dependence = dat$dependence[1],
    copula_parameter = dat$copula_parameter[1],
    n = dat$n[1],
    m = dat$m[1],
    estimator = estimator,
    integrated_bias_sq = trapz(pointwise_parts$p, pointwise_parts$bias_sq),
    integrated_variance = trapz(pointwise_parts$p, pointwise_parts$variance),
    mise = trapz(pointwise_parts$p, pointwise_parts$mse)
  )
}

case_splits <- split(
  pointwise_results,
  list(pointwise_results$family, pointwise_results$dependence, pointwise_results$n),
  drop = TRUE
)

bias_variance_decomposition <- do.call(rbind, lapply(case_splits, function(dat) {
  rbind(decompose_case(dat, "Empirical"), decompose_case(dat, "Bernstein"))
}))

bias_variance_decomposition$bias_sq_percent <- with(
  bias_variance_decomposition,
  100 * integrated_bias_sq / mise
)
bias_variance_decomposition$variance_percent <- with(
  bias_variance_decomposition,
  100 * integrated_variance / mise
)

bias_variance_decomposition <- bias_variance_decomposition[
  order(match(bias_variance_decomposition$family, family_values),
        bias_variance_decomposition$n,
        bias_variance_decomposition$dependence,
        bias_variance_decomposition$estimator),
]

## -------------------------------------------------------------------------
## Figures
## -------------------------------------------------------------------------

safe_name <- function(x) tolower(gsub("[^A-Za-z0-9]+", "_", x))

boxplot_panel_ylim <- function(box_values, cap_prob = 0.98) {
  pooled <- unlist(box_values, use.names = FALSE)
  upper <- as.numeric(quantile(pooled, probs = cap_prob, na.rm = TRUE))

  if (!is.finite(upper) || upper <= 0) {
    upper <- max(pooled, na.rm = TRUE)
  }

  c(0, 1.08 * upper)
}

plot_mise_summary <- function(summary_results, figure_file) {
  pdf(figure_file, width = 11, height = 8, pointsize = 12, bg = "white")
  op <- par(mfrow = c(length(family_values), length(n_values)),
            mar = c(4.8, 4.8, 3.2, 1),
            mgp = c(2.8, 0.85, 0),
            cex.lab = 1.1, cex.axis = 0.95, cex.main = 1.1,
            font.main = 1,
            bg = "white")
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  y_max <- max(summary_results$mise_empirical,
               summary_results$mise_bernstein,
               na.rm = TRUE)

  for (fam in family_values) {
    for (n_i in n_values) {
      dat <- summary_results[summary_results$family == fam &
                               summary_results$n == n_i, ]
      dat <- dat[order(dat$dependence), ]
      mat <- rbind(dat$mise_empirical, dat$mise_bernstein)

      barplot(
        mat,
        beside = TRUE,
        names.arg = dat$dependence,
        col = c("darkorange", "#1f77b4"),
        border = NA,
        ylim = c(0, 1.08 * y_max),
        xlab = if (fam == "FGM") expression(theta) else expression(tau),
        ylab = "MISE",
        main = sprintf("%s, n = %d", fam, n_i)
      )

      legend("topleft",
             legend = c("Empirical", "Bernstein"),
             fill = c("darkorange", "#1f77b4"),
             border = NA,
             bty = "n",
             cex = 0.9)
    }
  }
}

draw_ise_boxplot_panel <- function(dat, family_name, n_value) {
  dat <- dat[order(dat$dependence), ]
  dep_order <- sort(unique(dat$dependence))

  box_values <- list()
  at <- numeric(0)
  box_cols <- character(0)
  pos <- 1

  for (dep_i in dep_order) {
    dat_i <- dat[dat$dependence == dep_i, ]
    box_values[[length(box_values) + 1L]] <- dat_i$ise_empirical
    box_values[[length(box_values) + 1L]] <- dat_i$ise_bernstein
    at <- c(at, pos, pos + 0.35)
    box_cols <- c(box_cols, "darkorange", "#1f77b4")
    pos <- pos + 1.35
  }

  midpoints <- tapply(at, rep(seq_along(dep_order), each = 2), mean)

  boxplot(
    box_values,
    at = at,
    col = box_cols,
    border = "grey20",
    boxwex = 0.28,
    outline = TRUE,
    pch = 16,
    cex = 0.4,
    ylim = boxplot_panel_ylim(box_values),
    xaxt = "n",
    xlab = if (family_name == "FGM") expression(theta) else expression(tau),
    ylab = "ISE",
    main = sprintf("%s copula, n = %d", family_name, n_value)
  )

  axis(1, at = midpoints, labels = dep_order)

  legend("topleft",
         legend = c("Empirical copula estimator", "Bernstein estimator"),
         fill = c("darkorange", "#1f77b4"),
         border = NA,
         bty = "n",
         cex = 1.35)
}

plot_ise_boxplot <- function(ise_results, family_name, n_value, figure_file) {
  dat <- ise_results[ise_results$family == family_name &
                       ise_results$n == n_value, ]

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

  draw_ise_boxplot_panel(dat, family_name, n_value)
}

plot_combined_boxplot <- function(ise_results, n_value, figure_file) {
  pdf(figure_file, width = 18, height = 6.2, pointsize = 12, bg = "white")
  op <- par(mfrow = c(1, length(family_values)),
            mar = c(5.2, 5.8, 4.2, 1.4),
            mgp = c(3.3, 1.0, 0),
            cex.lab = 1.7, cex.axis = 1.4, cex.main = 1.6,
            font.main = 1,
            bg = "white")
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (fam in family_values) {
    dat <- ise_results[ise_results$family == fam & ise_results$n == n_value, ]
    draw_ise_boxplot_panel(dat, fam, n_value)
  }
}

plot_mise_summary(summary_results, figure_file)

boxplot_files <- character(0)
for (fam in family_values) {
  for (n_i in n_values) {
    file_i <- file.path(
      results_dir,
      sprintf("%s_ise_boxplot_n%d.pdf", safe_name(fam), n_i)
    )
    plot_ise_boxplot(ise_results, fam, n_i, file_i)
    boxplot_files <- c(boxplot_files, file_i)
  }
}

combined_boxplot_files <- character(0)
for (n_i in n_values) {
  file_i <- file.path(
    results_dir,
    sprintf("all_families_ise_boxplot_n%d.pdf", n_i)
  )
  plot_combined_boxplot(ise_results, n_i, file_i)
  combined_boxplot_files <- c(combined_boxplot_files, file_i)
}

## -------------------------------------------------------------------------
## Save outputs
## -------------------------------------------------------------------------

write.csv(ise_results, ise_file, row.names = FALSE)
write.csv(summary_results, summary_file, row.names = FALSE)
write.csv(pointwise_mse, pointwise_file, row.names = FALSE)
write.csv(bias_variance_decomposition, decomposition_file, row.names = FALSE)
write.csv(truth_curves, truth_file, row.names = FALSE)

saveRDS(
  list(
    config = list(
      family_values = family_values,
      theta_values_fgm = theta_values_fgm,
      tau_values = tau_values,
      n_values = n_values,
      p_grid = p_grid,
      K = K,
      base_seed = base_seed,
      truth_halton_size = truth_halton_size,
      ncores = ncores
    ),
    dependence_cases = dependence_cases,
    truth_curves = truth_curves,
    ise_results = ise_results,
    summary_results = summary_results,
    pointwise_mse = pointwise_mse,
    bias_variance_decomposition = bias_variance_decomposition,
    pointwise_results = pointwise_results,
    figure_files = list(
      mise_summary = figure_file,
      ise_boxplots = boxplot_files,
      combined_ise_boxplots = combined_boxplot_files
    )
  ),
  raw_file
)

invisible(NULL)
