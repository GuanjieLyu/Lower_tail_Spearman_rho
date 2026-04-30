## FGM confidence intervals and stable-window threshold selection
##
## This script conducts the Monte Carlo experiment for the Bernstein-smoothed
## lower-tail Spearman's rho estimator under the Farlie-Gumbel-Morgenstern
## (FGM) copula. For each Monte Carlo sample, it computes pointwise smooth
## Bernstein bootstrap confidence intervals over a grid of lower-tail
## thresholds and then applies the precision-constrained stable-window
## threshold selector.
##
## Parallel computing note:
##   This script uses parallel::mclapply(), which relies on process forking.
##   It is intended for macOS, Linux, and WSL/Unix-like environments. It is
##   not suitable for native Windows R, where forking is unavailable.
##
## Run:
##   Rscript fgm_ci_threshold_selector_simulation.R
##

options(stringsAsFactors = FALSE, warn = 1)

if (.Platform$OS.type == "windows") {
  stop(
    "This script uses parallel::mclapply() and is intended for macOS, ",
    "Linux, or WSL/Unix-like environments, not native Windows R."
  )
}

## -------------------------------------------------------------------------
## Simulation settings
## -------------------------------------------------------------------------

theta_values <- c(-1, -0.5, 0, 0.5, 1)
n_values <- c(50L, 200L)
p_grid <- seq(0.05, 1.00, length.out = 20L)

K <- 500L
B_boot <- 399L
alpha <- 0.05
base_seed <- 20260430L
ncores <- max(1L, parallel::detectCores() - 2L)

## Bernstein degree used in the paper's simulation.
bernstein_degree <- function(n) floor(n^(2 / 3))

## Stable-window selector tuning parameters.
selector_tuning <- list(
  L = 3L,
  lambda = 1,
  p_min = 0.10,
  h_max = 0.20,
  n_min = 5L
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- sprintf("fgm_ci_threshold_selector_results_%s", timestamp)
figure_dir <- file.path(results_dir, "figures_pdf")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

raw_file <- file.path(results_dir, "fgm_ci_raw_results.rds")
ci_summary_file <- file.path(results_dir, "fgm_ci_summary.csv")
selector_file <- file.path(results_dir, "fgm_threshold_selector_selected.csv")
selector_summary_file <- file.path(results_dir, "fgm_threshold_selector_summary.csv")

cat("FGM Bernstein CI and threshold-selector simulation\n")
cat(sprintf("K = %d Monte Carlo replications\n", K))
cat(sprintf("B = %d bootstrap replications\n", B_boot))
cat(sprintf("alpha = %.3f\n", alpha))
cat(sprintf("ncores = %d\n", ncores))
cat("theta values:", paste(theta_values, collapse = ", "), "\n")
cat("n values:", paste(n_values, collapse = ", "), "\n")
cat("p grid:", paste(sprintf("%.2f", p_grid), collapse = ", "), "\n\n")

## -------------------------------------------------------------------------
## Copula simulation and target functional
## -------------------------------------------------------------------------

rfgm <- function(n, theta) {
  u <- runif(n)
  z <- runif(n)
  a <- theta * (1 - 2 * u)

  v <- numeric(n)
  near_ind <- abs(a) < 1e-12
  v[near_ind] <- z[near_ind]

  disc <- (1 + a[!near_ind])^2 - 4 * a[!near_ind] * z[!near_ind]
  v[!near_ind] <- ((1 + a[!near_ind]) - sqrt(pmax(disc, 0))) /
    (2 * a[!near_ind])

  cbind(u = u, v = pmin(pmax(v, 0), 1))
}

true_rho_fgm <- function(theta, p) {
  d <- p^3 / 3 - p^4 / 4
  a <- p^2 / 2 - p^3 / 3
  theta * a^2 / d
}

## -------------------------------------------------------------------------
## Bernstein lower-tail Spearman estimator
## -------------------------------------------------------------------------

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

lower_tail_counts <- function(u, v, p_grid) {
  n <- length(u)
  rx <- rank(u, ties.method = "first")
  ry <- rank(v, ties.method = "first")

  vapply(p_grid, function(p) {
    cutoff <- floor(n * p + 1e-12)
    sum(rx <= cutoff & ry <= cutoff)
  }, numeric(1))
}

## -------------------------------------------------------------------------
## Smooth Bernstein bootstrap
## -------------------------------------------------------------------------

smooth_boot_sample <- function(u, v, m_boot) {
  n <- length(u)
  rx <- rank(u, ties.method = "first")
  ry <- rank(v, ties.method = "first")

  ii <- sample.int(n, n, replace = TRUE)
  ax <- ceiling(m_boot * rx[ii] / n)
  ay <- ceiling(m_boot * ry[ii] / n)
  ax <- pmin(pmax(ax, 1L), m_boot)
  ay <- pmin(pmax(ay, 1L), m_boot)

  cbind(
    u = rbeta(n, ax, m_boot + 1L - ax),
    v = rbeta(n, ay, m_boot + 1L - ay)
  )
}

bootstrap_ci <- function(u, v, p_grid, m, B, alpha) {
  n <- length(u)
  weights <- bernstein_weights(p_grid, m)
  rho_hat <- bernstein_rho(u, v, p_grid, m, weights)

  boot <- matrix(NA_real_, nrow = B, ncol = length(p_grid))
  for (b in seq_len(B)) {
    wb <- smooth_boot_sample(u, v, m_boot = m)
    boot[b, ] <- bernstein_rho(wb[, 1], wb[, 2], p_grid, m, weights)
  }

  t_abs <- abs(sqrt(n) * sweep(boot, 2, rho_hat, "-"))
  cutoff <- apply(t_abs, 2, quantile, probs = 1 - alpha,
                  type = 8, names = FALSE, na.rm = TRUE)
  half_width <- cutoff / sqrt(n)

  list(
    rho_hat = rho_hat,
    lower = rho_hat - half_width,
    upper = rho_hat + half_width,
    width = 2 * half_width
  )
}

## -------------------------------------------------------------------------
## Precision-constrained stable-window selector
## -------------------------------------------------------------------------

stable_window_selector <- function(p_grid, rho_hat, half_width, n_lower,
                                   L, lambda, p_min, h_max, n_min) {
  J <- length(p_grid)
  admissible <- which(
    p_grid >= p_min &
      half_width <= h_max &
      n_lower >= n_min
  )

  stable <- integer(0)
  for (j in seq_len(J - L + 1L)) {
    idx <- j:(j + L - 1L)
    if (!all(idx %in% admissible)) next

    omega <- max(abs(outer(rho_hat[idx], rho_hat[idx], "-")))
    uncert <- max(outer(half_width[idx], half_width[idx], "+"))
    if (omega <= lambda * uncert) stable <- c(stable, j)
  }

  if (length(stable) == 0L) {
    return(data.frame(
      selected = FALSE,
      j_stab = NA_integer_,
      p_stab = NA_real_,
      omega_stab = NA_real_,
      uncert_stab = NA_real_
    ))
  }

  j_stab <- min(stable)
  idx <- j_stab:(j_stab + L - 1L)

  data.frame(
    selected = TRUE,
    j_stab = j_stab,
    p_stab = p_grid[j_stab],
    omega_stab = max(abs(outer(rho_hat[idx], rho_hat[idx], "-"))),
    uncert_stab = max(outer(half_width[idx], half_width[idx], "+"))
  )
}

## -------------------------------------------------------------------------
## One Monte Carlo replication and case runner
## -------------------------------------------------------------------------

one_replication <- function(job) {
  set.seed(job$seed)

  uv <- rfgm(job$n, job$theta)
  m <- bernstein_degree(job$n)
  truth <- true_rho_fgm(job$theta, job$p_grid)
  ci <- bootstrap_ci(uv[, 1], uv[, 2], job$p_grid, m,
                     B = job$B_boot, alpha = job$alpha)
  n_lower <- lower_tail_counts(uv[, 1], uv[, 2], job$p_grid)

  half_width <- ci$width / 2
  selector <- stable_window_selector(
    p_grid = job$p_grid,
    rho_hat = ci$rho_hat,
    half_width = half_width,
    n_lower = n_lower,
    L = job$selector_tuning$L,
    lambda = job$selector_tuning$lambda,
    p_min = job$selector_tuning$p_min,
    h_max = job$selector_tuning$h_max,
    n_min = job$selector_tuning$n_min
  )

  raw <- data.frame(
    theta = job$theta,
    n = job$n,
    m = m,
    rep = job$rep,
    p = job$p_grid,
    truth = truth,
    rho_hat = ci$rho_hat,
    lower = ci$lower,
    upper = ci$upper,
    width = ci$width,
    half_width = half_width,
    n_lower = n_lower,
    cover = as.integer(ci$lower <= truth & truth <= ci$upper),
    lower_miss = as.integer(truth < ci$lower),
    upper_miss = as.integer(truth > ci$upper)
  )

  if (isTRUE(selector$selected)) {
    j <- selector$j_stab
    selector$rho_stab <- raw$rho_hat[j]
    selector$truth_stab <- raw$truth[j]
    selector$lower_stab <- raw$lower[j]
    selector$upper_stab <- raw$upper[j]
    selector$half_width_stab <- raw$half_width[j]
    selector$n_lower_stab <- raw$n_lower[j]
    selector$selected_cover <- raw$cover[j]
    selector$selected_abs_error <- abs(raw$rho_hat[j] - raw$truth[j])
  } else {
    selector$rho_stab <- NA_real_
    selector$truth_stab <- NA_real_
    selector$lower_stab <- NA_real_
    selector$upper_stab <- NA_real_
    selector$half_width_stab <- NA_real_
    selector$n_lower_stab <- NA_real_
    selector$selected_cover <- NA_real_
    selector$selected_abs_error <- NA_real_
  }

  selector <- cbind(
    theta = job$theta,
    n = job$n,
    m = m,
    rep = job$rep,
    L = job$selector_tuning$L,
    lambda = job$selector_tuning$lambda,
    p_min = job$selector_tuning$p_min,
    h_max = job$selector_tuning$h_max,
    n_min = job$selector_tuning$n_min,
    selector
  )

  list(raw = raw, selector = selector)
}

run_case <- function(theta, n, case_id) {
  m <- bernstein_degree(n)
  cat(sprintf("\nCase %02d: theta = %s, n = %d, m = %d\n",
              case_id, theta, n, m))

  jobs <- lapply(seq_len(K), function(rep_id) {
    list(
      theta = theta,
      n = n,
      rep = rep_id,
      p_grid = p_grid,
      B_boot = B_boot,
      alpha = alpha,
      selector_tuning = selector_tuning,
      seed = base_seed + 100000L * case_id + rep_id
    )
  })

  wrapped <- function(job) {
    ans <- one_replication(job)
    if (job$rep %% max(1L, floor(K / 10L)) == 0L) {
      cat(sprintf("[%s] theta = %s, n = %d, replicate %d/%d finished\n",
                  format(Sys.time(), "%H:%M:%S"),
                  job$theta, job$n, job$rep, K))
      flush.console()
    }
    ans
  }

  parallel::mclapply(jobs, wrapped, mc.cores = ncores, mc.set.seed = FALSE)
}

## -------------------------------------------------------------------------
## Run simulation
## -------------------------------------------------------------------------

case_grid <- expand.grid(theta = theta_values, n = n_values)
case_grid$case_id <- seq_len(nrow(case_grid))

case_results <- vector("list", nrow(case_grid))
for (i in seq_len(nrow(case_grid))) {
  case_results[[i]] <- run_case(
    theta = case_grid$theta[i],
    n = case_grid$n[i],
    case_id = case_grid$case_id[i]
  )
}

rep_results <- unlist(case_results, recursive = FALSE)
raw_results <- do.call(rbind, lapply(rep_results, `[[`, "raw"))
selector_results <- do.call(rbind, lapply(rep_results, `[[`, "selector"))
selector_results$selected <- as.logical(selector_results$selected)

## -------------------------------------------------------------------------
## Summaries
## -------------------------------------------------------------------------

ci_summary <- aggregate(
  cbind(rho_hat, lower, upper, width, half_width, n_lower,
        cover, lower_miss, upper_miss) ~ theta + n + m + p + truth,
  data = raw_results,
  FUN = mean
)

names(ci_summary)[names(ci_summary) == "rho_hat"] <- "mean_rho_hat"
names(ci_summary)[names(ci_summary) == "lower"] <- "mean_lower"
names(ci_summary)[names(ci_summary) == "upper"] <- "mean_upper"
names(ci_summary)[names(ci_summary) == "width"] <- "mean_width"
names(ci_summary)[names(ci_summary) == "half_width"] <- "mean_half_width"
names(ci_summary)[names(ci_summary) == "n_lower"] <- "mean_n_lower"
names(ci_summary)[names(ci_summary) == "cover"] <- "coverage"
names(ci_summary)[names(ci_summary) == "lower_miss"] <- "lower_miss_rate"
names(ci_summary)[names(ci_summary) == "upper_miss"] <- "upper_miss_rate"

ci_summary$bias <- ci_summary$mean_rho_hat - ci_summary$truth
ci_summary$abs_bias <- abs(ci_summary$bias)
ci_summary$K <- K
ci_summary$B_boot <- B_boot
ci_summary$alpha <- alpha

ci_summary <- ci_summary[order(ci_summary$theta, ci_summary$n, ci_summary$p), ]

mean_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}
median_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}
q25_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else
    unname(quantile(x, 0.25, na.rm = TRUE))
}
q75_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else
    unname(quantile(x, 0.75, na.rm = TRUE))
}

selector_keys <- unique(selector_results[
  c("theta", "n", "m", "L", "lambda", "p_min", "h_max", "n_min")
])

selector_summary_rows <- vector("list", nrow(selector_keys))
for (i in seq_len(nrow(selector_keys))) {
  dat <- selector_results[
    selector_results$theta == selector_keys$theta[i] &
      selector_results$n == selector_keys$n[i],
  ]
  selected_dat <- dat[dat$selected, ]

  selector_summary_rows[[i]] <- data.frame(
    theta = selector_keys$theta[i],
    n = selector_keys$n[i],
    m = selector_keys$m[i],
    L = selector_keys$L[i],
    lambda = selector_keys$lambda[i],
    p_min = selector_keys$p_min[i],
    h_max = selector_keys$h_max[i],
    n_min = selector_keys$n_min[i],
    selection_rate = mean(dat$selected),
    mean_p_selected = mean_na(selected_dat$p_stab),
    median_p_selected = median_na(selected_dat$p_stab),
    q25_p_selected = q25_na(selected_dat$p_stab),
    q75_p_selected = q75_na(selected_dat$p_stab),
    mean_abs_error = mean_na(selected_dat$selected_abs_error),
    selected_ci_coverage = mean_na(as.numeric(selected_dat$selected_cover)),
    mean_half_width = mean_na(selected_dat$half_width_stab),
    mean_n_lower = mean_na(selected_dat$n_lower_stab)
  )
}

selector_summary <- do.call(rbind, selector_summary_rows)
selector_summary <- selector_summary[
  order(selector_summary$theta, selector_summary$n),
]

## -------------------------------------------------------------------------
## Figures
## -------------------------------------------------------------------------

plot_case <- function(dat, selector_dat, theta, n, figure_dir) {
  dat <- dat[order(dat$p), ]
  mean_coverage <- mean(dat$coverage)
  m_value <- unique(dat$m)
  p_stab_median <- selector_dat$median_p_selected[1]

  y_range <- range(dat$truth, dat$mean_rho_hat, dat$mean_lower, dat$mean_upper,
                   finite = TRUE)
  pad <- 0.06 * diff(y_range)
  if (!is.finite(pad) || pad == 0) pad <- 0.05
  y_range <- y_range + c(-pad, pad)

  file_theta <- gsub("-", "neg", as.character(theta), fixed = TRUE)
  file_theta <- gsub("\\.", "p", file_theta)
  fig_file <- file.path(
    figure_dir,
    sprintf("fgm_theta_%s_n%d.pdf", file_theta, n)
  )

  pdf(fig_file, width = 8.5, height = 6.2, pointsize = 12, bg = "white")
  op <- par(mar = c(5.2, 5.8, 4.2, 1.4), mgp = c(3.3, 1.0, 0),
            cex.lab = 1.7, cex.axis = 1.4, cex.main = 1.6,
            bg = "white")
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  plot(dat$p, dat$mean_rho_hat,
       type = "l", lwd = 5, col = "#1f77b4",
       xlab = "p",
       ylab = expression(hat(rho)[m*","*n](p)),
       ylim = y_range,
       main = bquote("FGM copula: " * theta == .(theta) *
                       ", n = " * .(n) * ", m = " * .(m_value)))
  lines(dat$p, dat$truth, lwd = 5, col = "#222222", lty = 2)
  lines(dat$p, dat$mean_lower, lwd = 5, col = "#d62728", lty = 3)
  lines(dat$p, dat$mean_upper, lwd = 5, col = "#d62728", lty = 3)

  if (is.finite(p_stab_median)) {
    abline(v = p_stab_median, lwd = 5, col = "#2ca02c", lty = 4)
  }

  legend("bottomleft",
         legend = c("Bernstein estimator",
                    "True function",
                    "Confidence intervals",
                    "Selected threshold",
                    sprintf("Mean coverage = %.1f%%", 100 * mean_coverage)),
         col = c("#1f77b4", "#222222", "#d62728", "#2ca02c", NA),
         lty = c(1, 2, 3, 4, NA),
         lwd = c(5, 5, 5, 5, NA),
         bty = "n",
         cex = 1.35,
         y.intersp = 1.15)
}

case_keys <- unique(ci_summary[c("theta", "n")])
for (i in seq_len(nrow(case_keys))) {
  theta_i <- case_keys$theta[i]
  n_i <- case_keys$n[i]
  dat_i <- ci_summary[ci_summary$theta == theta_i & ci_summary$n == n_i, ]
  selector_i <- selector_summary[
    selector_summary$theta == theta_i & selector_summary$n == n_i,
  ]
  plot_case(dat_i, selector_i, theta_i, n_i, figure_dir)
}

