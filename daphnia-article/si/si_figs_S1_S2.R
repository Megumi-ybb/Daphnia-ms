# =============================================================================
# si_figs_S1_S2.R
# Stand-alone, editable code for the two SI figures:
#   * Figure S-1  (label fig:mcap-coverage-plots) -- MCAP coverage panels
#                  rn / sigP / sigF / ri  (four PDFs, arranged 2x2 in si.Rnw)
#   * Figure S-2  (label fig:lrt-sim)              -- LRT vs chi-square toy model
#
# The S-1 code is extracted from Mixed-species/SIRJPF2/coverage_study/collect_coverage.R
# (the plotting part only). The S-2 code is the verbatim toy-model code that lives
# in the si.Rnw knitr chunks <<lrt-test>> and <<fig-lrt-sim>>.
#
# Run from daphnia-article/si/ :  Rscript si_figs_S1_S2.R
#   S-1 -> writes coverage/coverage_plot_{rn,sigP,sigF,ri}.pdf  (where si.Rnw reads them)
#   S-2 -> writes fig_S2_lrt_sim.pdf
# Edit freely; nothing here is sourced by si.Rnw, so changes only affect the PDFs
# you regenerate.
# =============================================================================


# =============================================================================
# FIGURE S-1 : MCAP coverage panels
# =============================================================================
# Needs pomp (for mcap()) + ggplot2, and the profile_<param>_<b>.rds produced by
# the coverage sweep. EDIT coverage_dir to point at the coverage_study folder.
suppressMessages({ library(pomp); library(ggplot2) })

coverage_dir <- "../../Mixed-species/SIRJPF2/coverage_study"  # EDIT if needed
out_dir      <- "coverage"                                    # si.Rnw reads coverage/*.pdf
level        <- 0.95
params_S1    <- c("rn", "sigP", "sigF", "ri")                 # the four panels

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
true_params <- readRDS(file.path(coverage_dir, "simulated_data", "true_params.rds"))

suppressMessages(library(patchwork))

# Read one parameter's profiles -> MCAP CIs. Coverage is scored from the RAW mcap
# interval; for plotting, open (non-finite) bounds are clamped to the scanned
# profile grid so flat-profile / degenerate cases (notably ri) still draw a
# full-height segment instead of silently dropping the line.
coverage_data <- function(param, coverage_dir, true_params, level = 0.95) {
  true_log_val <- log(true_params[param]); log_col <- paste0("log_", param)
  files <- list.files(file.path(coverage_dir, "coverage_results"),
                      pattern = sprintf("^profile_%s_\\d+\\.rds$", param), full.names = TRUE)
  grid_lo <- Inf; grid_hi <- -Inf; rows <- list()
  for (f in files) {
    d <- readRDS(f)
    grid_lo <- min(grid_lo, d[[log_col]], na.rm = TRUE)
    grid_hi <- max(grid_hi, d[[log_col]], na.rm = TRUE)
    # mcap() emits "NaNs produced" from sqrt(se_stat_squared)/sqrt(se_total_squared)
    # whenever its weighted quadratic fit has curvature a <= 0 (flat / non-quadratic /
    # locally-inverted profile near the max). Those se_* fields are diagnostics ONLY:
    # the CI used here comes from delta = qchisq(.,1)*(a*se_mc_squared + 0.5) and
    # range(grid[logLik_diff < delta]), which never reference se_stat/se_total. The
    # warning is therefore harmless to coverage (most frequent for ri's flat profiles),
    # so mute just that NaN-from-sqrt warning; real mcap errors still propagate to NULL.
    m <- tryCatch(
      withCallingHandlers(
        mcap(d$loglik, d[[log_col]], level = level, span = 0.95, Ngrid = 1000),
        warning = function(w) {
          if (grepl("NaNs produced", conditionMessage(w))) invokeRestart("muffleWarning")
        }),
      error = function(e) NULL)
    if (is.null(m)) next
    rows[[length(rows) + 1L]] <- data.frame(
      ci_lo = m$ci[1], ci_hi = m$ci[2], mle = m$mle,
      covers = (true_log_val >= m$ci[1]) & (true_log_val <= m$ci[2]))
  }
  res <- do.call(rbind, rows)
  res$open  <- !is.finite(res$ci_lo) | !is.finite(res$ci_hi)   # open intervals (flat profile)
  res$ci_lo <- pmax(ifelse(is.finite(res$ci_lo), res$ci_lo, grid_lo), grid_lo)
  res$ci_hi <- pmin(ifelse(is.finite(res$ci_hi), res$ci_hi, grid_hi), grid_hi)
  res$status    <- factor(ifelse(res$covers, "covers", "misses"), levels = c("covers", "misses"))
  res$b_ordered <- rank(res$mle, ties.method = "first")
  attr(res, "true_log") <- true_log_val
  res
}

# One compact panel. show_x = FALSE drops the x-axis so the two panels in a column
# share it (only the bottom row carries the axis).
coverage_panel <- function(res, name, tag, show_x = TRUE, base = 8.5) {
  tl <- attr(res, "true_log"); n <- nrow(res); ncov <- sum(res$status == "covers")
  g <- ggplot(res, aes(x = b_ordered)) +
    geom_segment(aes(xend = b_ordered, y = ci_lo, yend = ci_hi, color = status),
                 linewidth = 0.35) +
    geom_point(aes(y = mle), size = 0.45, color = "grey20") +
    geom_hline(yintercept = tl, linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_color_manual(values = c(covers = "steelblue", misses = "red"),
                       drop = FALSE, name = NULL) +
    labs(x = "dataset (ordered by MLE)", y = paste0("log(", name, ")")) +
    ggtitle(sprintf("(%s) %s : %d/%d = %.0f%%", tag, name, ncov, n, 100 * ncov / n)) +
    theme_bw(base_size = base) +
    theme(plot.title  = element_text(size = base, hjust = 0, margin = margin(b = 1)),
          plot.margin = margin(2, 4, 2, 3),
          legend.position = "bottom")
  if (!show_x)
    g <- g + theme(axis.title.x = element_blank(),
                   axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  g
}

cat("Figure S-1 (MCAP coverage, combined):\n")
d_rn   <- coverage_data("rn",   coverage_dir, true_params, level)
d_sigP <- coverage_data("sigP", coverage_dir, true_params, level)
d_sigF <- coverage_data("sigF", coverage_dir, true_params, level)
d_ri   <- coverage_data("ri",   coverage_dir, true_params, level)

# layout:  (a) rn   (b) sigP   <- top row,    x-axis hidden (shared down each column)
#          (c) sigF (d) ri     <- bottom row, x-axis shown
fig_S1 <- ( coverage_panel(d_rn,   "rn",   "a", show_x = FALSE) |
            coverage_panel(d_sigP, "sigP", "b", show_x = FALSE) ) /
          ( coverage_panel(d_sigF, "sigF", "c", show_x = TRUE)  |
            coverage_panel(d_ri,   "ri",   "d", show_x = TRUE)  ) +
          plot_layout(guides = "collect") & theme(legend.position = "bottom")

out_S1 <- file.path(out_dir, "coverage_combined.pdf")
ggsave(out_S1, fig_S1, width = 6.5, height = 4.6)
cat("  saved -> ", out_S1, "\n")


# =============================================================================
# FIGURE S-2 : LRT vs chi-square toy model
# -----------------------------------------------------------------------------
# The function definitions below are verbatim from the si.Rnw <<lrt-test>> chunk;
# the plotting call at the very bottom is the <<fig-lrt-sim>> chunk, wrapped in a
# pdf() device. Edit U / N / sigma / tau / nsim to change the figure.
# =============================================================================

############################################################
## Random effects versus fixed effects LRT simulation
## Base R only
############################################################

############################################################
## Model:
##
## H0: Y_ij = mu + R_i + epsilon_ij
##     R_i ~ N(0, tau^2), epsilon_ij ~ N(0, sigma^2)
##
## H1: Y_ij = mu_i + R_i + epsilon_ij
##
## Balanced panel: i = 1,...,I, j = 1,...,J
##
## We use the marginal likelihood integrating out R_i.
##
## The text uses U in place of I and V in place of J, for
## consistency with PanelPOMP notation
############################################################


############################################################
## Sufficient statistics
############################################################

suff_stats <- function(y) {
  y <- as.matrix(y)
  I <- nrow(y)
  J <- ncol(y)
  
  if (I < 2) stop("Need I >= 2.")
  if (J < 2) stop("Need J >= 2 for H1 to have finite residual variation.")
  
  ybar_i <- rowMeans(y)
  ybar <- mean(ybar_i)
  
  SSW <- sum(sweep(y, 1, ybar_i, "-")^2)
  SSB <- sum((ybar_i - ybar)^2)
  
  list(
    I = I,
    J = J,
    ybar_i = ybar_i,
    ybar = ybar,
    SSW = SSW,
    SSB = SSB
  )
}


############################################################
## Marginal likelihood details
##
## For one unit, covariance is
##
##   Sigma = sigma^2 I_J + tau^2 11'
##
## Define
##
##   A = sigma^2
##   B = sigma^2 + J tau^2
##
## Then
##
##   |Sigma| = A^(J-1) B
##
## and the quadratic form decomposes into
##
##   SSW_i / A + J (ybar_i - mean)^2 / B.
##
## Under H0, mean = mu.
## Under H1, mean = mu_i, so the second term is zero
## after setting mu_i = ybar_i.
############################################################


############################################################
## H0 fit by closed-form marginal MLE
############################################################

fit_h0_stats <- function(I, J, SSW, SSB) {
  if (I < 2) stop("Need I >= 2.")
  if (J < 2) stop("Need J >= 2.")
  if (SSW <= 0) stop("Degenerate SSW <= 0.")
  
  ## Unconstrained MLEs in terms of
  ## A = sigma^2 and B = sigma^2 + J tau^2.
  A_unc <- SSW / (I * (J - 1))
  B_unc <- J * SSB / I
  
  ## Constraint is B >= A, equivalent to tau^2 >= 0.
  if (B_unc >= A_unc) {
    A <- A_unc
    B <- B_unc
    boundary <- FALSE
  } else {
    ## Boundary tau^2 = 0, so B = A.
    A <- (SSW + J * SSB) / (I * J)
    B <- A
    boundary <- TRUE
  }
  
  sigma2 <- A
  tau2 <- max((B - A) / J, 0)
  
  loglik <- -0.5 * (
    I * J * log(2 * pi) +
      I * (J - 1) * log(A) +
      I * log(B) +
      SSW / A +
      J * SSB / B
  )
  
  list(
    loglik = loglik,
    sigma2 = sigma2,
    tau2 = tau2,
    A = A,
    B = B,
    boundary_tau0 = boundary
  )
}


fit_h0 <- function(y) {
  s <- suff_stats(y)
  fit <- fit_h0_stats(
    I = s$I,
    J = s$J,
    SSW = s$SSW,
    SSB = s$SSB
  )
  fit$mu <- s$ybar
  fit$I <- s$I
  fit$J <- s$J
  fit$SSW <- s$SSW
  fit$SSB <- s$SSB
  fit
}


############################################################
## H1 fit by closed-form marginal MLE
##
## Under H1:
##
##   mu_i_hat = ybar_i
##
## and the likelihood is maximized at tau^2_hat = 0.
##
## Therefore:
##
##   sigma2_hat = SSW / (I J)
##
## This is the usual ML denominator, not the unbiased
## residual denominator I(J - 1).
############################################################

fit_h1_stats <- function(I, J, SSW) {
  if (I < 2) stop("Need I >= 2.")
  if (J < 2) stop("Need J >= 2.")
  if (SSW <= 0) {
    stop("Degenerate SSW <= 0. H1 likelihood is unbounded.")
  }
  
  sigma2 <- SSW / (I * J)
  tau2 <- 0
  
  loglik <- -0.5 * (
    I * J * log(2 * pi) +
      I * J * log(sigma2) +
      SSW / sigma2
  )
  
  list(
    loglik = loglik,
    sigma2 = sigma2,
    tau2 = tau2,
    boundary_tau0 = TRUE
  )
}


fit_h1 <- function(y) {
  s <- suff_stats(y)
  fit <- fit_h1_stats(
    I = s$I,
    J = s$J,
    SSW = s$SSW
  )
  fit$mu_i <- s$ybar_i
  fit$I <- s$I
  fit$J <- s$J
  fit$SSW <- s$SSW
  fit$SSB <- s$SSB
  fit
}


############################################################
## LRT statistic
############################################################

lrt_stats <- function(I, J, SSW, SSB) {
  h0 <- fit_h0_stats(I = I, J = J, SSW = SSW, SSB = SSB)
  h1 <- fit_h1_stats(I = I, J = J, SSW = SSW)
  
  LRT <- 2 * (h1$loglik - h0$loglik)
  
  ## Protect against tiny negative numerical roundoff.
  LRT <- max(LRT, 0)
  
  list(
    LRT = LRT,
    h0 = h0,
    h1 = h1
  )
}


lrt <- function(y) {
  s <- suff_stats(y)
  lrt_stats(
    I = s$I,
    J = s$J,
    SSW = s$SSW,
    SSB = s$SSB
  )
}


############################################################
## Simulate a full panel under H0
############################################################

simulate_panel_h0 <- function(I, J, mu = 0, sigma = 1, tau = 1) {
  R <- rnorm(I, mean = 0, sd = tau)
  eps <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)
  
  y <- matrix(mu + R, nrow = I, ncol = J) + eps
  y
}


############################################################
## Faster simulation using exact null distribution of
## sufficient statistics
##
## Under H0:
##
##   SSW / sigma^2 ~ chi-square_{I(J-1)}
##
## and
##
##   SSB / (tau^2 + sigma^2 / J) ~ chi-square_{I-1}
##
## independently.
############################################################

simulate_suff_h0 <- function(I, J, sigma = 1, tau = 1) {
  SSW <- sigma^2 * rchisq(1, df = I * (J - 1))
  SSB <- (tau^2 + sigma^2 / J) * rchisq(1, df = I - 1)
  
  list(
    I = I,
    J = J,
    SSW = SSW,
    SSB = SSB
  )
}


############################################################
## Main simulation routine
############################################################

simulate_lrt <- function(nsim = 10000,
                         I = 10,
                         J = 20,
                         mu = 0,
                         sigma = 1,
                         tau = 1,
                         seed = NULL,
                         method = c("suff", "panel")) {
  
  method <- match.arg(method)
  
  if (!is.null(seed)) set.seed(seed)
  
  out <- data.frame(
    sim = seq_len(nsim),
    LRT = numeric(nsim),
    p_chisq = numeric(nsim),
    h0_sigma2 = numeric(nsim),
    h0_tau2 = numeric(nsim),
    h0_boundary_tau0 = logical(nsim),
    h1_sigma2 = numeric(nsim),
    h1_tau2 = numeric(nsim)
  )
  
  df_chisq <- I - 1
  
  for (b in seq_len(nsim)) {
    
    if (method == "suff") {
      s <- simulate_suff_h0(I = I, J = J, sigma = sigma, tau = tau)
      fit <- lrt_stats(I = I, J = J, SSW = s$SSW, SSB = s$SSB)
    } else {
      y <- simulate_panel_h0(I = I, J = J, mu = mu, sigma = sigma, tau = tau)
      fit <- lrt(y)
    }
    
    out$LRT[b] <- fit$LRT
    out$p_chisq[b] <- pchisq(fit$LRT, df = df_chisq, lower.tail = FALSE)
    
    out$h0_sigma2[b] <- fit$h0$sigma2
    out$h0_tau2[b] <- fit$h0$tau2
    out$h0_boundary_tau0[b] <- fit$h0$boundary_tau0
    
    out$h1_sigma2[b] <- fit$h1$sigma2
    out$h1_tau2[b] <- fit$h1$tau2
  }
  
  attr(out, "I") <- I
  attr(out, "J") <- J
  attr(out, "sigma") <- sigma
  attr(out, "tau") <- tau
  attr(out, "df_chisq") <- df_chisq
  attr(out, "method") <- method
  
  out
}


############################################################
## Summaries
############################################################

summarize_lrt <- function(sim,
                          alpha = c(0.10, 0.05, 0.01),
                          probs = c(0.50, 0.90, 0.95, 0.99)) {
  
  df_chisq <- attr(sim, "df_chisq")
  
  empirical_quantiles <- quantile(sim$LRT, probs = probs, names = FALSE)
  chisq_quantiles <- qchisq(probs, df = df_chisq)
  
  qtab <- data.frame(
    prob = probs,
    empirical = empirical_quantiles,
    chisq = chisq_quantiles,
    ratio_empirical_to_chisq = empirical_quantiles / chisq_quantiles
  )
  
  rtab <- data.frame(
    alpha = alpha,
    nominal_chisq_cutoff = qchisq(1 - alpha, df = df_chisq),
    empirical_rejection_rate = sapply(alpha, function(a) mean(sim$p_chisq <= a))
  )
  
  list(
    settings = list(
      I = attr(sim, "I"),
      J = attr(sim, "J"),
      sigma = attr(sim, "sigma"),
      tau = attr(sim, "tau"),
      df_chisq = df_chisq,
      method = attr(sim, "method"),
      nsim = nrow(sim)
    ),
    mean_LRT = mean(sim$LRT),
    median_LRT = median(sim$LRT),
    sd_LRT = sd(sim$LRT),
    h0_boundary_tau0_frequency = mean(sim$h0_boundary_tau0),
    h1_tau0_frequency = mean(sim$h1_tau2 == 0),
    quantiles = qtab,
    rejection_rates = rtab
  )
}


############################################################
## Plots
############################################################

plot_lrt_hist <- function(sim, breaks = 60, main = NULL,...) {
  df_chisq <- attr(sim, "df_chisq")
  
  if (is.null(main)) {
    main <- paste0(
      "LRT under H0: I=", attr(sim, "I"),
      ", J=", attr(sim, "J"),
      ", tau=", attr(sim, "tau")
    )
  }
  
  hist(
    sim$LRT,
    breaks = breaks,
    freq = FALSE,
    col = "gray85",
    border = "white",
    main = main,
    xlab = "LRT statistic",
    ...
  )
  
  curve(
    dchisq(x, df = df_chisq),
    add = TRUE,
    col = "red",
    lwd = 2
  )
  
  legend(
    "topright",
    legend = c("Simulation", paste0("Chi-square df=", df_chisq)),
    col = c("gray50", "red"),
    lwd = c(8, 2),
    bty = "n"
  )
}


plot_lrt_qq <- function(sim, main = NULL) {
  df_chisq <- attr(sim, "df_chisq")
  
  if (is.null(main)) {
    main <- paste0(
      "QQ plot against chi-square df=", df_chisq
    )
  }
  
  n <- nrow(sim)
  theo <- qchisq(ppoints(n), df = df_chisq)
  emp <- sort(sim$LRT)
  
  plot(
    theo,
    emp,
    pch = 19,
    cex = 0.5,
    xlab = paste0("Chi-square(", df_chisq, ") quantiles"),
    ylab = "Empirical LRT quantiles",
    main = main
  )
  abline(0, 1, col = "red", lwd = 2)
}


############################################################
## Grid experiment
############################################################

run_grid <- function(nsim = 5000,
                     I = 10,
                     J_values = c(3, 5, 10, 20, 50),
                     tau_values = c(0, 0.25, 0.5, 1, 2),
                     sigma = 1,
                     seed = 1) {
  
  set.seed(seed)
  
  grid <- expand.grid(
    J = J_values,
    tau = tau_values
  )
  
  ans <- vector("list", nrow(grid))
  
  for (k in seq_len(nrow(grid))) {
    J <- grid$J[k]
    tau <- grid$tau[k]
    
    sim <- simulate_lrt(
      nsim = nsim,
      I = I,
      J = J,
      sigma = sigma,
      tau = tau,
      method = "suff",
      seed = NULL
    )
    
    df_chisq <- I - 1
    
    ans[[k]] <- data.frame(
      I = I,
      J = J,
      sigma = sigma,
      tau = tau,
      nsim = nsim,
      mean_LRT = mean(sim$LRT),
      median_LRT = median(sim$LRT),
      q90_LRT = as.numeric(quantile(sim$LRT, 0.90)),
      q95_LRT = as.numeric(quantile(sim$LRT, 0.95)),
      q99_LRT = as.numeric(quantile(sim$LRT, 0.99)),
      chisq_q90 = qchisq(0.90, df = df_chisq),
      chisq_q95 = qchisq(0.95, df = df_chisq),
      chisq_q99 = qchisq(0.99, df = df_chisq),
      reject_10 = mean(sim$p_chisq <= 0.10),
      reject_05 = mean(sim$p_chisq <= 0.05),
      reject_01 = mean(sim$p_chisq <= 0.01),
      h0_boundary_tau0_frequency = mean(sim$h0_boundary_tau0)
    )
  }
  
  do.call(rbind, ans)
}


############################################################
## Parametric bootstrap p-value for an observed panel
############################################################

bootstrap_lrt_pvalue <- function(y,
                                 nsim = 5000,
                                 seed = NULL,
                                 method = c("suff", "panel")) {
  
  method <- match.arg(method)
  
  if (!is.null(seed)) set.seed(seed)
  
  s <- suff_stats(y)
  fit0 <- fit_h0(y)
  obs <- lrt(y)$LRT
  
  sim <- simulate_lrt(
    nsim = nsim,
    I = s$I,
    J = s$J,
    mu = fit0$mu,
    sigma = sqrt(fit0$sigma2),
    tau = sqrt(fit0$tau2),
    method = method,
    seed = NULL
  )
  
  p_boot <- (1 + sum(sim$LRT >= obs)) / (nsim + 1)
  p_chisq <- pchisq(obs, df = s$I - 1, lower.tail = FALSE)
  
  list(
    observed_LRT = obs,
    p_bootstrap = p_boot,
    p_chisq = p_chisq,
    fit_h0 = fit0,
    sim = sim
  )
}


# ---- Figure S-2: produce the plot (was the <<fig-lrt-sim>> chunk) ----
# Compact size for the SI (was 7 x 4); tight margins + smaller text.
pdf("fig_S2_lrt_sim.pdf", width = 4.6, height = 2.9)
par(mar = c(3.6, 3.6, 0.6, 0.6), mgp = c(2.1, 0.6, 0), cex = 0.7, cex.lab = 0.95)

U = 10
N = 20
sigma = 1
tau = 0.5

sim <- simulate_lrt(
  nsim = 20000,
  I = U,
  J = N,
  sigma = sigma,
  tau = tau,
  seed = 123,
  method = "suff"
)

summ <- summarize_lrt(sim)
# print(summ)

plot_lrt_hist(sim,main="",ylim = c(0,0.1),xlim=c(0,42))
dev.off()
cat("Figure S-2 saved -> fig_S2_lrt_sim.pdf\n")
