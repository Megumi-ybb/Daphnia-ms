#!/usr/bin/env Rscript
# ============================================================================
# lrt_diagnostic_table.R -- the compact per-alternative SI diagnostic that
# rebuts "it's just the optimizer / just the filter" in one table, plus a
# Monte-Carlo (binomial) CI on each bootstrap p-value (B=100 adequacy check).
# Columns: df | Lambda_obs | boot null mean | var/2df | alt_se/null_se |
#          neg-fraction | p_chisq | p_boot [95% MC CI] | implied rho* |
#          sigP & sigF collapse (alt/null)  -- the boundary-mechanism fingerprint.
# Run from this directory.  Writes lrt_diagnostic_table.csv.
# ============================================================================
B <- 100
I <- 8; J <- 10
alts <- c("xi","theta_Si_theta_Sn","theta_P","theta_Ii_theta_In","ri_rn","probi_probn","f_Si_f_Sn")
rmap <- list(xi=5, theta_Si_theta_Sn=6, theta_P=4, theta_Ii_theta_In=2, ri_rn=8, probi_probn=3, f_Si_f_Sn=7)

## --- exact-MLE toy (for the implied random-effect strength rho*) ---
lambda_from_SS <- function(SSW, SSB, I, J) {
  N <- I*J; ll1 <- -N/2*(log(2*pi)+log(SSW/N)+1)
  s2 <- SSW/(I*(J-1)); tau2 <- SSB/I; interior <- tau2 >= s2
  ll0 <- numeric(length(SSW))
  ll0[interior]  <- -N/2*log(2*pi) - I/2*log(tau2[interior]) - I*(J-1)/2*log(s2[interior]) - N/2
  sp <- (SSW+SSB)/N; ll0[!interior] <- -N/2*(log(2*pi)+log(sp[!interior])+1)
  2*(ll1-ll0)
}
sim_block <- function(Bn, phi2) {
  SS <- vapply(seq_len(Bn), function(.) {
    Y <- matrix(rnorm(I*J), I, J) + rnorm(I, 0, sqrt(phi2)); um <- rowMeans(Y)
    c(sum((Y-um)^2), J*sum((um-mean(Y))^2)) }, numeric(2))
  lambda_from_SS(SS[1,], SS[2,], I, J)
}
set.seed(11)
toy_mean <- function(rho, nb) { phi2 <- (rho-1)/J; m <- 0; for (k in 1:nb) m <- m + sim_block(4000, phi2); mean(m) }
implied_rho <- function(target, nb) {           # rho >= 1 (phi^2 >= 0); NA if deflated below the toy floor
  if (target <= toy_mean(1.001, nb)) return(NA_real_)
  tryCatch(uniroot(function(r) toy_mean(r, nb) - target, c(1.001, 80))$root, error = function(e) NA_real_)
}

## --- observed log-liks ---
load("../../../data/Target_dynamics/para/Target_para_loglik_df.rds")  # target_para_parameter_table
Tg <- target_para_parameter_table
nll <- Tg["all_shared","ll"]; kN <- (Tg["all_shared","AIC"] + 2*nll)/2

## --- null bootstrap fits (loglik, se, sigP, sigF) ---
nf <- lapply(1:B, function(b){f<-sprintf("results_null/lrt_null_%d.rds",b); if(file.exists(f)) readRDS(f)})
nll_v <- sapply(nf, function(x) if(is.null(x)) NA else x$ll)
nse_v <- sapply(nf, function(x) if(is.null(x)) NA else x$se)
nP    <- sapply(nf, function(x) if(is.null(x)) NA else unname(x$coef["sigP"]))
nF    <- sapply(nf, function(x) if(is.null(x)) NA else unname(x$coef["sigF"]))

rows <- lapply(alts, function(a) {
  ri <- rmap[[a]]; alo <- Tg[ri,"block_ll"]; kA <- (Tg[ri,"block_AIC"] + 2*alo)/2
  df <- round(kA - kN); Lobs <- 2*(alo - nll)
  af <- lapply(1:B, function(b){f<-sprintf("results_alt/lrt_%s_%d.rds",a,b); if(file.exists(f)) readRDS(f)})
  al_v <- sapply(af, function(x) if(is.null(x)) NA else x$ll)
  ase_v <- sapply(af, function(x) if(is.null(x)) NA else x$se)
  aP <- sapply(af, function(x) if(is.null(x)) NA else unname(x$coef["sigP"]))
  aF <- sapply(af, function(x) if(is.null(x)) NA else unname(x$coef["sigF"]))
  v <- !is.na(nll_v) & !is.na(al_v); L <- 2*(al_v[v] - nll_v[v]); Bv <- sum(v)
  k <- sum(L >= Lobs); pboot <- (1+k)/(1+Bv)
  ci <- binom.test(k, Bv)$conf.int        # Clopper-Pearson on exceedance prob = MC CI for p_boot
  data.frame(alt=a, df=df, Lambda_obs=round(Lobs,2),
             boot_mean=round(mean(L),2), var_2df=round(var(L)/(2*df),2),
             altse_nullse=round(mean(ase_v[v])/mean(nse_v[v]),2),
             negfrac=round(mean(L<0),3),
             p_chisq=signif(pchisq(Lobs,df,lower.tail=FALSE),3),
             p_boot=round(pboot,3),
             p_boot_CI=sprintf("[%.3f,%.3f]", ci[1], ci[2]),
             rho_star=round(implied_rho(mean(L), df/7),2),
             sigP_alt_null=round(mean(aP[v])/mean(nP[v]),2),
             sigF_alt_null=round(mean(aF[v])/mean(nF[v]),2),
             stringsAsFactors=FALSE)
})
tab <- do.call(rbind, rows)
cat("\n=== SIRJPF2 LRT diagnostic table ===\n"); print(tab, row.names=FALSE)
write.csv(tab, "lrt_diagnostic_table.csv", row.names=FALSE)

cat("\n=== correlation: process-noise collapse vs inflation ===\n")
infl <- tab$boot_mean/tab$df; coll <- pmin(tab$sigP_alt_null, tab$sigF_alt_null, na.rm=TRUE)
cat(sprintf("Pearson(inflation, min sig-collapse) = %.2f ; Spearman = %.2f\n",
            cor(infl, coll), cor(infl, coll, method="spearman")))

cat("\n=== toy calibration controls (exact MLE, no filter, no optimizer) ===\n")
for (rho in c(1, 1.93)) for (nb in 1:2) {
  df <- 7*nb; phi2 <- (rho-1)/J
  L <- sim_block(15000, phi2); if (nb==2) L <- L + sim_block(15000, phi2)
  cat(sprintf("  rho=%.2f df=%2d : mean=%.2f  true size of nominal-5%% chi^2 = %.3f\n",
              rho, df, mean(L), mean(L > qchisq(0.95, df))))
}
cat("\nSaved lrt_diagnostic_table.csv\n")
