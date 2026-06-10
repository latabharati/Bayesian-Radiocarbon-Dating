rm(list = ls())
#install.packages("rcarbon")
#install.packages("coda")
library(rcarbon)
library(coda)

# Question 3: Metropolis-within-Gibbs Sampler for Bayesian Radiocarbon Dating
# Application: Ordered archaeological samples from one phase
set.seed(372)

# Load the IntCal20 data
intcal20 <- read.csv(
  "https://intcal.org/curves/intcal20.14c",
  header = FALSE,
  comment.char = "#",
  col.names = c("CALBP", "C14BP", "C14ERR", "Delta14C", "Delta14CERR")
)

# Check key columns are numeric and print basic info
stopifnot(is.numeric(intcal20$CALBP), is.numeric(intcal20$C14BP))
cat("Calibration curve loaded:", nrow(intcal20), "data points\n")
cat("Calendar age range:", range(intcal20$CALBP), "cal BP\n\n")

# Create interpolation function to evaluate f(theta) at any calendar age
cal_curve <- approxfun(
  x = intcal20$CALBP,
  y = intcal20$C14BP,
  rule = 2
)

# Test the function
cat("f(3000 cal BP) =", cal_curve(3000), "C14 BP\n\n")

# Setting up the data
n_samples <- 5
c14_ages <- c(2980, 3050, 3080, 3120, 3200)
c14_sigma <- c(35, 40, 35, 30, 45)

# Sort by radiocarbon age to match expected order
order_idx <- order(c14_ages)
c14_ages <- c14_ages[order_idx]
c14_sigma <- c14_sigma[order_idx]

theta_min <- 2800
theta_max <- 3600

cat("Data after sorting:\n")
print(data.frame(Sample = 1:n_samples, C14_Age = c14_ages, Sigma = c14_sigma))

# Log-likelihood function
log_likelihood <- function(theta_val, x_obs, sigma) {
  mu <- cal_curve(theta_val)
  -0.5 * ((x_obs - mu) / sigma)^2
}

# Metropolis-within-Gibbs sampler for ordered calendar ages
gibbs_radiocarbon <- function(c14_ages,
                              c14_sigma,
                              theta_min,
                              theta_max,
                              n_iter = 60000,
                              burnin = 10000,
                              proposal_sd = 50,
                              init_theta = NULL) {
  n <- length(c14_ages)
  n_store <- n_iter - burnin
  
  # Starting values
  if (is.null(init_theta)) {
    theta <- numeric(n)
    
    for (i in 1:n) {
      best_idx <- which.min(abs(intcal20$C14BP - c14_ages[i]))
      theta[i] <- intcal20$CALBP[best_idx]
    }
    
    theta <- sort(theta)
  } else {
    theta <- sort(init_theta)
  }
  
  cat("Initial theta values:", round(theta, 1), "\n\n")
  
  # Store samples after burn-in
  samples <- matrix(NA, nrow = n_store, ncol = n)
  colnames(samples) <- paste0("theta_", 1:n)
  
  accept <- rep(0, n)
  
  # Main MCMC loop
  for (iter in 1:n_iter) {
    for (i in 1:n) {
      lower <- if (i == 1) theta_min else theta[i - 1]
      upper <- if (i == n) theta_max else theta[i + 1]
      
      theta_star <- rnorm(1, mean = theta[i], sd = proposal_sd)
      
      # Reject if ordering constraint is broken
      if (theta_star <= lower || theta_star >= upper) {
        next
      }
      
      log_alpha <- log_likelihood(theta_star, c14_ages[i], c14_sigma[i]) -
        log_likelihood(theta[i], c14_ages[i], c14_sigma[i])
      
      if (log(runif(1)) < log_alpha) {
        theta[i] <- theta_star
        
        if (iter > burnin) {
          accept[i] <- accept[i] + 1
        }
      }
    }
    
    if (iter > burnin) {
      samples[iter - burnin, ] <- theta
    }
    
    if (iter %% 10000 == 0) {
      cat("Iteration", iter, "/", n_iter, "- current theta:", round(theta, 1), "\n")
    }
  }
  
  cat("\nPost-burn-in acceptance rates:\n")
  acc_rates <- round(accept / n_store, 3)
  names(acc_rates) <- paste0("theta_", 1:n)
  print(acc_rates)
  cat("(Aim: 0.20 - 0.50)\n\n")
  
  return(samples)
}

# Run the sampler
samples <- gibbs_radiocarbon(
  c14_ages = c14_ages,
  c14_sigma = c14_sigma,
  theta_min = theta_min,
  theta_max = theta_max,
  n_iter = 60000,
  burnin = 10000,
  proposal_sd = 50
)

cat("========== Posterior Summary ==========\n")

results <- data.frame(
  Sample = paste0("Sample ", 1:n_samples),
  C14_Age_BP = c14_ages,
  Sigma = c14_sigma,
  Post_Mean = round(colMeans(samples), 1),
  Post_SD = round(apply(samples, 2, sd), 1),
  CI_2.5 = round(apply(samples, 2, quantile, 0.025), 1),
  CI_97.5 = round(apply(samples, 2, quantile, 0.975), 1),
  CI_Width = round(
    apply(samples, 2, quantile, 0.975) -
      apply(samples, 2, quantile, 0.025),
    1
  )
)

print(results)

# Monte Carlo Error
mcmc_obj <- as.mcmc(samples)
ess <- effectiveSize(mcmc_obj)
mc_error <- apply(samples, 2, sd) / sqrt(ess)

cat("========== Monte Carlo Error Assessment ==========\n")

mc_table <- data.frame(
  Sample = paste0("Sample ", 1:n_samples),
  N_Iterations = nrow(samples),
  ESS = round(ess, 1),
  ESS_pct = paste0(round(100 * ess / nrow(samples), 1), "%"),
  MC_Error_yrs = round(mc_error, 3)
)

print(mc_table)
cat("\n")

# Naive MC error assumes independent samples, so it is too small
naive_mc_error <- apply(samples, 2, sd) / sqrt(nrow(samples))

cat("Naive MC error (ignoring autocorrelation - incorrect):\n")
print(round(naive_mc_error, 3))

# Figure 1: IntCal20 Calibration Curve with Posterior Results
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

cal_idx <- intcal20$CALBP >= theta_min & intcal20$CALBP <= theta_max

plot(
  intcal20$CALBP[cal_idx],
  intcal20$C14BP[cal_idx],
  type = "l",
  lwd = 2,
  col = "black",
  main = "Figure 1: IntCal20 Calibration Curve with Posterior Means and 95% CIs",
  xlab = "Calendar Age (cal BP)",
  ylab = "Radiocarbon Age (C14 BP)",
  cex.main = 1.1,
  cex.lab = 1.0
)

lines(
  intcal20$CALBP[cal_idx],
  intcal20$C14BP[cal_idx] + 2 * intcal20$C14ERR[cal_idx],
  col = "grey70",
  lty = 2
)

lines(
  intcal20$CALBP[cal_idx],
  intcal20$C14BP[cal_idx] - 2 * intcal20$C14ERR[cal_idx],
  col = "grey70",
  lty = 2
)

points(colMeans(samples), c14_ages, pch = 19, col = "red", cex = 1.2)

label_offset_x <- c(-60, -30, 60, 80, 60)
label_offset_y <- c(20, -28, 20, -28, 20)

for (i in 1:n_samples) {
  arrows(
    x0 = quantile(samples[, i], 0.025),
    y0 = c14_ages[i],
    x1 = quantile(samples[, i], 0.975),
    y1 = c14_ages[i],
    angle = 90,
    code = 3,
    length = 0.06,
    col = "red",
    lwd = 1.5
  )
  
  text(
    x = colMeans(samples)[i] + label_offset_x[i],
    y = c14_ages[i] + label_offset_y[i],
    labels = paste0("S", i),
    col = "red",
    cex = 0.8
  )
}

legend(
  "topleft",
  legend = c("IntCal20 curve", "Curve ±2 sigma", "Posterior mean + 95% CI"),
  col = c("black", "grey70", "red"),
  lty = c(1, 2, NA),
  pch = c(NA, NA, 19),
  bty = "n",
  cex = 0.9
)

# Figure 2: Trace Plots
par(mfrow = c(1, 1))

for (i in 1:n_samples) {
  y <- samples[1:5000, i]
  
  plot(
    y,
    type = "l",
    col = "darkblue",
    ylim = c(min(y) - 10, max(y) + 10),
    main = bquote("Figure 2." * .(i) * ": Trace Plot for " * theta[.(i)]),
    xlab = "Iteration (post burn-in)",
    ylab = "Calendar Age (cal BP)"
  )
  
  abline(h = mean(samples[, i]), col = "red", lty = 2, lwd = 2)
  
  legend(
    "topright",
    legend = paste0("Mean = ", round(mean(samples[, i]), 1), " cal BP"),
    col = "red",
    lty = 2,
    bty = "n",
    cex = 0.9
  )
}

# Figure 3: Running Mean Plots
par(mfrow = c(1, 1))

for (i in 1:n_samples) {
  running_mean <- cumsum(samples[, i]) / seq_along(samples[, i])
  true_mean <- mean(samples[, i])
  
  plot(
    running_mean,
    type = "l",
    col = "darkblue",
    ylim = range(running_mean) + c(-10, 10),
    main = bquote("Figure 3." * .(i) * ": Running Mean for " * theta[.(i)]),
    xlab = "Iteration (post burn-in)",
    ylab = "Running Mean (cal BP)",
    cex.main = 1.0
  )
  
  abline(h = true_mean, col = "red", lty = 2, lwd = 2)
  
  legend(
    "topright",
    legend = paste0("Final mean = ", round(true_mean, 1), " cal BP"),
    col = "red",
    lty = 2,
    bty = "n",
    cex = 0.9
  )
}

# Figure 4: Posterior Histograms
par(mfrow = c(1, 1))

for (i in 1:n_samples) {
  post_mean <- mean(samples[, i])
  ci_lo <- quantile(samples[, i], 0.025)
  ci_hi <- quantile(samples[, i], 0.975)
  
  hist(
    samples[, i],
    breaks = 50,
    probability = TRUE,
    col = "lightsteelblue",
    border = "white",
    main = bquote("Figure 4." * .(i) * ": Posterior of " * theta[.(i)]),
    xlab = "Calendar Age (cal BP)",
    ylab = "Density"
  )
  
  lines(density(samples[, i]), col = "darkgreen", lwd = 2)
  abline(v = post_mean, col = "red", lwd = 2)
  abline(v = ci_lo, col = "blue", lty = 2, lwd = 1.5)
  abline(v = ci_hi, col = "blue", lty = 2, lwd = 1.5)
  
  legend(
    "topright",
    legend = c(
      paste0("Mean: ", round(post_mean, 1)),
      paste0("95% CI: [", round(ci_lo, 1), ", ", round(ci_hi, 1), "]")
    ),
    col = c("red", "blue"),
    lty = c(1, 2),
    bty = "n",
    cex = 0.9
  )
}

# Figure 5: Autocorrelation Plots
par(mfrow = c(1, 1))

for (i in 1:n_samples) {
  acf(
    samples[, i],
    lag.max = 60,
    main = bquote("Figure 5." * .(i) * ": Autocorrelation Function for " * theta[.(i)]),
    xlab = "Lag",
    col = "darkblue",
    lwd = 2
  )
}

# Figure 6: Summary Plot
par(mfrow = c(1, 1), mar = c(5, 7, 4, 2))

post_means <- colMeans(samples)
ci_lo_all <- apply(samples, 2, quantile, 0.025)
ci_hi_all <- apply(samples, 2, quantile, 0.975)

plot(
  post_means,
  1:n_samples,
  xlim = range(c(ci_lo_all, ci_hi_all)) + c(-30, 30),
  ylim = c(0.5, n_samples + 0.5),
  pch = 19,
  col = "red",
  cex = 1.3,
  main = "Figure 6: Posterior Means and 95% Credible Intervals for All Calendar Ages",
  xlab = "Calendar Age (cal BP)",
  ylab = "",
  yaxt = "n",
  cex.main = 1.0
)

axis(
  2,
  at = 1:n_samples,
  labels = paste0("Sample ", 1:n_samples, "  (C14 = ", c14_ages, " BP)"),
  las = 1,
  cex.axis = 0.85
)

for (i in 1:n_samples) {
  lines(c(ci_lo_all[i], ci_hi_all[i]), c(i, i), col = "red", lwd = 2)
  lines(c(ci_lo_all[i], ci_lo_all[i]), c(i - 0.1, i + 0.1), col = "red", lwd = 2)
  lines(c(ci_hi_all[i], ci_hi_all[i]), c(i - 0.1, i + 0.1), col = "red", lwd = 2)
}

abline(v = post_means, col = "grey85", lty = 3)

legend(
  "bottomright",
  legend = c("Posterior mean", "95% Credible Interval"),
  col = "red",
  pch = c(19, NA),
  lty = c(NA, 1),
  lwd = c(NA, 2),
  bty = "n",
  cex = 0.9
)

# Final Results Summary
cat("Posterior calendar ages (cal BP):\n\n")

for (i in 1:n_samples) {
  cat(sprintf(
    "  Sample %d: Mean = %.1f cal BP | 95%% CI = [%.1f, %.1f] | ESS = %.0f | MC Error = %.3f yrs\n",
    i,
    post_means[i],
    ci_lo_all[i],
    ci_hi_all[i],
    ess[i],
    mc_error[i]
  ))
}

cat(
  "\nOrdering constraint theta_1 < ... < theta_5 satisfied:",
  all(diff(post_means) > 0),
  "\n"
)

print(round(naive_mc_error, 3))

# Multiple-chain diagnostic using Gelman-Rubin R-hat
init_list <- list(
  c(2850, 3050, 3200, 3350, 3500),
  c(2900, 3100, 3250, 3400, 3550),
  c(3000, 3150, 3300, 3420, 3580),
  c(2950, 3120, 3280, 3450, 3590)
)

multi_chains <- vector("list", 4)

for (ch in 1:4) {
  set.seed(372 + ch)
  
  multi_chains[[ch]] <- gibbs_radiocarbon(
    c14_ages = c14_ages,
    c14_sigma = c14_sigma,
    theta_min = theta_min,
    theta_max = theta_max,
    n_iter = 60000,
    burnin = 10000,
    proposal_sd = 50,
    init_theta = init_list[[ch]]
  )
}

mcmc_list <- mcmc.list(lapply(multi_chains, as.mcmc))
rhat_result <- gelman.diag(mcmc_list, autoburnin = FALSE)

cat("========== Gelman-Rubin R-hat Diagnostic ==========\n")
print(rhat_result)

cat("R-hat values:\n")
print(round(rhat_result$psrf[, 1], 3))