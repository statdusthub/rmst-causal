## RMST-based causal inference with pseudo-values: RA / IPW / DR
## Simulation: exponential survival times, independent censoring
##
## Estimand: Delta_tau = E[min{T(1), tau}] - E[min{T(0), tau}]

library(survival)
library(pseudo)

set.seed(123)

## 1. True RMST for an exponential distribution with the given rate
rmst_exp <- function(rate, tau) {
  (1 - exp(-rate * tau)) / rate
}

## 2. One simulation run
## assignment = "randomized": A ~ Bernoulli(0.5)
## assignment = "confounded": logit P(A = 1 | X) = alpha[1] + alpha[2]*x1 + alpha[3]*x2
one_run <- function(n = 300, tau = 10,
                    assignment = c("randomized", "confounded"),
                    alpha = c(0, 0.5, 0.5)) {
  assignment <- match.arg(assignment)

  ## ---- covariates & treatment ----
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  a <- if (assignment == "randomized") {
    rbinom(n, 1, 0.5)
  } else {
    rbinom(n, 1, plogis(alpha[1] + alpha[2] * x1 + alpha[3] * x2))
  }

  ## ---- exponential hazards ----
  rate  <- exp(-(x1 + x2 + a * (x1 + x2)))  # observed arm
  rate0 <- exp(-(x1 + x2))                  # control potential outcome
  rate1 <- exp(-2 * (x1 + x2))              # treated potential outcome

  ## ---- survival & independent censoring ----
  time_event <- rexp(n, rate)
  time_cens  <- rexp(n, exp(-3))
  delta <- as.numeric(time_event < time_cens)
  y <- pmin(time_event, time_cens)

  ## ---- true ATE of RMST (closed form, averaged over sampled covariates) ----
  true_ate <- mean(rmst_exp(rate1, tau) - rmst_exp(rate0, tau))

  ## ---- RMST pseudo-values ----
  pmean <- pseudomean(y, delta, tau)

  ## ---- propensity score model ----
  ps_fit <- glm(a ~ x1 + x2, family = binomial())
  ps <- predict(ps_fit, type = "response")

  ## ---- outcome regression on pseudo-values ----
  reg_fit <- lm(pmean ~ a * (x1 + x2))
  q0 <- predict(reg_fit, newdata = data.frame(a = 0, x1, x2))
  q1 <- predict(reg_fit, newdata = data.frame(a = 1, x1, x2))

  ## ---- estimators ----
  est_ipw <- mean(a * pmean / ps - (1 - a) * pmean / (1 - ps))   # IPW
  est_reg <- mean(q1 - q0)                                       # Regression adjustment
  est_dr  <- mean((q1 - q0) +                                    # Doubly robust (AIPW)
                    a * (pmean - q1) / ps -
                    (1 - a) * (pmean - q0) / (1 - ps))

  c(true = true_ate, reg = est_reg, ipw = est_ipw, dr = est_dr)
}

## 3. Repeat the simulation
B <- 100
res <- as.data.frame(t(replicate(B, one_run(assignment = "randomized"))))

## 4. Summary statistics
est <- c("reg", "ipw", "dr")
summary_table <- data.frame(
  Method = toupper(est),
  Mean = unname(sapply(est, function(m) mean(res[[m]]))),
  SD   = unname(sapply(est, function(m) sd(res[[m]]))),
  Bias = unname(sapply(est, function(m) mean(res[[m]] - res$true))),
  MC_SE_of_Bias = unname(sapply(est, function(m) sd(res[[m]] - res$true) / sqrt(B)))
)

## 5. Output
cat("True ATE (average across simulations):", mean(res$true), "\n\n")
cat("Summary table (B =", B, "):\n")
summary_table[, -1] <- round(summary_table[, -1], 4)
print(summary_table, row.names = FALSE)
