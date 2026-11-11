## Real-data example: effect of tamoxifen on 5-year RMST (RA / IPW / DR)
## Data: `gse` from the drRML package (Choi et al., 2022), i.e. the GSE6532
## breast cancer cohort (Loi et al., 2007; NCBI GEO). It is loaded from the
## package and is not copied into this repo.
## remotes::install_github("taehwa015/drRML")

library(survival)
library(pseudo)
library(drRML)

set.seed(123)

## 1. Data
data(gse)
str(gse)

tau  <- 365 * 5                 # 5 years (time is in days)
covs <- names(gse)[3:6]         # Age, Grade, Size, Er
dat  <- data.frame(y = gse$Y, delta = gse$Censor, a = gse$trt, gse[, 3:6])

stopifnot(all(dat$a %in% c(0, 1)))
stopifnot(tau < max(dat$y))

ps_formula  <- as.formula(paste("a ~", paste(covs, collapse = " + ")))
reg_formula <- as.formula(paste("pm ~ a * (", paste(covs, collapse = " + "), ")"))

## 2. Estimators
est_rmst <- function(d, tau) {
  d$pm <- pseudomean(d$y, d$delta, tau)

  ps <- predict(glm(ps_formula, family = binomial(), data = d), type = "response")

  reg_fit <- lm(reg_formula, data = d)
  q1 <- predict(reg_fit, newdata = transform(d, a = 1))
  q0 <- predict(reg_fit, newdata = transform(d, a = 0))

  c(RA  = mean(q1 - q0),
    IPW = mean(d$a * d$pm / ps - (1 - d$a) * d$pm / (1 - ps)),
    DR  = mean((q1 - q0) +
                 d$a * (d$pm - q1) / ps -
                 (1 - d$a) * (d$pm - q0) / (1 - ps)))
}

## 3. Estimates and bootstrap SE
est <- est_rmst(dat, tau)

ps_hat <- predict(glm(ps_formula, family = binomial(), data = dat), type = "response")
cat("Propensity score range:", round(range(ps_hat), 3), "\n")
cat("Treated:", sum(dat$a == 1), " Control:", sum(dat$a == 0),
    " Event rate:", round(mean(dat$delta), 3), "\n\n")

B <- 200
boot <- t(replicate(B, est_rmst(dat[sample(nrow(dat), replace = TRUE), ], tau)))

se <- apply(boot, 2, sd)
result <- data.frame(
  Method   = names(est),
  Estimate = unname(est),
  SE       = unname(se),
  Lower    = unname(est - 1.96 * se),
  Upper    = unname(est + 1.96 * se)
)
result[, -1] <- round(result[, -1], 2)
cat("Difference in 5-year RMST (treated - control), days; bootstrap B =", B, "\n")
print(result, row.names = FALSE)

## 4. Compare with drRML::dr_rml()
fit_pkg <- dr_rml(Y = dat$y, Censor = dat$delta, A = dat$a,
                  Xps = gse[, 3:6], Xreg = gse[, 3:6],
                  L = tau, PS = "logit", Reg = "lm", nboot = 10)
cat("\ndr_rml() difference in 5-year RMST (days):", round(fit_pkg$ace, 2),
    " (SE with nboot = 10:", round(fit_pkg$se, 2), ")\n")
