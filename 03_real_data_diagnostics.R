## Checks for the real-data example: overlap, balance, sensitivity
## Data: `gse` from the drRML package (GSE6532, Loi et al., 2007), loaded with
## data(gse). Point estimates only, no bootstrap.

library(survival)
library(pseudo)
library(drRML)

set.seed(123)
data(gse)

tau  <- 365 * 5
covs <- c("Age", "Size", "Grade", "Er")
dat  <- data.frame(y = gse$Y, delta = gse$Censor, a = gse$trt, gse[, 3:6])
dat$pm <- pseudomean(dat$y, dat$delta, tau)   # computed once on the full sample

base_ps  <- a  ~ Age + Size + Grade + Er
base_reg <- pm ~ a * (Age + Size + Grade + Er)

## 1. Data
cat("---- Data ----\n")
cat("\nEr by treatment:\n");    print(table(Er = dat$Er, trt = dat$a))
cat("\nGrade by treatment:\n"); print(table(Grade = dat$Grade, trt = dat$a))
cat("\nMean of covariates by treatment:\n")
print(aggregate(dat[, covs], by = list(trt = dat$a), FUN = mean))
cat("\nPseudo-values:\n"); print(summary(dat$pm))
fit_base <- lm(base_reg, data = dat)
cat("\nOutcome model rank:", fit_base$rank, "of", length(coef(fit_base)), "\n")

## 2. Propensity score and weights
cat("\n---- Propensity score and weights ----\n")
ps_fit <- glm(base_ps, family = binomial(), data = dat)
dat$ps <- predict(ps_fit, type = "response")
cat("\nPropensity model:\n")
print(round(summary(ps_fit)$coefficients, 3))

cat("\nPS by treatment group:\n")
print(tapply(dat$ps, dat$a, summary))

w <- ifelse(dat$a == 1, 1 / dat$ps, 1 / (1 - dat$ps))
ess <- function(x) sum(x)^2 / sum(x^2)    # effective sample size
cat("\nLargest 10 IPW weights:", round(sort(w, decreasing = TRUE)[1:10], 1), "\n")
cat("ESS, treated:  ", round(ess(w[dat$a == 1]), 1), "of", sum(dat$a == 1), "\n")
cat("ESS, untreated:", round(ess(w[dat$a == 0]), 1), "of", sum(dat$a == 0), "\n")
cat("PS < 0.05 or > 0.95:", sum(dat$ps < 0.05 | dat$ps > 0.95), "of", nrow(dat), "\n")

par(mfrow = c(1, 2))
hist(dat$ps[dat$a == 1], breaks = 20, xlim = c(0, 1),
     main = "Treated", xlab = "Estimated propensity score")
hist(dat$ps[dat$a == 0], breaks = 20, xlim = c(0, 1),
     main = "Untreated", xlab = "Estimated propensity score")
par(mfrow = c(1, 1))
dev.copy(png, "ps_overlap.png", width = 900, height = 400)
dev.off()

## 3. Covariate balance (standardized mean difference, treated - untreated)
cat("\n---- Covariate balance (SMD) ----\n")
smd <- function(x, a, w = rep(1, length(x))) {
  m1 <- weighted.mean(x[a == 1], w[a == 1])
  m0 <- weighted.mean(x[a == 0], w[a == 0])
  v1 <- sum(w[a == 1] * (x[a == 1] - m1)^2) / sum(w[a == 1])
  v0 <- sum(w[a == 0] * (x[a == 0] - m0)^2) / sum(w[a == 0])
  (m1 - m0) / sqrt((v1 + v0) / 2)
}
bal <- data.frame(
  Covariate  = covs,
  Unweighted = sapply(covs, function(v) smd(dat[[v]], dat$a)),
  IPW        = sapply(covs, function(v) smd(dat[[v]], dat$a, w))
)
bal[, -1] <- round(bal[, -1], 3)
print(bal, row.names = FALSE)

## 4. Sensitivity (point estimates, days)
## clip: cap the PS at [0.05, 0.95]; trim: drop patients outside that range
## (trimming changes the target population)
cat("\n---- Sensitivity ----\n")
est_fun <- function(d, ps_formula, reg_formula, clip = NULL, trim = NULL) {
  ps <- predict(glm(ps_formula, family = binomial(), data = d), type = "response")
  if (!is.null(clip)) ps <- pmin(pmax(ps, clip[1]), clip[2])
  keep <- rep(TRUE, nrow(d))
  if (!is.null(trim)) keep <- ps >= trim[1] & ps <= trim[2]
  d <- d[keep, ]; ps <- ps[keep]
  fit <- lm(reg_formula, data = d)
  q1 <- predict(fit, newdata = transform(d, a = 1))
  q0 <- predict(fit, newdata = transform(d, a = 0))
  c(RA  = mean(q1 - q0),
    IPW = mean(d$a * d$pm / ps - (1 - d$a) * d$pm / (1 - ps)),
    DR  = mean((q1 - q0) + d$a * (d$pm - q1) / ps -
                 (1 - d$a) * (d$pm - q0) / (1 - ps)),
    n   = nrow(d))
}

rich_ps  <- a  ~ Age + Size + Grade + Er + I(Age^2) + I(Size^2)
rich_reg <- pm ~ a * (Age + Size + Grade + Er + I(Age^2) + I(Size^2))
fac_reg  <- pm ~ a * (Age + Size + factor(Grade) + Er)

sens <- rbind(
  "Baseline"                              = est_fun(dat, base_ps, base_reg),
  "PS clipped to [0.05, 0.95]"            = est_fun(dat, base_ps, base_reg, clip = c(0.05, 0.95)),
  "PS trimmed (drop outside 0.05-0.95)"   = est_fun(dat, base_ps, base_reg, trim = c(0.05, 0.95)),
  "Quadratic Age and Size (PS + outcome)" = est_fun(dat, rich_ps, rich_reg),
  "Grade as factor (outcome model)"       = est_fun(dat, base_ps, fac_reg)
)
print(round(sens, 1))
