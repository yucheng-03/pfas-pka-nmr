# ==============================================================
# fit.R -- the two estimation routes, written once for both models.
#
# Route A (joint / aggregated): one non-linear least squares fit of
#   all titration points, with pKa_k = theta' d_k + offset_k imposed
#   inside the Henderson-Hasselbalch mean function.  8 parameters:
#   theta plus the three pairs (delta_A_k, delta_HA_k).
#
# Route B (two-step): one 3-parameter fit per mixture, then ordinary
#   least squares of the three pKa estimates on d_k.
#
# The covariance of theta from Route A is computed ANALYTICALLY as
#   V = sigma^2 M^{-1},   M = sum_k w_k d_k d_k' ,
#   w_k = || (I - P_k) u_k ||^2 ,
# where u_k holds the sensitivities of the fitted curve to pKa_k and
# P_k projects onto the two limiting-shift directions.  The
# optimizer's finite-difference vcov() is NOT used: on these
# ill-conditioned fits (corr(theta_1, theta_2) < -0.98) it was found
# to be unreliable.  vcov() is still reported by analyze() as a
# numerical cross-check.
# ==============================================================
suppressPackageStartupMessages(library(minpack.lm))

# Henderson-Hasselbalch mean shift
hh <- function(pH, dA, dHA, pKa) (dA * 10^pH + dHA * 10^pKa) / (10^pH + 10^pKa)

# ---- Route B, stage 1: one mixture, three free parameters ----
fit_mixture <- function(sub)
  nlsLM(ChemShift ~ hh(pH, dA, dHA, pKa), data = sub,
        start = list(dA  = unname(sub$ChemShift[which.max(sub$pH)]),
                     dHA = unname(sub$ChemShift[which.min(sub$pH)]),
                     pKa = unname(mean(sub$pH))),
        control = nls.lm.control(maxiter = 500))

# ---- Route A: all mixtures at once, solvent model imposed ----
# The fit is carried out in the CENTERED covariate (x_k - x0), i.e. in
# the coordinates (aq, slope) with
#     pKa_k = aq + slope * (x_k - x0) + offset_k ,
# so that the aqueous value is a parameter of the fit.  This is the
# same model, only better conditioned: in the natural coordinates the
# two coefficients are correlated at about -0.99 and the optimizer
# stops at slightly different points along that ridge depending on how
# the covariate is written (A/eps vs A*(1/eps)).  Centering removes the
# ridge; across the datasets used here it attains the lowest residual
# sum of squares of the parameterizations tried.  coef_natural() maps
# the result back to the reported coefficients.
fit_joint <- function(data, spec) {
  data$xc  <- spec$x[data$grp] - spec$x0
  data$off <- spec$offset[data$grp]
  mean_fun <- function(pH, grp, xc, off, aq, slope,
                       dA_1, dHA_1, dA_2, dHA_2, dA_3, dHA_3) {
    dA  <- ifelse(grp == 1, dA_1,  ifelse(grp == 2, dA_2,  dA_3))
    dHA <- ifelse(grp == 1, dHA_1, ifelse(grp == 2, dHA_2, dHA_3))
    hh(pH, dA, dHA, aq + slope * xc + off)
  }
  # starting values: per-mixture fits, then OLS on the offset-removed values
  pre  <- sapply(1:3, function(k)
    unname(coef(fit_mixture(data[data$grp == k, ]))["pKa"]))
  lin  <- lm(y ~ xc, data = data.frame(y = pre - spec$offset,
                                       xc = spec$x - spec$x0))
  sdA  <- sapply(1:3, function(k) { s <- data[data$grp == k, ]
    unname(s$ChemShift[which.max(s$pH)]) })
  sdHA <- sapply(1:3, function(k) { s <- data[data$grp == k, ]
    unname(s$ChemShift[which.min(s$pH)]) })
  nlsLM(ChemShift ~ mean_fun(pH, grp, xc, off, aq, slope,
                             dA_1, dHA_1, dA_2, dHA_2, dA_3, dHA_3),
        data = data,
        start = list(aq = unname(coef(lin)[1]), slope = unname(coef(lin)[2]),
                     dA_1 = unname(sdA[1]), dHA_1 = unname(sdHA[1]),
                     dA_2 = unname(sdA[2]), dHA_2 = unname(sdHA[2]),
                     dA_3 = unname(sdA[3]), dHA_3 = unname(sdHA[3])),
        control = nls.lm.control(maxiter = 1024, ftol = 1e-14, ptol = 1e-14))
}

# (aq, slope) -> the reported coefficients theta = (theta_1, theta_2):
#   pKa_k = theta_1 + theta_2 x_k + offset_k , theta_1 = aq - slope * x0
coef_natural <- function(fit, spec) {
  cf <- coef(fit)
  c(unname(cf["aq"]) - unname(cf["slope"]) * spec$x0, unname(cf["slope"]))
}

# ---- effective information weights w_k of the joint fit ----
joint_weights <- function(fit, data, spec) {
  th  <- coef_natural(fit, spec)
  cf  <- coef(fit)
  pKa <- th[1] + th[2] * spec$x[data$grp] + spec$offset[data$grp]
  alp <- 1 / (1 + 10^(pKa - data$pH))          # degree of dissociation
  sapply(1:3, function(k) {
    i   <- data$grp == k
    Del <- unname(cf[paste0("dHA_", k)] - cf[paste0("dA_", k)])
    u   <- log(10) * Del * alp[i] * (1 - alp[i])   # sensitivity to pKa_k
    V   <- cbind(alp[i], 1 - alp[i])               # limiting-shift directions
    as.numeric(t(u) %*% u - t(u) %*% V %*% solve(t(V) %*% V) %*% t(V) %*% u)
  })
}

# ---- analytic covariance of theta from the joint fit ----
vcov_joint <- function(fit, data, spec) {
  w   <- joint_weights(fit, data, spec)
  sig <- summary(fit)$sigma
  W   <- sum(w); xbar <- sum(w * spec$x) / W
  S   <- sum(w * (spec$x - xbar)^2)
  sig^2 * matrix(c(S + W * xbar^2, -W * xbar, -W * xbar, W), 2, 2) / (W * S)
}
