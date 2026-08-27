# ==============================================================
# inference.R -- from the fits to the aqueous pKa and its standard
# error, plus the lack-of-fit diagnostics.
#
# Both routes end in the same quadratic form (constants.R):
#     pKa_aqueous = g' theta ,   Var = g' V g .
#
# Route A: V is the analytic covariance of fit.R.
# Route B: the three per-mixture estimates have variances nu_k, and
#     Var(theta~) = (Z'Z)^-1 Z' [diag(nu) + tau2 I] Z (Z'Z)^-1 .
#   tau2_hat = max(0, SS_res - sum (1-h_k) nu_k) is the lack-of-fit
#   component (a DerSimonian-Laird moment estimator): the scatter of
#   the three points about the fitted line beyond what their own
#   measurement variances explain.  It is 0 whenever the points are
#   consistent with the line, and then SE_c reduces to the pure
#   measurement propagation SE_a = sqrt(sum ell_k^2 nu_k).
#
#   SE_a : measurement propagation only  (tau2 forced to 0)
#   SE_b : ordinary-least-squares mean-response SE, i.e. the
#          equal-variance convention -- reported for reference only
#   SE_c : the convention used in the report (measurement + misfit)
#
# Q = SS_res / sum (1-h_k) nu_k compares observed scatter with the
# scatter the measurement variances alone would produce; under exact
# adherence to the solvent model Q ~ chi^2_1.  tau2_hat > 0 iff Q > 1.
# ==============================================================

analyze <- function(data, spec, label = "") {
  # ---------------- Route A: joint fit ----------------
  fit  <- fit_joint(data, spec)
  th   <- coef_natural(fit, spec)          # reported coefficients
  V_j  <- vcov_joint(fit, data, spec)      # analytic covariance of theta
  sig  <- summary(fit)$sigma
  aq_j <- as.numeric(t(spec$g) %*% th)
  se_j <- aqueous_se(V_j, spec$g)
  # cross-check: the fit is done in centered coordinates, where the
  # aqueous value is the first parameter, so the optimizer reports its
  # standard error directly.
  se_j_num <- unname(summary(fit)$coefficients["aq", "Std. Error"])

  # ---------------- Route B: two-step ----------------
  ind   <- lapply(1:3, function(k) fit_mixture(data[data$grp == k, ]))
  pKa_k <- sapply(ind, function(f) unname(coef(f)["pKa"]))
  se_k  <- sapply(ind, function(f)
    unname(summary(f)$coefficients["pKa", "Std. Error"]))
  sig_k <- sapply(ind, function(f) summary(f)$sigma)
  nu_k  <- se_k^2
  y     <- pKa_k - spec$offset            # offset-removed, linear in x
  ols   <- lm(y ~ x, data = data.frame(y = y, x = spec$x))
  th_t  <- unname(coef(ols))
  aq_t  <- as.numeric(t(spec$g) %*% th_t)

  SS_res  <- sum(residuals(ols)^2)
  noise_E <- sum((1 - spec$h) * nu_k)     # E[SS_res] when tau2 = 0
  tau2    <- max(0, SS_res - noise_E)
  Q       <- SS_res / noise_E
  p_Q     <- pchisq(Q, df = 1, lower.tail = FALSE)

  V_a <- spec$ZtZi %*% t(spec$Z) %*% diag(nu_k)        %*% spec$Z %*% spec$ZtZi
  V_c <- spec$ZtZi %*% t(spec$Z) %*% diag(nu_k + tau2) %*% spec$Z %*% spec$ZtZi
  se_a <- aqueous_se(V_a, spec$g)
  se_c <- aqueous_se(V_c, spec$g)
  se_b <- unname(predict(ols, newdata = data.frame(x = spec$x0),
                         se.fit = TRUE)$se.fit)

  R2_ols   <- summary(ols)$r.squared
  R2_joint <- 1 - sum((y - (th[1] + th[2] * spec$x))^2) / sum((y - mean(y))^2)
  n_k <- sapply(1:3, function(k) sum(data$grp == k))

  list(
    label = label, spec = spec, data = data, ols = ols,
    pKa_k = pKa_k, se_k = se_k, y = y,
    theta_joint = th, V_joint = V_j,
    theta_two = th_t, V_two_c = V_c,
    aq_joint = aq_j, se_joint = se_j, aq_two = aq_t, se_c = se_c,
    tau2 = tau2, Q = Q, p_Q = p_Q, R2_joint = R2_joint, R2_ols = R2_ols,
    row = data.frame(
      dataset = label, model = spec$model,
      n1 = n_k[1], n2 = n_k[2], n3 = n_k[3],
      pKa1 = pKa_k[1], pKa2 = pKa_k[2], pKa3 = pKa_k[3],
      SE1 = se_k[1], SE2 = se_k[2], SE3 = se_k[3],
      sig1 = sig_k[1], sig2 = sig_k[2], sig3 = sig_k[3],
      theta1 = th[1], theta2 = th[2],
      SE_theta1 = sqrt(V_j[1, 1]), SE_theta2 = sqrt(V_j[2, 2]),
      corr_theta = V_j[1, 2] / sqrt(V_j[1, 1] * V_j[2, 2]),
      sigma = sig,
      theta1_two = th_t[1], theta2_two = th_t[2],
      SE_theta1_two = sqrt(V_c[1, 1]), SE_theta2_two = sqrt(V_c[2, 2]),
      corr_theta_two = V_c[1, 2] / sqrt(V_c[1, 1] * V_c[2, 2]),
      aq_joint = aq_j, SE_joint = se_j, SE_joint_numeric = se_j_num,
      aq_twostep = aq_t, SE_a = se_a, SE_b = se_b, SE_c = se_c,
      SS_res = SS_res, E0_SS_res = noise_E, tau2_hat = tau2,
      Q = Q, p_Q = p_Q, R2_joint = R2_joint, R2_ols = R2_ols,
      stringsAsFactors = FALSE))
}
