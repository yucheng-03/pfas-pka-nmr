# Plain model code for the five-compound, 40/50/60% ACN analysis.
# This dataset is separate from the PFBA/PFPrA 60/70/80% analysis.
# Function order: get_constants(), fit_ys(), fit_concentration(),
#                 fit_ys_separated(), fit_concentration_separated().

# =========== Fixed solvent constants ===========
# Values follow R/constants.R in this repository,
# which cites the collaborator email dated 2026-07-15.
# water_activity is a_H2O, not the water mole fraction in the workbooks.
# The pure-water row supplies the aqueous read-off point for both models.
# Compute reciprocals and logarithms in the model equations to avoid rounding.

get_constants <- function() {
  constants <- data.frame(
    condition      = c("40% ACN", "50% ACN", "60% ACN", "pure water"),
    acn_fraction   = c(0.40,      0.50,      0.60,      0.00),
    epsilon        = c(63.50,     58.71,     53.91,     78.33),
    water_activity = c(0.913,     0.908,     0.901,     1.000)
  )

  return(constants)
}

# =========== YS model: joint fit of all three mixtures ===========
# Input: one compound, with columns pH, ChemShift, and Condition.
# pH must already contain the intended values (currently Adjusted pH).
# Condition must be "40% ACN", "50% ACN", or "60% ACN".
#
# YS equation: pKa = A / epsilon + B - log10(water_activity).
# In pure water, pka_water = A / epsilon_water + B.
# We fit pka_water and A, then recover B = pka_water - A / epsilon_water.
# This is the same equation in the coordinates used by the current analysis.

fit_ys <- function(data) {
  # 1. Match each observation to its solvent constants by condition name.
  constants <- get_constants()
  mixtures <- constants[constants$condition != "pure water", ]
  epsilon_water <- constants$epsilon[constants$condition == "pure water"]

  stopifnot(all(c("pH", "ChemShift", "Condition") %in% names(data)))
  if (any(!is.finite(data$pH)) || any(!is.finite(data$ChemShift))) {
    stop("pH and ChemShift must contain finite numbers; no rows are dropped.")
  }
  data$grp <- match(data$Condition, mixtures$condition)
  if (anyNA(data$grp)) stop("Condition must be 40% ACN, 50% ACN, or 60% ACN.")
  data$epsilon <- mixtures$epsilon[data$grp]
  data$water_activity <- mixtures$water_activity[data$grp]

  # 2. Write the full mean equation, with a separate plateau pair per mixture.
  ys_mean <- function(pH, grp, epsilon, water_activity, pka_water, A,
                      dA_40, dHA_40, dA_50, dHA_50, dA_60, dHA_60) {
    dA <- ifelse(grp == 1, dA_40, ifelse(grp == 2, dA_50, dA_60))
    dHA <- ifelse(grp == 1, dHA_40, ifelse(grp == 2, dHA_50, dHA_60))
    pKa <- pka_water + A * (1 / epsilon - 1 / epsilon_water) -
      log10(water_activity)
    return((dA * 10^pH + dHA * 10^pKa) / (10^pH + 10^pKa))
  }

  # 3. Separate curve fits provide starting values only.
  # All eight parameters will be estimated together in step 4.
  pka_start <- numeric(3)
  dA_start <- numeric(3)
  dHA_start <- numeric(3)
  for (k in 1:3) {
    group_data <- data[data$grp == k, ]
    if (nrow(group_data) < 4) stop("Each mixture needs at least four observations.")
    dA_start[k] <- group_data$ChemShift[which.max(group_data$pH)]
    dHA_start[k] <- group_data$ChemShift[which.min(group_data$pH)]
    curve_fit <- minpack.lm::nlsLM(
      ChemShift ~ (dA * 10^pH + dHA * 10^pKa) / (10^pH + 10^pKa),
      data = group_data,
      start = list(dA = unname(dA_start[k]), dHA = unname(dHA_start[k]),
                   pKa = unname(mean(group_data$pH))),
      control = minpack.lm::nls.lm.control(maxiter = 500)
    )
    if (!isTRUE(curve_fit$convInfo$isConv)) stop("A starting curve fit did not converge.")
    pka_start[k] <- unname(coef(curve_fit)["pKa"])
  }
  start_data <- data.frame(
    corrected_pka = pka_start + log10(mixtures$water_activity),
    inverse_epsilon_centered = 1 / mixtures$epsilon - 1 / epsilon_water
  )
  line_start <- lm(corrected_pka ~ inverse_epsilon_centered, data = start_data)

  # 4. One joint fit: pka_water, A, and six limiting shifts are all free.
  fit <- minpack.lm::nlsLM(
    ChemShift ~ ys_mean(pH, grp, epsilon, water_activity, pka_water, A,
                       dA_40, dHA_40, dA_50, dHA_50, dA_60, dHA_60),
    data = data,
    start = list(
      pka_water = unname(coef(line_start)[1]), A = unname(coef(line_start)[2]),
      dA_40 = unname(dA_start[1]), dHA_40 = unname(dHA_start[1]),
      dA_50 = unname(dA_start[2]), dHA_50 = unname(dHA_start[2]),
      dA_60 = unname(dA_start[3]), dHA_60 = unname(dHA_start[3])
    ),
    control = minpack.lm::nls.lm.control(maxiter = 1024, ftol = 1e-14, ptol = 1e-14)
  )
  if (!isTRUE(fit$convInfo$isConv)) stop("The joint YS fit did not converge.")
  estimates <- coef(fit)
  pka_water <- unname(estimates["pka_water"])
  A <- unname(estimates["A"])
  B <- pka_water - A / epsilon_water

  # 5. Analytic derivatives of the mean with respect to ALL eight parameters.
  # alpha is the fraction in the high-pH form; 1 - alpha is the other fraction.
  pKa <- pka_water + A * (1 / data$epsilon - 1 / epsilon_water) -
    log10(data$water_activity)
  alpha <- 1 / (1 + 10^(pKa - data$pH))
  dA <- unname(estimates[c("dA_40", "dA_50", "dA_60")])[data$grp]
  dHA <- unname(estimates[c("dHA_40", "dHA_50", "dHA_60")])[data$grp]
  derivative_pka <- log(10) * (dHA - dA) * alpha * (1 - alpha)

  jacobian <- cbind(
    pka_water = derivative_pka,
    A = derivative_pka * (1 / data$epsilon - 1 / epsilon_water),
    dA_40 = alpha * (data$grp == 1), dHA_40 = (1 - alpha) * (data$grp == 1),
    dA_50 = alpha * (data$grp == 2), dHA_50 = (1 - alpha) * (data$grp == 2),
    dA_60 = alpha * (data$grp == 3), dHA_60 = (1 - alpha) * (data$grp == 3)
  )
  sigma_squared <- sum(residuals(fit)^2) / (nrow(data) - 8)
  covariance <- sigma_squared * solve(t(jacobian) %*% jacobian)
  se_water <- sqrt(covariance["pka_water", "pka_water"])

  # covariance uses the named coordinates (pka_water, A, dA_40, dHA_40, ...).
  # The nls fit object retains the six plateau estimates and fitted shifts.
  return(list(A = A, B = B, pka_water = pka_water, se_water = se_water,
              covariance = covariance, fit = fit))
}

# =========== ACN concentration model: joint fit of all three mixtures ===========
# Input: one compound, with columns pH, ChemShift, and Condition.
# pH must already contain the intended values (currently Adjusted pH).
# Condition must be "40% ACN", "50% ACN", or "60% ACN".
#
# Concentration equation: pKa = beta0 + beta1 * acn_fraction.
# acn_fraction is 0.40, 0.50, or 0.60; pure water has acn_fraction = 0.
# Therefore beta0 is the aqueous pKa itself.

fit_concentration <- function(data) {
  # 1. Match each observation to its ACN volume fraction by condition name.
  constants <- get_constants()
  mixtures <- constants[constants$condition != "pure water", ]

  stopifnot(all(c("pH", "ChemShift", "Condition") %in% names(data)))
  if (any(!is.finite(data$pH)) || any(!is.finite(data$ChemShift))) {
    stop("pH and ChemShift must contain finite numbers; no rows are dropped.")
  }
  data$grp <- match(data$Condition, mixtures$condition)
  if (anyNA(data$grp)) stop("Condition must be 40% ACN, 50% ACN, or 60% ACN.")
  data$acn_fraction <- mixtures$acn_fraction[data$grp]

  # 2. Write the full mean equation, with a separate plateau pair per mixture.
  concentration_mean <- function(pH, grp, acn_fraction, beta0, beta1,
                                 dA_40, dHA_40, dA_50, dHA_50, dA_60, dHA_60) {
    dA <- ifelse(grp == 1, dA_40, ifelse(grp == 2, dA_50, dA_60))
    dHA <- ifelse(grp == 1, dHA_40, ifelse(grp == 2, dHA_50, dHA_60))
    pKa <- beta0 + beta1 * acn_fraction
    return((dA * 10^pH + dHA * 10^pKa) / (10^pH + 10^pKa))
  }

  # 3. Separate curve fits provide starting values only.
  # Regress their pKa estimates directly on ACN fraction to start beta0/beta1.
  pka_start <- numeric(3)
  dA_start <- numeric(3)
  dHA_start <- numeric(3)
  for (k in 1:3) {
    group_data <- data[data$grp == k, ]
    if (nrow(group_data) < 4) stop("Each mixture needs at least four observations.")
    dA_start[k] <- group_data$ChemShift[which.max(group_data$pH)]
    dHA_start[k] <- group_data$ChemShift[which.min(group_data$pH)]
    curve_fit <- minpack.lm::nlsLM(
      ChemShift ~ (dA * 10^pH + dHA * 10^pKa) / (10^pH + 10^pKa),
      data = group_data,
      start = list(dA = unname(dA_start[k]), dHA = unname(dHA_start[k]),
                   pKa = unname(mean(group_data$pH))),
      control = minpack.lm::nls.lm.control(maxiter = 500)
    )
    if (!isTRUE(curve_fit$convInfo$isConv)) stop("A starting curve fit did not converge.")
    pka_start[k] <- unname(coef(curve_fit)["pKa"])
  }
  start_data <- data.frame(pKa = pka_start, acn_fraction = mixtures$acn_fraction)
  line_start <- lm(pKa ~ acn_fraction, data = start_data)

  # 4. One joint fit: beta0, beta1, and six limiting shifts are all free.
  fit <- minpack.lm::nlsLM(
    ChemShift ~ concentration_mean(pH, grp, acn_fraction, beta0, beta1,
                                  dA_40, dHA_40, dA_50, dHA_50, dA_60, dHA_60),
    data = data,
    start = list(
      beta0 = unname(coef(line_start)[1]), beta1 = unname(coef(line_start)[2]),
      dA_40 = unname(dA_start[1]), dHA_40 = unname(dHA_start[1]),
      dA_50 = unname(dA_start[2]), dHA_50 = unname(dHA_start[2]),
      dA_60 = unname(dA_start[3]), dHA_60 = unname(dHA_start[3])
    ),
    control = minpack.lm::nls.lm.control(maxiter = 1024, ftol = 1e-14, ptol = 1e-14)
  )
  if (!isTRUE(fit$convInfo$isConv)) stop("The joint concentration fit did not converge.")
  estimates <- coef(fit)
  beta0 <- unname(estimates["beta0"])
  beta1 <- unname(estimates["beta1"])

  # 5. Analytic derivatives of the mean with respect to ALL eight parameters.
  # The global derivatives are d(mean)/d(pKa) times 1 and acn_fraction.
  pKa <- beta0 + beta1 * data$acn_fraction
  alpha <- 1 / (1 + 10^(pKa - data$pH))
  dA <- unname(estimates[c("dA_40", "dA_50", "dA_60")])[data$grp]
  dHA <- unname(estimates[c("dHA_40", "dHA_50", "dHA_60")])[data$grp]
  derivative_pka <- log(10) * (dHA - dA) * alpha * (1 - alpha)

  jacobian <- cbind(
    beta0 = derivative_pka,
    beta1 = derivative_pka * data$acn_fraction,
    dA_40 = alpha * (data$grp == 1), dHA_40 = (1 - alpha) * (data$grp == 1),
    dA_50 = alpha * (data$grp == 2), dHA_50 = (1 - alpha) * (data$grp == 2),
    dA_60 = alpha * (data$grp == 3), dHA_60 = (1 - alpha) * (data$grp == 3)
  )
  sigma_squared <- sum(residuals(fit)^2) / (nrow(data) - 8)
  covariance <- sigma_squared * solve(t(jacobian) %*% jacobian)
  se_water <- sqrt(covariance["beta0", "beta0"])

  # covariance uses the named coordinates (beta0, beta1, dA_40, dHA_40, ...).
  return(list(beta0 = beta0, beta1 = beta1, pka_water = beta0, se_water = se_water,
              covariance = covariance, fit = fit))
}

# =========== YS model: separate curve fits, then an OLS line ===========
# Input columns are pH, ChemShift, and Condition, as in fit_ys().
# Stage 1 estimates a separate pKa and its SE for each mixture.
# Stage 2 fits pKa + log10(water_activity) = B + A / epsilon by ordinary LS.
# The current report uses SE_c (measurement error + lack-of-fit) for the
# aqueous value. SE_a and SE_b are also returned, explicitly as references.

fit_ys_separated <- function(data) {
  # 1. Identify the three mixtures and the pure-water read-off point.
  constants <- get_constants()
  mixtures <- constants[constants$condition != "pure water", ]
  epsilon_water <- constants$epsilon[constants$condition == "pure water"]
  n_mixtures <- nrow(mixtures)

  stopifnot(all(c("pH", "ChemShift", "Condition") %in% names(data)))
  if (any(!is.finite(data$pH)) || any(!is.finite(data$ChemShift))) {
    stop("pH and ChemShift must contain finite numbers; no rows are dropped.")
  }
  data$grp <- match(data$Condition, mixtures$condition)
  if (anyNA(data$grp)) stop("Condition must be 40% ACN, 50% ACN, or 60% ACN.")

  # 2. Fit all THREE parameters within each mixture, including both plateaus.
  # These fits are the final separated curves, not merely starting values.
  mixture_results <- mixtures
  mixture_results$n <- integer(n_mixtures)
  mixture_results$pka <- numeric(n_mixtures)
  mixture_results$se_pka <- numeric(n_mixtures)
  mixture_fits <- vector("list", n_mixtures)
  names(mixture_fits) <- mixtures$condition

  for (k in 1:n_mixtures) {
    group_data <- data[data$grp == k, ]
    if (nrow(group_data) < 4) stop("Each mixture needs at least four observations.")
    dA_start <- group_data$ChemShift[which.max(group_data$pH)]
    dHA_start <- group_data$ChemShift[which.min(group_data$pH)]
    curve_fit <- minpack.lm::nlsLM(
      ChemShift ~ (dA * 10^pH + dHA * 10^pKa) / (10^pH + 10^pKa),
      data = group_data,
      start = list(dA = unname(dA_start), dHA = unname(dHA_start),
                   pKa = unname(mean(group_data$pH))),
      control = minpack.lm::nls.lm.control(maxiter = 500)
    )
    if (!isTRUE(curve_fit$convInfo$isConv)) stop("A separated curve fit did not converge.")
    mixture_results$n[k] <- nrow(group_data)
    mixture_results$pka[k] <- unname(coef(curve_fit)["pKa"])
    mixture_results$se_pka[k] <- unname(summary(curve_fit)$coefficients["pKa", "Std. Error"])
    mixture_fits[[k]] <- curve_fit
  }

  # 3. Ordinary, unweighted YS regression on the three fitted pKa values.
  mixture_results$inverse_epsilon <- 1 / mixtures$epsilon
  mixture_results$corrected_pka <- mixture_results$pka + log10(mixtures$water_activity)
  line_fit <- lm(corrected_pka ~ inverse_epsilon, data = mixture_results)
  B <- unname(coef(line_fit)[1])
  A <- unname(coef(line_fit)[2])
  pka_water <- B + A / epsilon_water

  # 4. OLS extrapolation weights: pka_water = sum(weight * corrected_pka).
  # They show explicitly how each mixture's uncertainty enters the water value.
  x <- mixture_results$inverse_epsilon
  x_water <- 1 / epsilon_water
  x_mean <- mean(x)
  sum_squares_x <- sum((x - x_mean)^2)
  water_weights <- 1 / n_mixtures + (x_water - x_mean) * (x - x_mean) / sum_squares_x
  leverage <- 1 / n_mixtures + (x - x_mean)^2 / sum_squares_x
  pka_variance <- mixture_results$se_pka^2

  # The expected line RSS due to measurement alone is sum((1-h_k) * variance_k).
  # Three mixtures minus two line coefficients leave one residual degree of freedom.
  line_rss <- sum(residuals(line_fit)^2)
  expected_rss <- sum((1 - leverage) * pka_variance)
  tau_squared <- max(0, (line_rss - expected_rss) / (n_mixtures - 2))
  total_variance <- pka_variance + tau_squared
  se_water <- sqrt(sum(water_weights^2 * total_variance))           # SE_c: report
  se_measurement_only <- sqrt(sum(water_weights^2 * pka_variance))  # SE_a: reference
  se_ols <- unname(predict(line_fit,
                          newdata = data.frame(inverse_epsilon = x_water),
                          se.fit = TRUE)$se.fit)                  # SE_b: reference

  # 5. The same propagated variances give the covariance of the line's (B, A).
  # Each coefficient is a weighted sum of the three corrected pKa estimates.
  slope_weights <- (x - x_mean) / sum_squares_x
  intercept_weights <- 1 / n_mixtures - x_mean * slope_weights
  variance_A <- sum(slope_weights^2 * total_variance)
  variance_B <- sum(intercept_weights^2 * total_variance)
  covariance_AB <- sum(slope_weights * intercept_weights * total_variance)
  covariance <- matrix(c(variance_B, covariance_AB, covariance_AB, variance_A),
                       nrow = 2, dimnames = list(c("B", "A"), c("B", "A")))
  q_statistic <- line_rss / expected_rss
  p_value <- pchisq(q_statistic, df = 1, lower.tail = FALSE)

  # fit is the OLS line; mixture_fits holds the three separate titration fits.
  return(list(A = A, B = B, pka_water = pka_water, se_water = se_water,
              covariance = covariance, se_measurement_only = se_measurement_only,
              se_ols = se_ols, tau_squared = tau_squared,
              q_statistic = q_statistic, p_value = p_value,
              r_squared = summary(line_fit)$r.squared,
              mixture_results = mixture_results, mixture_fits = mixture_fits,
              fit = line_fit))
}

# =========== ACN concentration model: separate curve fits, then an OLS line ===========
# Input columns are pH, ChemShift, and Condition, as in fit_concentration().
# Stage 1 estimates a separate pKa and its SE for each mixture.
# Stage 2 fits pKa = beta0 + beta1 * acn_fraction by ordinary LS.
# Pure water has acn_fraction = 0, so beta0 is the aqueous pKa.
# The reported aqueous SE is SE_c; SE_a and SE_b are returned as references.

fit_concentration_separated <- function(data) {
  # 1. Identify the three mixtures; ACN fractions are 0.40, 0.50, and 0.60.
  constants <- get_constants()
  mixtures <- constants[constants$condition != "pure water", ]
  n_mixtures <- nrow(mixtures)

  stopifnot(all(c("pH", "ChemShift", "Condition") %in% names(data)))
  if (any(!is.finite(data$pH)) || any(!is.finite(data$ChemShift))) {
    stop("pH and ChemShift must contain finite numbers; no rows are dropped.")
  }
  data$grp <- match(data$Condition, mixtures$condition)
  if (anyNA(data$grp)) stop("Condition must be 40% ACN, 50% ACN, or 60% ACN.")

  # 2. Fit all THREE parameters within each mixture, including both plateaus.
  # These are the same separate titration fits used by the separated YS method.
  mixture_results <- mixtures
  mixture_results$n <- integer(n_mixtures)
  mixture_results$pka <- numeric(n_mixtures)
  mixture_results$se_pka <- numeric(n_mixtures)
  mixture_fits <- vector("list", n_mixtures)
  names(mixture_fits) <- mixtures$condition

  for (k in 1:n_mixtures) {
    group_data <- data[data$grp == k, ]
    if (nrow(group_data) < 4) stop("Each mixture needs at least four observations.")
    dA_start <- group_data$ChemShift[which.max(group_data$pH)]
    dHA_start <- group_data$ChemShift[which.min(group_data$pH)]
    curve_fit <- minpack.lm::nlsLM(
      ChemShift ~ (dA * 10^pH + dHA * 10^pKa) / (10^pH + 10^pKa),
      data = group_data,
      start = list(dA = unname(dA_start), dHA = unname(dHA_start),
                   pKa = unname(mean(group_data$pH))),
      control = minpack.lm::nls.lm.control(maxiter = 500)
    )
    if (!isTRUE(curve_fit$convInfo$isConv)) stop("A separated curve fit did not converge.")
    mixture_results$n[k] <- nrow(group_data)
    mixture_results$pka[k] <- unname(coef(curve_fit)["pKa"])
    mixture_results$se_pka[k] <- unname(summary(curve_fit)$coefficients["pKa", "Std. Error"])
    mixture_fits[[k]] <- curve_fit
  }

  # 3. Ordinary, unweighted regression of the three pKa values on ACN fraction.
  line_fit <- lm(pka ~ acn_fraction, data = mixture_results)
  beta0 <- unname(coef(line_fit)[1])
  beta1 <- unname(coef(line_fit)[2])
  pka_water <- beta0

  # 4. At C = 0 the extrapolation weights are the intercept weights.
  x <- mixture_results$acn_fraction
  x_mean <- mean(x)
  sum_squares_x <- sum((x - x_mean)^2)
  water_weights <- 1 / n_mixtures - x_mean * (x - x_mean) / sum_squares_x
  leverage <- 1 / n_mixtures + (x - x_mean)^2 / sum_squares_x
  pka_variance <- mixture_results$se_pka^2

  # Subtract expected measurement scatter before estimating extra line variation.
  # Three mixtures minus two line coefficients leave one residual degree of freedom.
  line_rss <- sum(residuals(line_fit)^2)
  expected_rss <- sum((1 - leverage) * pka_variance)
  tau_squared <- max(0, (line_rss - expected_rss) / (n_mixtures - 2))
  total_variance <- pka_variance + tau_squared
  se_water <- sqrt(sum(water_weights^2 * total_variance))           # SE_c: report
  se_measurement_only <- sqrt(sum(water_weights^2 * pka_variance))  # SE_a: reference
  se_ols <- unname(predict(line_fit,
                          newdata = data.frame(acn_fraction = 0),
                          se.fit = TRUE)$se.fit)                  # SE_b: reference

  # 5. Propagate the same variances to both line coefficients, beta0 and beta1.
  slope_weights <- (x - x_mean) / sum_squares_x
  intercept_weights <- water_weights
  variance_beta0 <- sum(intercept_weights^2 * total_variance)
  variance_beta1 <- sum(slope_weights^2 * total_variance)
  covariance_beta0_beta1 <- sum(intercept_weights * slope_weights * total_variance)
  covariance <- matrix(c(variance_beta0, covariance_beta0_beta1,
                         covariance_beta0_beta1, variance_beta1),
                       nrow = 2, dimnames = list(c("beta0", "beta1"), c("beta0", "beta1")))
  q_statistic <- line_rss / expected_rss
  p_value <- pchisq(q_statistic, df = 1, lower.tail = FALSE)

  # fit is the OLS line; mixture_fits holds the three separate titration fits.
  return(list(beta0 = beta0, beta1 = beta1, pka_water = pka_water, se_water = se_water,
              covariance = covariance, se_measurement_only = se_measurement_only,
              se_ols = se_ols, tau_squared = tau_squared,
              q_statistic = q_statistic, p_value = p_value,
              r_squared = summary(line_fit)$r.squared,
              mixture_results = mixture_results, mixture_fits = mixture_fits,
              fit = line_fit))
}
