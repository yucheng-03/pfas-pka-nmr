# ==============================================================
# constants.R -- solvent constants and the two model specifications.
#
# Both solvent models have the same algebraic form
#
#     pKa_k = theta' d_k + offset_k ,      theta = (theta_1, theta_2)'
#     pKa_aqueous = theta' g
#
# and differ only in the covariate d_k, the known offset, and the
# read-off point g.  Every downstream function takes a `spec` built
# here, so the two models share one implementation.
#
#   Yasuda-Shedlovsky : d_k = (1, 1/eps_k)',  offset_k = -log10(a_k),
#                       g = (1, 1/78.33)'
#   ACN concentration : d_k = (1, C_k)',      offset_k = 0,
#                       g = (1, 0)'
#
# Solvent constants are fixed here and are NOT read from the
# spreadsheets (the 'water' column in those files is outdated).
# Source: collaborator email, 2026-07-15 (literature dielectric
# constants and water activities).
# ==============================================================

CONDITIONS <- c("40% ACN", "50% ACN", "60% ACN")

EPSILON       <- c(63.50, 58.71, 53.91)   # dielectric constant per mixture
WATER_ACT     <- c(0.913, 0.908, 0.901)   # water activity per mixture
ACN_FRACTION  <- c(0.40, 0.50, 0.60)      # ACN volume fraction per mixture
EPSILON_WATER <- 78.33                    # pure water
WATER_ACT_W   <- 1.00                     # pure water

CONDITION_COLORS <- c("40% ACN" = "#00A087",
                      "50% ACN" = "#4DBBD5",
                      "60% ACN" = "#E64B35")

# --------------------------------------------------------------
# model_spec(): assemble one model.  Fields used downstream:
#   name    label for titles
#   x       covariate value per mixture (length 3)
#   x0      covariate value of pure water (the read-off point)
#   offset  known additive term of pKa_k per mixture (length 3)
#   g       read-off vector c(1, x0)
#   Z       3 x 2 design matrix of rows d_k'
#   ell     design coefficients: aqueous = sum_k ell_k * y_k
#   h       leverages of Z (diag of the hat matrix)
#   xlab    axis label for the extrapolation figure
#   ylab    y-axis label; the plotted quantity is pKa_tilde_k - offset_k
# --------------------------------------------------------------
model_spec <- function(model = c("ys", "conc")) {
  model <- match.arg(model)
  if (model == "ys") {
    x      <- 1 / EPSILON
    x0     <- 1 / EPSILON_WATER
    offset <- -log10(WATER_ACT)          # pKa_k = B + A x_k - log10 a_k
    name   <- "Yasuda-Shedlovsky"
    xlab   <- expression("1 / "*epsilon)
    ylab   <- expression(widetilde(pKa)[k]*" + log"[10]*" "*italic(a)["H"[2]*"O"])
  } else {
    x      <- ACN_FRACTION
    x0     <- 0
    offset <- rep(0, 3)                  # pKa_k = beta0 + beta1 C_k
    name   <- "ACN concentration"
    xlab   <- "ACN volume fraction C"
    ylab   <- expression(widetilde(pKa)[k])
  }
  Z    <- cbind(1, x)
  ZtZi <- solve(t(Z) %*% Z)
  g    <- c(1, x0)
  list(model = model, name = name, x = x, x0 = x0, offset = offset,
       g = g, Z = Z, ZtZi = ZtZi,
       ell = as.numeric(t(g) %*% ZtZi %*% t(Z)),   # aqueous = sum ell_k y_k
       h   = diag(Z %*% ZtZi %*% t(Z)),
       xlab = xlab, ylab = ylab)
}

# Quadratic form g' V g -- every standard error in this project is
# an instance of it (V = covariance of the two model coefficients).
aqueous_se <- function(V, g) sqrt(as.numeric(t(g) %*% V %*% g))
