# ==============================================================
# verify.R -- self-checks.  Run after run_all.R:
#
#   Rscript verify.R
#
# 0. Model configuration: each spec's covariate, offset and read-off
#    point are what the analysis claims they are (a wrong x0 or a
#    wrong offset would otherwise propagate consistently and stay
#    invisible to the identity checks below).
# 1. Algebraic identities that must hold, for every dataset and both
#    models:
#      a) the aqueous estimate equals g' theta;
#      b) SE_c equals sqrt(g' V_c g);
#      c) under the concentration model (read-off at C = 0) the
#         aqueous SE equals the SE of the intercept;
#      d) tau2_hat > 0 exactly when Q > 1;
#      e) the two models share the same stage-1 estimates (those are
#         fitted per mixture and cannot depend on the solvent model);
#      f) the analytic covariance agrees with the optimizer's own
#         standard error for the aqueous parameter.  The joint fit is
#         run in coordinates centered at the read-off point, so that
#         parameter IS the aqueous value and its reported SE is an
#         independent estimate of g' V g -- this is what catches a
#         transposed covariance or a mis-contracted read-off vector.
# 2. Agreement with the reference CSV of the delivered report.
# Every check stops the script on failure.
# ==============================================================
for (f in c("constants", "data", "fit", "inference"))
  source(file.path("R", paste0(f, ".R")))

stopifnot(file.exists(file.path("output", "summary_all.csv")))
res <- read.csv(file.path("output", "summary_all.csv"), stringsAsFactors = FALSE)
ok  <- function(name, value, tol = 1e-9) {
  cat(sprintf("%-58s %.2e %s\n", name, value, if (value <= tol) "OK" else "FAIL"))
  if (value > tol) stop("check failed: ", name)
}

specs <- list(ys = model_spec("ys"), conc = model_spec("conc"))

# (0) the model configuration is what the analysis claims
ok("ys spec: covariate equals 1/eps",
   max(abs(specs$ys$x - 1 / EPSILON)))
ok("ys spec: read-off point equals 1/eps_water",
   abs(specs$ys$x0 - 1 / EPSILON_WATER))
ok("ys spec: offset equals -log10(water activity)",
   max(abs(specs$ys$offset + log10(WATER_ACT))))
ok("conc spec: covariate equals the ACN volume fraction",
   max(abs(specs$conc$x - ACN_FRACTION)))
ok("conc spec: read-off point equals 0 and offset equals 0",
   max(abs(specs$conc$x0), max(abs(specs$conc$offset))))
for (m in names(specs))
  ok(sprintf("%s spec: g equals (1, x0) and sum(ell) = 1, sum(ell*x) = x0", m),
     max(abs(specs[[m]]$g - c(1, specs[[m]]$x0)),
         abs(sum(specs[[m]]$ell) - 1),
         abs(sum(specs[[m]]$ell * specs[[m]]$x) - specs[[m]]$x0)))

# (a) aqueous estimate = g' theta, both routes
d <- max(sapply(seq_len(nrow(res)), function(i) {
  sp <- specs[[res$model[i]]]
  max(abs(res$aq_joint[i]   - sum(sp$g * c(res$theta1[i],     res$theta2[i]))),
      abs(res$aq_twostep[i] - sum(sp$g * c(res$theta1_two[i], res$theta2_two[i]))))
}))
ok("aqueous estimate equals g' theta", d)

# (b) SE_c = sqrt(g' V_c g) rebuilt from the reported parameter SEs
d <- max(sapply(seq_len(nrow(res)), function(i) {
  sp <- specs[[res$model[i]]]
  V <- matrix(c(res$SE_theta1_two[i]^2,
                res$corr_theta_two[i]*res$SE_theta1_two[i]*res$SE_theta2_two[i],
                res$corr_theta_two[i]*res$SE_theta1_two[i]*res$SE_theta2_two[i],
                res$SE_theta2_two[i]^2), 2, 2)
  abs(res$SE_c[i] - aqueous_se(V, sp$g))
}))
ok("SE_c equals sqrt(g' V_c g) from the reported entries", d, 1e-8)

# (c) concentration model: aqueous SE is the intercept SE
cc <- res[res$model == "conc", ]
ok("concentration model: SE_joint equals SE(theta1)",
   max(abs(cc$SE_joint - cc$SE_theta1)))
ok("concentration model: SE_c equals SE(theta1_two)",
   max(abs(cc$SE_c - cc$SE_theta1_two)))

# (d) tau2 > 0 exactly when Q > 1
ok("tau2_hat > 0 exactly when Q > 1",
   as.numeric(any((res$tau2_hat > 0) != (res$Q > 1))))

# (e) stage 1 does not depend on the solvent model
ys <- res[res$model == "ys", ]
stopifnot(identical(ys$dataset, cc$dataset))
ok("stage-1 estimates identical across the two models",
   max(abs(as.matrix(ys[, c("pKa1","pKa2","pKa3","SE1","SE2","SE3")]) -
           as.matrix(cc[, c("pKa1","pKa2","pKa3","SE1","SE2","SE3")]))))

# (f) analytic covariance vs the optimizer's own SE for the aqueous
# parameter.  Independent quantities: one comes from sigma^2 M^-1 with
# analytic weights, the other from the fit's own Jacobian in centered
# coordinates.  A transposed V or a wrong g breaks this immediately.
ok("analytic vs optimizer covariance (relative difference)",
   max(abs(res$SE_joint - res$SE_joint_numeric) / res$SE_joint), 0.05)

# ---- agreement with the delivered report ----
# summary_all_from_code.csv is the copy of this pipeline's output that
# accompanies the delivered report; the two must be identical.
ref <- file.path("..", "Report_2026-08-26", "summary_all_from_code.csv")
if (file.exists(ref)) {
  r <- read.csv(ref, stringsAsFactors = FALSE)
  key <- function(d) paste(d$dataset, d$model)
  n <- res[match(key(r), key(res)), ]
  cols <- intersect(names(r), names(res))
  cols <- cols[sapply(r[cols], is.numeric)]
  cat("\n-- agreement with the delivered report --\n")
  ok("report reference: max |difference| over all numeric columns",
     max(abs(as.matrix(r[, cols]) - as.matrix(n[, cols]))), 1e-10)
} else {
  cat("\n(reference CSV of the delivered report not reachable; skipped)\n")
}
cat("\nall checks passed\n")
