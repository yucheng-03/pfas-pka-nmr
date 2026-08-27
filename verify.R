# ==============================================================
# verify.R -- self-checks.  Run after run_all.R:
#
#   Rscript verify.R
#
# 1. Algebraic identities that must hold exactly, for every dataset
#    and both models:
#      a) the aqueous estimate equals g' theta;
#      b) SE_c equals sqrt(g' V_c g);
#      c) under the concentration model (read-off at C = 0) the
#         aqueous SE equals the SE of the intercept;
#      d) tau2_hat > 0 exactly when Q > 1;
#      e) the two models share the same stage-1 estimates (those are
#         fitted per mixture and cannot depend on the solvent model).
# 2. Agreement with the delivered report numbers, if the reference
#    CSV of the report folder is reachable (optional).
# Any failure stops the script with an error.
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

# (f) analytic vs numerical joint covariance (informational, loose)
rel <- max(abs(res$SE_joint - res$SE_joint_numeric) / res$SE_joint)
cat(sprintf("%-58s %.2e %s\n", "analytic vs optimizer covariance (relative)",
            rel, if (rel < 0.05) "OK" else "CHECK"))

# ---- optional: agreement with the delivered report ----
ref <- file.path("..", "Report_2026-08-26", "summary_part1_YS.csv")
if (file.exists(ref)) {
  cmp <- function(path, model, keys) {
    r <- read.csv(path, stringsAsFactors = FALSE)
    n <- res[res$model == model & res$dataset %in% r$dataset, ]
    n <- n[match(r$dataset, n$dataset), ]
    cols <- intersect(c("pKa1","pKa2","pKa3","SE1","SE2","SE3",
                        "aq_joint","SE_joint","aq_twostep","SE_a","SE_b","SE_c",
                        "tau2_hat","Q","p_Q","R2_joint","R2_ols"), names(r))
    max(abs(as.matrix(r[, cols]) - as.matrix(n[, cols])))
  }
  cat("\n-- agreement with the delivered report numbers --\n")
  for (p in 1:2) for (m in c("YS", "Conc")) {
    f <- file.path("..", "Report_2026-08-26", sprintf("summary_part%d_%s.csv", p, m))
    if (file.exists(f))
      cat(sprintf("%-58s %.2e\n", sprintf("part %d, %s: max |difference|", p, m),
                  cmp(f, tolower(ifelse(m == "YS", "ys", "conc")))))
  }
} else {
  cat("\n(reference report CSVs not found; skipping that comparison)\n")
}
cat("\nall checks passed\n")
