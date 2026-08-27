# ==============================================================
# tables.R -- LaTeX table fragments for the report.
#
# Every fragment prints the model coefficients, their standard
# errors and their correlation, so that each aqueous standard error
# quoted can be recomputed by the reader from the entries shown via
#     Var(aqueous) = g' V g .
# Symbols are model-dependent: (B, A) for Yasuda-Shedlovsky,
# (beta_0, beta_1) for the concentration model.
# ==============================================================

.fmt_p <- function(p) if (p >= 0.01) sprintf("$%.2f$", p) else {
  e <- floor(log10(p)); sprintf("$%.1f\\times 10^{%d}$", p / 10^e, e) }
.fmt_tau <- function(v) if (v == 0) "0" else sprintf("%.4f", v)
.tex_esc <- function(x) gsub("_", "\\\\_", x)

# symbol names and the read-off expression per model
.symbols <- function(model) {
  if (model == "ys") {
  list(t1 = "\\hat B", t2 = "\\hat A", u1 = "\\tilde B", u2 = "\\tilde A",
       corr = "\\operatorname{corr}(\\hat A,\\hat B)",
       corr2 = "\\operatorname{corr}(\\tilde A,\\tilde B)",
       readoff = "\\hat B + \\hat A x_w", readoff2 = "\\tilde B + \\tilde A x_w",
       gtext = "g = (1,\\ 0.012767)'")
  } else {
  list(t1 = "\\hat\\beta_0", t2 = "\\hat\\beta_1",
       u1 = "\\tilde\\beta_0", u2 = "\\tilde\\beta_1",
       corr = "\\operatorname{corr}(\\hat\\beta_0,\\hat\\beta_1)",
       corr2 = "\\operatorname{corr}(\\tilde\\beta_0,\\tilde\\beta_1)",
       readoff = "\\hat\\beta_0", readoff2 = "\\tilde\\beta_0",
       gtext = "g = (1,\\ 0)'")
  }
}

# per-dataset, per-model table
table_dataset <- function(r, path) {
  s <- .symbols(r$model)
  dig <- if (r$model == "ys") c(4, 2) else c(3, 3)   # theta1, theta2 digits
  f2 <- function(v, d) sprintf(paste0("%.", d, "f"), v)
  writeLines(c(
    "\\begin{center}\\small",
    "\\begin{tabular}{lccc}", "\\toprule",
    " & 40\\% ACN & 50\\% ACN & 60\\% ACN\\\\", "\\midrule",
    sprintf("Green points $\\widetilde{\\mathrm{pKa}}_k$ (SE) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)\\\\",
            r$pKa1, r$SE1, r$pKa2, r$SE2, r$pKa3, r$SE3),
    sprintf("Titration noise $\\hat\\sigma_k$ (ppm) & %.4f & %.4f & %.4f\\\\",
            r$sig1, r$sig2, r$sig3),
    sprintf("Points per mixture $n_k$ & %d & %d & %d\\\\", r$n1, r$n2, r$n3),
    "\\midrule",
    sprintf("\\multicolumn{4}{l}{\\textbf{Joint fit:} $%s = %s\\ (%s)$, $%s = %s\\ (%s)$, $%s = %.3f$, $\\hat\\sigma = %.4f$ ppm}\\\\",
            s$t1, f2(r$theta1, dig[1]), f2(r$SE_theta1, dig[1]),
            s$t2, f2(r$theta2, dig[2]), f2(r$SE_theta2, dig[2]),
            s$corr, r$corr_theta, r$sigma),
    sprintf("\\multicolumn{4}{l}{\\quad aqueous $\\mathrm{pKa} = %s = %.3f$;\\quad SE $= \\sqrt{g'\\hat V g} = %.3f$;\\quad $R^2 = %.3f$}\\\\",
            s$readoff, r$aq_joint, r$SE_joint, r$R2_joint),
    "\\midrule",
    sprintf("\\multicolumn{4}{l}{\\textbf{Two-step:} $%s = %s\\ (%s)$, $%s = %s\\ (%s)$, $%s = %.3f$}\\\\",
            s$u1, f2(r$theta1_two, dig[1]), f2(r$SE_theta1_two, dig[1]),
            s$u2, f2(r$theta2_two, dig[2]), f2(r$SE_theta2_two, dig[2]),
            s$corr2, r$corr_theta_two),
    sprintf("\\multicolumn{4}{l}{\\quad aqueous $\\mathrm{pKa} = %s = %.3f$;\\quad $\\mathrm{SE}_c = \\sqrt{g'\\tilde V_c g} = %.3f$;\\quad $R^2 = %.3f$}\\\\",
            s$readoff2, r$aq_twostep, r$SE_c, r$R2_ols),
    sprintf("\\multicolumn{4}{l}{\\quad $\\hat\\tau^2 = %s$;\\quad $Q = %.2f$ (%s)}\\\\",
            .fmt_tau(r$tau2_hat), r$Q, .fmt_p(r$p_Q)),
    "\\bottomrule", "\\end{tabular}", "\\end{center}"), path)
}

# per-dataset comparison of the two models
table_models <- function(ry, rc, path) writeLines(c(
  "\\begin{center}\\small", "\\begin{tabular}{lcc}", "\\toprule",
  "Model & aqueous $\\mathrm{pKa}$, joint (SE) & aqueous $\\mathrm{pKa}$, two-step (SE$_c$)\\\\",
  "\\midrule",
  sprintf("Yasuda--Shedlovsky & %.3f (%.3f) & %.3f (%.3f)\\\\",
          ry$aq_joint, ry$SE_joint, ry$aq_twostep, ry$SE_c),
  sprintf("ACN concentration & %.3f (%.3f) & %.3f (%.3f)\\\\",
          rc$aq_joint, rc$SE_joint, rc$aq_twostep, rc$SE_c),
  "\\bottomrule", "\\end{tabular}", "\\end{center}"), path)

# overview of one part: all datasets, both models
table_overview <- function(ys, cc, files, path) writeLines(c(
  "\\begin{center}\\small", "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llcccc}", "\\toprule",
  " & & \\multicolumn{2}{c}{Yasuda--Shedlovsky} & \\multicolumn{2}{c}{ACN concentration}\\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  "Dataset & Source spreadsheet & joint (SE) & two-step (SE$_c$) & joint (SE) & two-step (SE$_c$)\\\\",
  "\\midrule",
  sapply(seq_len(nrow(ys)), function(i) sprintf(
    "%s & \\texttt{%s} & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)\\\\",
    ys$dataset[i], .tex_esc(files[i]),
    ys$aq_joint[i], ys$SE_joint[i], ys$aq_twostep[i], ys$SE_c[i],
    cc$aq_joint[i], cc$SE_joint[i], cc$aq_twostep[i], cc$SE_c[i])),
  "\\bottomrule", "\\end{tabular}}", "\\end{center}"), path)

# lack-of-fit diagnostics of one part, both models
table_diagnostics <- function(ys, cc, path) writeLines(c(
  "\\begin{center}\\small", "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{3}{c}{Yasuda--Shedlovsky} & \\multicolumn{3}{c}{ACN concentration}\\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "Dataset & $\\hat\\tau^2$ & $Q$ & $p$ & $\\hat\\tau^2$ & $Q$ & $p$\\\\", "\\midrule",
  sapply(seq_len(nrow(ys)), function(i) sprintf(
    "%s & %s & %.2f & %s & %s & %.2f & %s\\\\", ys$dataset[i],
    .fmt_tau(ys$tau2_hat[i]), ys$Q[i], .fmt_p(ys$p_Q[i]),
    .fmt_tau(cc$tau2_hat[i]), cc$Q[i], .fmt_p(cc$p_Q[i]))),
  "\\bottomrule", "\\end{tabular}", "\\end{center}"), path)
