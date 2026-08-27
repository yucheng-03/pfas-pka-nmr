# ==============================================================
# figures.R -- the two figure types, written once for both models.
#
# fig_titration()    raw chemical shift vs Adjusted pH, one curve
#                    per mixture.
# fig_extrapolation() two panels sharing the same axes:
#   left  = joint fit line + aqueous star with its 95% CI from the
#           joint parameter covariance;
#   right = OLS line through the three per-mixture estimates +
#           aqueous star with its 95% CI from SE_c.
#   Green squares (both panels) are the per-mixture estimates with
#   +/- 1.96 SE bars; on Model-1 axes the known offset is removed so
#   the plotted quantity is linear in the covariate.
# All intervals use the normal factor 1.96.
# ==============================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

fig_titration <- function(data, title, xlab = "Adjusted pH") {
  ggplot(data[order(data$pH), ], aes(x = pH, y = ChemShift, color = Condition)) +
    geom_line(linewidth = 0.6, alpha = 0.5) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = CONDITION_COLORS, drop = FALSE) +
    labs(title = title, x = xlab,
         y = expression(""^19 * F ~ "Chemical Shift (ppm)"), color = "Solvent") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14))
}

fig_extrapolation <- function(res, title) {
  sp <- res$spec
  ci_pts <- 1.96 * res$se_k
  ci_j   <- 1.96 * res$se_joint
  ci_c   <- 1.96 * res$se_c
  span   <- diff(range(c(sp$x, sp$x0)))
  xlim   <- c(min(sp$x, sp$x0) - 0.17 * span, max(sp$x, sp$x0) + 0.17 * span)
  lo <- c(res$y - ci_pts, res$aq_joint - ci_j, res$aq_two - ci_c)
  hi <- c(res$y + ci_pts, res$aq_joint + ci_j, res$aq_two + ci_c)
  ylim <- c(min(lo) - 0.1, max(hi) + 0.3)
  xs   <- seq(xlim[1], xlim[2], length.out = 200)
  cap  <- 0.09 * span                          # error-bar cap width
  pts  <- data.frame(x = sp$x, y = res$y,
                     ymin = res$y - ci_pts, ymax = res$y + ci_pts)
  line_j <- data.frame(x = xs, y = res$theta_joint[1] + res$theta_joint[2] * xs)
  line_t <- data.frame(x = xs, y = res$theta_two[1]   + res$theta_two[2]   * xs)

  base <- function(p, panel_title) p +
    geom_errorbar(data = pts, aes(x = x, ymin = ymin, ymax = ymax),
                  width = cap, color = "#2CA02C", linewidth = 0.7) +
    geom_point(data = pts, aes(x, y), shape = 15, size = 3, color = "#2CA02C") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = panel_title, x = sp$xlab, y = sp$ylab) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          panel.grid.minor = element_blank())

  p_left <- base(ggplot() +
      geom_line(data = line_j, aes(x, y), color = "#D62728", linewidth = 0.9),
      "Joint fit (aggregated)") +
    annotate("errorbar", x = sp$x0, ymin = res$aq_joint - ci_j,
             ymax = res$aq_joint + ci_j, width = cap, color = "#1F77B4") +
    annotate("point", x = sp$x0, y = res$aq_joint, shape = 8, size = 5,
             color = "#1F77B4") +
    annotate("label", x = xlim[1] + 0.02 * span, y = ylim[2] - 0.05,
             label = sprintf("R^2 = %.3f (points vs joint line)\nAqueous pKa = %.3f +/- %.3f (95%% CI)",
                             res$R2_joint, res$aq_joint, ci_j),
             hjust = 0, vjust = 1, size = 3.2, lineheight = 1.05,
             fill = scales::alpha("white", 0.75))

  p_right <- base(ggplot() +
      geom_line(data = line_t, aes(x, y), color = "#1F77B4",
                linewidth = 0.9, linetype = "dashed"),
      "Two-step fit (OLS)") +
    annotate("errorbar", x = sp$x0, ymin = res$aq_two - ci_c,
             ymax = res$aq_two + ci_c, width = cap, color = "black",
             linewidth = 0.9) +
    annotate("point", x = sp$x0, y = res$aq_two, shape = 8, size = 5,
             color = "#9467BD") +
    annotate("label", x = xlim[1] + 0.02 * span, y = ylim[2] - 0.05,
             label = sprintf("R^2 = %.3f (OLS)\nAqueous pKa = %.3f +/- %.3f (95%% CI)\ntau_hat = %.3f,  Q p = %.2g",
                             res$R2_ols, res$aq_two, ci_c, sqrt(res$tau2), res$p_Q),
             hjust = 0, vjust = 1, size = 3.0, lineheight = 1.05,
             fill = scales::alpha("white", 0.75))

  p_left + p_right + plot_annotation(
    title = title,
    caption = paste0(
      "Point error bars (both panels) = 95% CI from each mixture's own NLS SE.  ",
      "Left: joint line; blue star = aqueous extrapolation with the joint-fit CI.\n",
      "Right: OLS line through the three mixture estimates; the star carries the 95% CI from ",
      "SE_c^2 = g'(Z'Z)^-1 Z'[diag(nu) + tau2 I]Z(Z'Z)^-1 g -- measurement + lack-of-fit."),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.caption = element_text(size = 8.5, color = "grey30",
                                              hjust = 0.5)))
}
