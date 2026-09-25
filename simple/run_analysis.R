# Run from the pfas-pka-nmr repository root:
# Rscript simple/run_analysis.R
# This script runs the eight workbooks below, saving four results and three figures each.
# All four fits include the supplied per-observation f; results go to simple/output/with_f/.
# Edit the dataset table to choose files, display labels, and chemical-shift columns.
# Required packages: readxl, minpack.lm, ggplot2, patchwork.
# Put the workbooks in data/ as described in data/README.md.

# =========== 1. List the workbooks and columns ===========
library(readxl)
library(ggplot2)
library(patchwork)

analysis_dir <- "simple"
if (!file.exists(file.path(analysis_dir, "models.R"))) {
  stop("Run this script from the pfas-pka-nmr repository root.")
}
source(file.path(analysis_dir, "models.R"))

data_dir <- "data"
output_dir <- file.path(analysis_dir, "output", "with_f")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# One row per distinct workbook; IDs distinguish different peaks and data versions.
datasets <- data.frame(
  dataset_id = c("5_3_FTCA_main", "7_3_FTCA_newest", "8_2_FTCA", "PFBA", "PFOA",
                 "5_3_FTCA_otherCF2", "7_3_FTCA_old", "7_3_FTCA_otherCF2"),
  compound = c("5:3 FTCA (main)", "7:3 FTCA (newest)", "8:2 FTCA", "PFBA", "PFOA",
               "5:3 FTCA (other CF2)", "7:3 FTCA (old)", "7:3 FTCA (other CF2)"),
  source_file = c("5_3_FTCA_pKa_Dielectric.xlsx",
                  "7_3_FTCA_pKa_Dielectric_Newest_Data.xlsx",
                  "8_2_FTCA_pKa_Dielectric.xlsx",
                  "PFBA_pKa_Dielectric.xlsx",
                  "PFOA_pKa_Dielectric.xlsx",
                  "5_3_FTCA_pKa_other_CF2.xlsx",
                  "7_3_FTCA_pKa_Dielectric_Old_Data.xlsx",
                  "7_3_FTCA_pKa_other_CF2.xlsx"),
  shift_column = c("chemical shift (ppm)", "chemical shift (ppm)",
                   "chemical shift (ppm)", "chemical shift (ppm)",
                   "chemical shift (ppm)", "deltai_calc",
                   "chemical shift (ppm)", "chemical shift (ppm)"),
  data_note = c("", "", "", "", "", "",
                "Old workbook pH correction needs confirmation; Adjusted pH used as supplied.",
                "")
)

# Names below use lowercase, matching the header normalization during reading.
# Pass workbook Adjusted pH unchanged; models.R includes f explicitly in the mean.
# Do not also subtract log10(f) here, which would apply f twice.
ph_column <- "adjusted ph"
f_column <- "f"
constants <- get_constants()
conditions <- constants$condition[constants$condition != "pure water"]
ci_multiplier <- 1.96
all_results <- data.frame()

if (anyDuplicated(datasets$dataset_id) || anyDuplicated(datasets$source_file)) {
  stop("Each dataset ID and source workbook must appear only once in the dataset table.")
}

for (dataset_index in seq_len(nrow(datasets))) {
  dataset_id <- datasets$dataset_id[dataset_index]
  compound <- datasets$compound[dataset_index]
  data_file <- file.path(data_dir, datasets$source_file[dataset_index])
  shift_column <- datasets$shift_column[dataset_index]
  data_note <- datasets$data_note[dataset_index]
  cat("\n==========", dataset_index, "/", nrow(datasets), ":", compound, "==========\n")
  if (!file.exists(data_file)) stop("Data file not found: ", data_file)

  # =========== 2. Read 40%, 50%, and 60% ACN separately ===========
  titration_data <- data.frame()

  for (condition in conditions) {
    # Keep original headers so duplicate selected columns can be detected explicitly.
    sheet_data <- read_excel(data_file, sheet = condition, .name_repair = "minimal")
    names(sheet_data) <- tolower(trimws(names(sheet_data)))

    if (sum(names(sheet_data) == ph_column) != 1) {
      stop(condition, ": expected exactly one column named '", ph_column, "'.")
    }
    if (sum(names(sheet_data) == shift_column) != 1) {
      stop(condition, ": expected exactly one column named '", shift_column, "'.")
    }
    if (sum(names(sheet_data) == f_column) != 1) {
      stop(condition, ": expected exactly one column named '", f_column, "'.")
    }
    if (nrow(sheet_data) == 0) stop(condition, ": no observations found.")

    pH <- as.numeric(sheet_data[[ph_column]])
    chemical_shift <- as.numeric(sheet_data[[shift_column]])
    f <- as.numeric(sheet_data[[f_column]])
    if (any(!is.finite(pH)) || any(!is.finite(chemical_shift))) {
      stop(condition, ": missing or invalid pH/shift values; inspect the source sheet.")
    }
    if (any(!is.finite(f)) || any(f <= 0)) {
      stop(condition, ": f must contain positive finite numbers; inspect the source sheet.")
    }

    # Preserve every observation, including repeated pH values, in its original order.
    group_data <- data.frame(pH = pH, ChemShift = chemical_shift, Condition = condition, f = f)
    titration_data <- rbind(titration_data, group_data)
  }
  titration_data$Condition <- factor(titration_data$Condition, levels = conditions)

  # =========== 3. Inspect the data supplied to the model functions ===========
  cat("\nCompound:", compound, "\nWorkbook:", data_file, "\n")
  cat("Selected columns:", ph_column, "/", shift_column, "/", f_column, "\n\n")
  print(table(titration_data$Condition))
  print(head(titration_data))

  # =========== 4. Fit both models by both methods and compare results ===========
  # Each call receives the same observations and returns its own aqueous estimate and SE.
  # f is a fixed input in every fit; reported SEs do not include uncertainty in f.
  ys_aggregated <- fit_ys(titration_data)
  ys_separated <- fit_ys_separated(titration_data)
  concentration_aggregated <- fit_concentration(titration_data)
  concentration_separated <- fit_concentration_separated(titration_data)

  # Aggregated uses the full joint-NLS covariance; separated uses the report's SE_c.
  # Keep full precision in this table; the digits argument below only controls printing.
  comparison_table <- data.frame(
    model = c("YS", "YS", "Concentration", "Concentration"),
    method = c("Aggregated", "Separated", "Aggregated", "Separated"),
    pka_water = c(ys_aggregated$pka_water, ys_separated$pka_water,
                  concentration_aggregated$pka_water, concentration_separated$pka_water),
    se_water = c(ys_aggregated$se_water, ys_separated$se_water,
                 concentration_aggregated$se_water, concentration_separated$se_water),
    se_type = c("SE_joint", "SE_c", "SE_joint", "SE_c")
  )

  cat("\nAqueous pKa comparison:", compound, "\n")
  cat("All fits include the supplied f, treated as fixed.\n")
  cat("SE_c includes measurement uncertainty and estimated lack-of-fit.\n\n")
  print(comparison_table, row.names = FALSE, digits = 6)
  if (nzchar(data_note)) cat("Data note:", data_note, "\n")

  # Save full-precision values with their source and the same intervals used in the plots.
  comparison_table$ci_lower <- comparison_table$pka_water - ci_multiplier * comparison_table$se_water
  comparison_table$ci_upper <- comparison_table$pka_water + ci_multiplier * comparison_table$se_water
  comparison_table <- cbind(
    data.frame(dataset_id = dataset_id, compound = compound,
               source_file = datasets$source_file[dataset_index],
               ph_column = ph_column, shift_column = shift_column,
               f_column = f_column, f_included = TRUE,
               n = nrow(titration_data), data_note = data_note),
    comparison_table
  )
  write.csv(comparison_table, file.path(output_dir, paste0("Results_", dataset_id, ".csv")),
            row.names = FALSE)
  all_results <- rbind(all_results, comparison_table)

  # =========== 5. Plot the observed titration data ===========
  # Sort a copy for plotting; every original observation remains in titration_data.
  # The lines connect observed points within each mixture, as in the current report.
  plot_data <- titration_data[order(titration_data$Condition, titration_data$pH), ]
  ph_label <- if (ph_column == "adjusted ph") "Adjusted pH" else "pH"
  condition_colors <- c("40% ACN" = "#00A087", "50% ACN" = "#4DBBD5", "60% ACN" = "#E64B35")

  titration_plot <- ggplot(plot_data,
                           aes(x = pH, y = ChemShift, color = Condition, group = Condition)) +
    geom_line(linewidth = 0.6, alpha = 0.5) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = condition_colors, drop = FALSE) +
    labs(title = paste(compound, "- Titration data"),
         x = ph_label, y = expression(""^19 * F ~ "Chemical Shift (ppm)"),
         color = "Solvent",
         caption = paste0("Points are observations; lines connect observations within each mixture.",
                          if (nzchar(data_note)) paste0("\n", data_note) else "")) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14),
          plot.caption = element_text(size = 9, hjust = 0),
          panel.grid.minor = element_blank())

  titration_file <- file.path(output_dir, paste0("Fig_titration_", dataset_id, ".png"))
  ggsave(titration_file, plot = titration_plot, width = 7, height = 5, units = "in", dpi = 200)
  if (interactive()) print(titration_plot)
  cat("\nTitration figure saved:", titration_file, "\n")

  # =========== 6. YS plots: aggregated and separated side by side ===========
  # Both panels use the same independent mixture estimates as green squares.
  # The YS ordinate is pKa + log10(water activity); its fitted line is B + A / epsilon.
  # Keep the current report's normal-approximation intervals: estimate +/- 1.96 * SE.
  epsilon_water <- constants$epsilon[constants$condition == "pure water"]
  ys_points <- data.frame(
    x = ys_separated$mixture_results$inverse_epsilon,
    y = ys_separated$mixture_results$corrected_pka,
    se = ys_separated$mixture_results$se_pka
  )
  ys_points$lower <- ys_points$y - ci_multiplier * ys_points$se
  ys_points$upper <- ys_points$y + ci_multiplier * ys_points$se

  ys_aggregated_water <- data.frame(
    x = 1 / epsilon_water, y = ys_aggregated$pka_water,
    lower = ys_aggregated$pka_water - ci_multiplier * ys_aggregated$se_water,
    upper = ys_aggregated$pka_water + ci_multiplier * ys_aggregated$se_water
  )
  ys_separated_water <- data.frame(
    x = 1 / epsilon_water, y = ys_separated$pka_water,
    lower = ys_separated$pka_water - ci_multiplier * ys_separated$se_water,
    upper = ys_separated$pka_water + ci_multiplier * ys_separated$se_water
  )

  # Shared axes include all three mixture intervals and both aqueous intervals.
  ys_x_range <- range(c(ys_points$x, 1 / epsilon_water))
  ys_x_span <- diff(ys_x_range)
  ys_x_limits <- ys_x_range + c(-0.12, 0.12) * ys_x_span
  ys_y_range <- range(c(ys_points$lower, ys_points$upper,
                        ys_aggregated_water$lower, ys_aggregated_water$upper,
                        ys_separated_water$lower, ys_separated_water$upper))
  ys_y_limits <- ys_y_range + c(-0.10, 0.10) * diff(ys_y_range)
  ys_cap_width <- 0.07 * ys_x_span

  ys_lines <- data.frame(x = seq(ys_x_limits[1], ys_x_limits[2], length.out = 200))
  ys_lines$aggregated <- ys_aggregated$B + ys_aggregated$A * ys_lines$x
  ys_lines$separated <- ys_separated$B + ys_separated$A * ys_lines$x

  # Left panel: joint line and the joint-fit aqueous confidence interval.
  ys_aggregated_plot <- ggplot() +
    geom_line(data = ys_lines, aes(x = x, y = aggregated), color = "#D62728", linewidth = 0.9) +
    geom_errorbar(data = ys_points, aes(x = x, ymin = lower, ymax = upper),
                  width = ys_cap_width, color = "#2CA02C", linewidth = 0.7) +
    geom_point(data = ys_points, aes(x = x, y = y), color = "#2CA02C", shape = 15, size = 3) +
    geom_errorbar(data = ys_aggregated_water, aes(x = x, ymin = lower, ymax = upper),
                  width = ys_cap_width, color = "#1F77B4", linewidth = 0.9) +
    geom_point(data = ys_aggregated_water, aes(x = x, y = y), color = "#1F77B4", shape = 8, size = 4.5) +
    coord_cartesian(xlim = ys_x_limits, ylim = ys_y_limits) +
    labs(title = "Aggregated (joint NLS)",
         subtitle = sprintf("Aqueous pKa = %.3f; 95%% CI [%.3f, %.3f]",
                            ys_aggregated$pka_water, ys_aggregated_water$lower, ys_aggregated_water$upper),
         x = expression(1 / epsilon), y = expression(pK[a] + log[10](a[H[2]*O]))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          panel.grid.minor = element_blank())

  # Right panel: OLS line and the separated aqueous confidence interval from SE_c.
  ys_separated_plot <- ggplot() +
    geom_line(data = ys_lines, aes(x = x, y = separated),
              color = "#1F77B4", linewidth = 0.9, linetype = "dashed") +
    geom_errorbar(data = ys_points, aes(x = x, ymin = lower, ymax = upper),
                  width = ys_cap_width, color = "#2CA02C", linewidth = 0.7) +
    geom_point(data = ys_points, aes(x = x, y = y), color = "#2CA02C", shape = 15, size = 3) +
    geom_errorbar(data = ys_separated_water, aes(x = x, ymin = lower, ymax = upper),
                  width = ys_cap_width, color = "#9467BD", linewidth = 0.9) +
    geom_point(data = ys_separated_water, aes(x = x, y = y), color = "#9467BD", shape = 8, size = 4.5) +
    coord_cartesian(xlim = ys_x_limits, ylim = ys_y_limits) +
    labs(title = "Separated (OLS)",
         subtitle = sprintf("Aqueous pKa = %.3f; 95%% CI [%.3f, %.3f]",
                            ys_separated$pka_water, ys_separated_water$lower, ys_separated_water$upper),
         x = expression(1 / epsilon), y = expression(pK[a] + log[10](a[H[2]*O]))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          panel.grid.minor = element_blank())

  ys_plot <- ys_aggregated_plot + ys_separated_plot + plot_annotation(
    title = paste(compound, "- Yasuda-Shedlovsky"),
    caption = paste0(
      "Fits include each observation's supplied f; intervals treat f as fixed.\n",
      "Green squares: independent mixture pKa estimates + log10(water activity), with 95% intervals.\n",
      "Stars: aqueous extrapolations. Left interval: joint-NLS SE; right interval: SE_c (measurement + lack-of-fit).",
      if (nzchar(data_note)) paste0("\n", data_note) else ""),
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.caption = element_text(size = 9, hjust = 0))
  )
  ys_file <- file.path(output_dir, paste0("Fig_YS_", dataset_id, ".png"))
  ggsave(ys_file, plot = ys_plot, width = 12, height = 5.7, units = "in", dpi = 200)
  if (interactive()) print(ys_plot)
  cat("\nYS figure saved:", ys_file, "\n")

  # =========== 7. Concentration plots: aggregated and separated side by side ===========
  # Both panels use independent mixture pKa estimates, now plotted directly against C.
  # The fitted lines are beta0 + beta1 * C; the pure-water point is at C = 0.
  acn_fraction_water <- constants$acn_fraction[constants$condition == "pure water"]
  concentration_points <- data.frame(
    x = concentration_separated$mixture_results$acn_fraction,
    y = concentration_separated$mixture_results$pka,
    se = concentration_separated$mixture_results$se_pka
  )
  concentration_points$lower <- concentration_points$y - ci_multiplier * concentration_points$se
  concentration_points$upper <- concentration_points$y + ci_multiplier * concentration_points$se

  concentration_aggregated_water <- data.frame(
    x = acn_fraction_water, y = concentration_aggregated$pka_water,
    lower = concentration_aggregated$pka_water - ci_multiplier * concentration_aggregated$se_water,
    upper = concentration_aggregated$pka_water + ci_multiplier * concentration_aggregated$se_water
  )
  concentration_separated_water <- data.frame(
    x = acn_fraction_water, y = concentration_separated$pka_water,
    lower = concentration_separated$pka_water - ci_multiplier * concentration_separated$se_water,
    upper = concentration_separated$pka_water + ci_multiplier * concentration_separated$se_water
  )

  # Shared axes include all three mixture intervals and both aqueous intervals.
  concentration_x_range <- range(c(concentration_points$x, acn_fraction_water))
  concentration_x_span <- diff(concentration_x_range)
  concentration_x_limits <- concentration_x_range + c(-0.12, 0.12) * concentration_x_span
  concentration_y_range <- range(c(concentration_points$lower, concentration_points$upper,
                                   concentration_aggregated_water$lower, concentration_aggregated_water$upper,
                                   concentration_separated_water$lower, concentration_separated_water$upper))
  concentration_y_limits <- concentration_y_range + c(-0.10, 0.10) * diff(concentration_y_range)
  concentration_cap_width <- 0.07 * concentration_x_span

  concentration_lines <- data.frame(
    x = seq(concentration_x_limits[1], concentration_x_limits[2], length.out = 200)
  )
  concentration_lines$aggregated <- concentration_aggregated$beta0 + concentration_aggregated$beta1 * concentration_lines$x
  concentration_lines$separated <- concentration_separated$beta0 + concentration_separated$beta1 * concentration_lines$x

  # Left panel: joint concentration line and the joint-fit aqueous confidence interval.
  concentration_aggregated_plot <- ggplot() +
    geom_line(data = concentration_lines, aes(x = x, y = aggregated), color = "#D62728", linewidth = 0.9) +
    geom_errorbar(data = concentration_points, aes(x = x, ymin = lower, ymax = upper),
                  width = concentration_cap_width, color = "#2CA02C", linewidth = 0.7) +
    geom_point(data = concentration_points, aes(x = x, y = y), color = "#2CA02C", shape = 15, size = 3) +
    geom_errorbar(data = concentration_aggregated_water, aes(x = x, ymin = lower, ymax = upper),
                  width = concentration_cap_width, color = "#1F77B4", linewidth = 0.9) +
    geom_point(data = concentration_aggregated_water, aes(x = x, y = y), color = "#1F77B4", shape = 8, size = 4.5) +
    coord_cartesian(xlim = concentration_x_limits, ylim = concentration_y_limits) +
    labs(title = "Aggregated (joint NLS)",
         subtitle = sprintf("Aqueous pKa = %.3f; 95%% CI [%.3f, %.3f]",
                            concentration_aggregated$pka_water,
                            concentration_aggregated_water$lower, concentration_aggregated_water$upper),
         x = expression("ACN volume fraction " * C), y = expression(pK[a])) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          panel.grid.minor = element_blank())

  # Right panel: OLS concentration line and the aqueous confidence interval from SE_c.
  concentration_separated_plot <- ggplot() +
    geom_line(data = concentration_lines, aes(x = x, y = separated),
              color = "#1F77B4", linewidth = 0.9, linetype = "dashed") +
    geom_errorbar(data = concentration_points, aes(x = x, ymin = lower, ymax = upper),
                  width = concentration_cap_width, color = "#2CA02C", linewidth = 0.7) +
    geom_point(data = concentration_points, aes(x = x, y = y), color = "#2CA02C", shape = 15, size = 3) +
    geom_errorbar(data = concentration_separated_water, aes(x = x, ymin = lower, ymax = upper),
                  width = concentration_cap_width, color = "#9467BD", linewidth = 0.9) +
    geom_point(data = concentration_separated_water, aes(x = x, y = y), color = "#9467BD", shape = 8, size = 4.5) +
    coord_cartesian(xlim = concentration_x_limits, ylim = concentration_y_limits) +
    labs(title = "Separated (OLS)",
         subtitle = sprintf("Aqueous pKa = %.3f; 95%% CI [%.3f, %.3f]",
                            concentration_separated$pka_water,
                            concentration_separated_water$lower, concentration_separated_water$upper),
         x = expression("ACN volume fraction " * C), y = expression(pK[a])) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          panel.grid.minor = element_blank())

  concentration_plot <- concentration_aggregated_plot + concentration_separated_plot + plot_annotation(
    title = paste(compound, "- ACN concentration model"),
    caption = paste0(
      "Fits include each observation's supplied f; intervals treat f as fixed.\n",
      "Green squares: independent mixture pKa estimates with 95% intervals. ACN fraction 0.4 means 40% ACN.\n",
      "Stars: aqueous extrapolations at C = 0. Left interval: joint-NLS SE; right interval: SE_c (measurement + lack-of-fit).",
      if (nzchar(data_note)) paste0("\n", data_note) else ""),
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.caption = element_text(size = 9, hjust = 0))
  )
  concentration_file <- file.path(output_dir, paste0("Fig_concentration_", dataset_id, ".png"))
  ggsave(concentration_file, plot = concentration_plot, width = 12, height = 5.7, units = "in", dpi = 200)
  if (interactive()) print(concentration_plot)
  cat("\nConcentration figure saved:", concentration_file, "\n")
}

# =========== 8. Save the combined results after every workbook has finished ===========
summary_file <- file.path(output_dir, "summary_all.csv")
write.csv(all_results, summary_file, row.names = FALSE)
cat("\nDone:", nrow(datasets), "workbooks,", nrow(all_results), "model/method results, and",
    3 * nrow(datasets), "figures.\n")
cat("Summary table:", summary_file, "\nFigures and individual tables:", output_dir, "\n")
