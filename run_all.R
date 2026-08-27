# ==============================================================
# run_all.R -- single entry point.
#
#   Rscript run_all.R
#
# Runs every dataset in data.R under both solvent models, writes all
# figures, LaTeX table fragments and summary CSVs into output/, and
# then compiles the report if pdflatex is available.  Everything in
# the report can be regenerated from this one command.
# ==============================================================
for (f in c("constants", "data", "fit", "inference", "figures", "tables"))
  source(file.path("R", paste0(f, ".R")))

OUT <- "output"
dir.create(OUT, showWarnings = FALSE)
specs <- list(ys = model_spec("ys"), conc = model_spec("conc"))

rows <- list(); body <- list(`1` = character(0), `2` = character(0))

for (ds in DATASETS) {
  data <- read_titration(ds$file, ph_col = "adjusted")
  res  <- lapply(specs, function(sp) analyze(data, sp, ds$label))
  rows[[ds$key]] <- do.call(rbind, lapply(res, function(r) r$row))

  ggsave(file.path(OUT, sprintf("Fig_raw_%s.pdf", ds$key)),
         fig_titration(data, ds$label), width = 6, height = 4.4)
  for (m in names(specs))
    ggsave(file.path(OUT, sprintf("Fig_%s_%s.pdf", m, ds$key)),
           fig_extrapolation(res[[m]],
             sprintf("%s: %s extrapolation", ds$label, specs[[m]]$name)),
           width = 13, height = 5.5)

  table_dataset(res$ys$row,   file.path(OUT, sprintf("tab_ys_%s.tex", ds$key)))
  table_dataset(res$conc$row, file.path(OUT, sprintf("tab_cc_%s.tex", ds$key)))
  table_models(res$ys$row, res$conc$row,
               file.path(OUT, sprintf("tab_cmp_%s.tex", ds$key)))

  p <- as.character(ds$part)
  body[[p]] <- c(body[[p]],
    "\\clearpage",
    sprintf("\\subsection*{%s}", ds$label),
    sprintf("Source spreadsheet: \\texttt{%s}\\quad (\\emph{Adjusted pH} column)",
            .tex_esc(ds$file)),
    "", "\\textbf{Titration data}",
    sprintf("\\begin{center}\\includegraphics[width=0.52\\textwidth]{Fig_raw_%s.pdf}\\end{center}", ds$key),
    "\\textbf{Model 1: Yasuda--Shedlovsky}",
    sprintf("\\begin{center}\\includegraphics[width=\\textwidth]{Fig_ys_%s.pdf}\\end{center}", ds$key),
    sprintf("\\input{tab_ys_%s.tex}", ds$key),
    "", "\\clearpage",
    sprintf("\\subsection*{%s \\ (continued)}", ds$label),
    "\\textbf{Model 2: ACN concentration}",
    sprintf("\\begin{center}\\includegraphics[width=\\textwidth]{Fig_conc_%s.pdf}\\end{center}", ds$key),
    sprintf("\\input{tab_cc_%s.tex}", ds$key),
    "", "\\textbf{The two models side by side}",
    sprintf("\\input{tab_cmp_%s.tex}", ds$key), "")

  cat(sprintf("%-22s YS %7.3f (%.3f) | Conc %7.3f (%.3f)\n", ds$label,
              res$ys$aq_joint, res$ys$se_joint,
              res$conc$aq_joint, res$conc$se_joint))
}

all_rows <- do.call(rbind, rows)
write.csv(all_rows, file.path(OUT, "summary_all.csv"), row.names = FALSE)

for (p in c(1, 2)) {
  keys  <- sapply(Filter(function(d) d$part == p, DATASETS), `[[`, "key")
  files <- sapply(Filter(function(d) d$part == p, DATASETS), `[[`, "file")
  ys <- all_rows[all_rows$model == "ys"   & rownames(all_rows) %in%
                   paste0(keys, ".ys"), ]
  cc <- all_rows[all_rows$model == "conc" & rownames(all_rows) %in%
                   paste0(keys, ".conc"), ]
  table_overview(ys, cc, files, file.path(OUT, sprintf("tab_overview_part%d.tex", p)))
  table_diagnostics(ys, cc, file.path(OUT, sprintf("tab_diag_part%d.tex", p)))
  writeLines(body[[as.character(p)]], file.path(OUT, sprintf("body_part%d.tex", p)))
}

file.copy("report.tex", file.path(OUT, "report.tex"), overwrite = TRUE)
latex <- Sys.which("pdflatex")
if (nzchar(latex) || file.exists("/Library/TeX/texbin/pdflatex")) {
  bin <- if (file.exists("/Library/TeX/texbin/pdflatex"))
    "/Library/TeX/texbin/pdflatex" else latex
  wd <- setwd(OUT)
  for (i in 1:2) system2(bin, c("-interaction=nonstopmode", "report.tex"),
                         stdout = FALSE, stderr = FALSE)
  setwd(wd)
  cat("\nreport written: output/report.pdf\n")
} else {
  cat("\npdflatex not found; figures/tables written, PDF not compiled.\n")
}
