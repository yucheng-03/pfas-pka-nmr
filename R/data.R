# ==============================================================
# data.R -- reading the titration spreadsheets.
#
# One workbook per dataset, one sheet per mixture ("40% ACN",
# "50% ACN", "60% ACN").  Columns are selected BY NAME, not by
# position, because the workbooks carry different extra columns and
# spell the headers differently:
#   pH column     : "pH" (probe reading) or "Adjusted pH" (corrected)
#                   -- the header appears in both capitalisations
#   shift column  : "Chemical Shift (ppm)" in most files,
#                   "deltai_calc" in the other-CF2 files
# The analysis uses the "Adjusted pH" column throughout; the raw
# "pH" column is readable via ph_col = "raw" for comparison only.
# ==============================================================
suppressPackageStartupMessages(library(readxl))

DATASETS <- list(
  # the five compounds (Part 1 of the report)
  list(key = "5_3_FTCA", label = "5:3 FTCA", part = 1,
       file = "5_3_FTCA_pKa_Dielectric.xlsx"),
  list(key = "7_3_FTCA", label = "7:3 FTCA", part = 1,
       file = "7_3_FTCA_pKa_Dielectric_Newest_Data.xlsx"),
  list(key = "8_2_FTCA", label = "8:2 FTCA", part = 1,
       file = "8_2_FTCA_pKa_Dielectric.xlsx"),
  list(key = "PFBA",     label = "PFBA",     part = 1,
       file = "PFBA_pKa_Dielectric.xlsx"),
  list(key = "PFOA",     label = "PFOA",     part = 1,
       file = "PFOA_pKa_Dielectric.xlsx"),
  # every spreadsheet shared for 5:3 and 7:3 (Part 2 of the report)
  list(key = "v5_3_main",     label = "5:3 FTCA (main CF2)",  part = 2,
       file = "5_3_FTCA_pKa_Dielectric.xlsx"),
  list(key = "v5_3_otherCF2", label = "5:3 FTCA (other CF2)", part = 2,
       file = "5_3_FTCA_pKa_other_CF2.xlsx"),
  list(key = "v7_3_newest",   label = "7:3 FTCA (newest)",    part = 2,
       file = "7_3_FTCA_pKa_Dielectric_Newest_Data.xlsx"),
  list(key = "v7_3_old",      label = "7:3 FTCA (old)",       part = 2,
       file = "7_3_FTCA_pKa_Dielectric_Old_Data.xlsx"),
  list(key = "v7_3_otherCF2", label = "7:3 FTCA (other CF2)", part = 2,
       file = "7_3_FTCA_pKa_other_CF2.xlsx"))

read_titration <- function(file, data_dir = "data",
                           ph_col = c("adjusted", "raw")) {
  ph_col <- match.arg(ph_col)
  want_ph <- switch(ph_col, adjusted = "adjusted ph", raw = "ph")
  do.call(rbind, lapply(seq_along(CONDITIONS), function(k) {
    df <- suppressMessages(read_excel(file.path(data_dir, file),
                                      sheet = CONDITIONS[k]))
    names(df) <- tolower(trimws(names(df)))
    if (!want_ph %in% names(df))
      stop(sprintf("column '%s' not found in %s / %s", want_ph, file, CONDITIONS[k]))
    shift <- intersect(c("chemical shift (ppm)", "deltai_calc"), names(df))
    if (length(shift) == 0)
      stop(sprintf("no chemical-shift column in %s / %s", file, CONDITIONS[k]))
    pH <- suppressWarnings(as.numeric(df[[want_ph]]))
    y  <- suppressWarnings(as.numeric(df[[shift[1]]]))
    ok <- !is.na(pH) & !is.na(y)
    data.frame(pH = pH[ok], ChemShift = y[ok], grp = k,
               Condition = factor(CONDITIONS[k], levels = CONDITIONS))
  }))
}
