# Data directory

The titration spreadsheets are **not** distributed with this repository: they are the
collaborators' unpublished experimental measurements (19F NMR chemical-shift titrations,
NMR Facility, Oregon State University). Obtain them from the project's shared Box folder
and place them here.

Expected files (as shared on 2026-08-26):

```
5_3_FTCA_pKa_Dielectric.xlsx
5_3_FTCA_pKa_other_CF2.xlsx
7_3_FTCA_pKa_Dielectric_Newest_Data.xlsx
7_3_FTCA_pKa_Dielectric_Old_Data.xlsx
7_3_FTCA_pKa_other_CF2.xlsx
8_2_FTCA_pKa_Dielectric.xlsx
PFBA_pKa_Dielectric.xlsx
PFOA_pKa_Dielectric.xlsx
```

Each workbook has one sheet per mixture, named `40% ACN`, `50% ACN`, `60% ACN`, and
carries at least a pH column (`pH` and `Adjusted pH`) and a chemical-shift column, named
either `Chemical Shift (ppm)` or `deltai_calc` (both appear among the shared files). Columns are located by
name, so additional columns and header capitalisation do not matter.

To read the spreadsheets from somewhere else instead, pass `data_dir` to
`read_titration()` (see `R/data.R`).
