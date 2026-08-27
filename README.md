# Aqueous pKa of PFAS compounds from 19F NMR titration

Estimation of the aqueous acid dissociation constant of five PFAS compounds from
19F NMR chemical-shift titrations measured in three ACN/water mixtures, under two
solvent models and by two estimation routes.

Running `Rscript run_all.R` regenerates every figure, table and number in
`output/report.pdf` from the spreadsheets in `data/`. `Rscript verify.R` re-checks the
algebraic identities the analysis relies on.

---

## 1. Data

Five compounds (5:3 FTCA, 7:3 FTCA, 8:2 FTCA, PFBA, PFOA) titrated at 40%, 50% and
60% (v/v) ACN; one workbook per dataset, one sheet per mixture. Measured by
Dr. Isabelle Logan (NMR Facility, Oregon State University).

**The spreadsheets are not distributed with this repository** — they are unpublished
experimental measurements. Obtain them from the project's shared Box folder and place
them in `data/`; `data/README.md` lists the expected file names and sheet structure, and
`read_titration(..., data_dir = ...)` reads them from any other location.

For 5:3 FTCA and 7:3 FTCA more than one spreadsheet exists — a second CF2 resonance for
both, and an earlier measurement session for 7:3. All of them are analyzed; the ten
`DATASETS` entries in `R/data.R` list which file each result comes from.

Columns are selected **by name**, never by position: the workbooks carry different extra
columns, the pH header appears as `pH` or `Adjusted pH` in either capitalisation, and the
shift column is `Chemical Shift (ppm)` in most sheets but `deltai_calc` in some (the
5:3 other-CF2 file); whichever of the two is present is the resonance that file reports. **All fits use the `Adjusted pH` column** (the probe reading corrected by the
experimentally measured offset). The raw `pH` column can be read with
`read_titration(..., ph_col = "raw")` for comparison.

Solvent constants are **fixed in `R/constants.R` and are not read from the spreadsheets**
(the `water` column in the files is outdated):

| Mixture | dielectric constant | water activity | ACN volume fraction |
|---|---|---|---|
| 40% ACN | 63.50 | 0.913 | 0.40 |
| 50% ACN | 58.71 | 0.908 | 0.50 |
| 60% ACN | 53.91 | 0.901 | 0.60 |
| pure water | 78.33 | 1 | 0 |

Source: collaborator email of 2026-07-15 (literature dielectric constants and water
activities).

## 2. Model

Each titration follows the Henderson–Hasselbalch relation between chemical shift and pH,

```
y_ki = ( dA_k + dHA_k * 10^(pKa_k - pH_ki) ) / ( 1 + 10^(pKa_k - pH_ki) ) + e_ki
```

with `dA_k`, `dHA_k` the two limiting shifts of mixture `k` and `pKa_k` its effective
pKa. The two solvent models share one algebraic form,

```
pKa_k = theta' d_k + offset_k ,        pKa_aqueous = theta' g
```

and differ only in the covariate, the known offset and the read-off point:

| | covariate `d_k` | `offset_k` | read-off `g` |
|---|---|---|---|
| Yasuda–Shedlovsky | `(1, 1/eps_k)` | `-log10(a_k)` | `(1, 1/78.33)` |
| ACN concentration | `(1, C_k)` | `0` | `(1, 0)` |

They are therefore one implementation with two configurations, built by
`model_spec("ys")` and `model_spec("conc")`. Under the concentration model the second
component of `g` is zero, so the aqueous pKa is the intercept itself.

## 3. Estimation

**Route A — joint (aggregated) fit.** All titration points are fitted at once by
non-linear least squares with the solvent model imposed inside the mean function; the
parameters are `theta` plus the three pairs `(dA_k, dHA_k)`. The fit is carried out in
the centered covariate `x_k - x0`, so the aqueous value is itself a parameter; in the
natural coordinates the two coefficients correlate at about −0.99 and the optimizer's
stopping point depends on how the covariate is written. Centering removes that ridge and
attained the lowest residual sum of squares of the parameterizations tried.
`coef_natural()` maps the result back to the reported coefficients.

The covariance of `theta` is computed **analytically** as

```
V = sigma^2 * M^(-1) ,   M = sum_k w_k d_k d_k' ,   w_k = || (I - P_k) u_k ||^2
```

where `u_k` collects the sensitivities of the fitted curve to `pKa_k` and `P_k` projects
onto the two limiting-shift directions, i.e. `w_k` is the information mixture `k` retains
about `theta` after its limiting shifts are profiled out. The optimizer's
finite-difference `vcov()` is not used for reporting; `verify.R` compares the two.

**Route B — two-step.** Each mixture is fitted separately (three parameters), giving
`pKa_k` with variance `nu_k` — the green squares and their error bars in the figures.
Ordinary least squares of those three values on `d_k` then gives `theta~`, with

```
V_c = (Z'Z)^-1 Z' [ diag(nu) + tau2 * I ] Z (Z'Z)^-1
tau2 = max( 0 , SS_res - sum_k (1 - h_k) nu_k )
```

`tau2` is the lack-of-fit component (a DerSimonian–Laird moment estimator): scatter of
the three points about the fitted line beyond what their own measurement variances
explain. It is zero whenever the points are consistent with the line, and `V_c` then
reduces to pure propagation of the `nu_k`.

**Standard errors.** Every standard error in the report is the same quadratic form
`Var(aqueous) = g' V g` evaluated at a different `V`. The report uses `SE_joint` for
Route A and `SE_c` for Route B; the CSV additionally carries three reference columns:

| column | `V` used | meaning |
|---|---|---|
| `SE_joint` | analytic covariance of Route A | reported for the joint route |
| `SE_c` | `V_c` above | reported for the two-step route |
| `SE_a` | `V_c` with `tau2` forced to 0 | measurement propagation only, reference |
| `SE_b` | equal-variance OLS mean-response SE | the earlier convention, reference |
| `SE_joint_numeric` | the optimizer's own SE for the aqueous parameter | independent cross-check of `SE_joint`, compared in `verify.R` |

**Lack-of-fit statistic.** `Q = SS_res / sum_k (1 - h_k) nu_k` compares the observed
scatter of the three points with the scatter their measurement variances alone would
produce; under exact adherence to the solvent model `Q ~ chi-square with 1 df`, and
`tau2 > 0` exactly when `Q > 1`.

## 4. Layout

```
run_all.R          single entry point: all datasets x both models -> output/
verify.R           algebraic self-checks (see section 5)
report.tex         source of the report; run_all.R copies it into output/ and compiles
R/constants.R      solvent constants; model_spec(); the quadratic form g'Vg
R/data.R           dataset registry; column-name-based spreadsheet reader
R/fit.R            Henderson-Hasselbalch; per-mixture fit; joint fit; analytic covariance
R/inference.R      analyze(): aqueous estimates, the four SEs, tau2, Q
R/figures.R        titration figure; two-panel extrapolation figure
R/tables.R         LaTeX table fragments
data/              the spreadsheets (not tracked; see data/README.md)
output/            everything generated (not tracked)
```

## 5. Outputs and where they appear in the report

| output file | contents | used in |
|---|---|---|
| `summary_all.csv` | one row per dataset x model, all fitted quantities at full precision | source of every number below |
| `Fig_raw_<key>.pdf` | titration curves vs Adjusted pH | first page of each dataset |
| `Fig_ys_<key>.pdf` | Yasuda–Shedlovsky extrapolation, two panels | first page of each dataset |
| `Fig_conc_<key>.pdf` | concentration extrapolation, two panels | second page of each dataset |
| `tab_ys_<key>.tex`, `tab_cc_<key>.tex` | per-model coefficient and diagnostic table | under the corresponding figure |
| `tab_cmp_<key>.tex` | the two models side by side for that dataset | second page of each dataset |
| `tab_overview_part<N>.tex` | aqueous pKa for all datasets of a part, both models | opening page of each part |
| `tab_diag_part<N>.tex` | `tau2`, `Q`, `p` for all datasets, both models | opening page of each part |
| `report.pdf` | the assembled report | the deliverable |

`<key>` is the dataset key of `R/data.R` (`5_3_FTCA`, `PFOA`, `v7_3_old`, …); keys
beginning with `v` are the Part 2 spreadsheets.

## 6. Checks performed by `verify.R`

0. the model configuration (see below);
1. the aqueous estimate equals `g' theta` for both routes;
2. `SE_c` equals `sqrt(g' V_c g)` rebuilt from the reported parameter standard errors and
   correlation — i.e. a reader can reproduce it from the printed table;
3. under the concentration model the aqueous standard error equals the intercept's
   standard error, for both routes;
4. `tau2 > 0` exactly when `Q > 1`;
5. the two models share identical stage-1 estimates (those are fitted per mixture and
   cannot depend on the solvent model);
6. the analytic joint covariance agrees with the optimizer's own standard error for the
   aqueous parameter, to 5% relative (currently ~1e-4). Because the joint fit is run in
   coordinates centered at the read-off point, that parameter *is* the aqueous value, so
   its standard error is an independent estimate of `g' V g`: this check is what would
   catch a transposed covariance matrix or a mis-contracted read-off vector;
7. agreement with `summary_all_from_code.csv` of the report folder, if reachable.

Check 0 additionally asserts the model configuration itself — each spec's covariate,
read-off point and offset, and the design identities `sum(ell) = 1`,
`sum(ell * x) = x0` — so that a mis-specified model cannot pass by propagating
consistently through the identities.

Any failure stops the script.

## 7. Requirements

R with `readxl`, `minpack.lm`, `ggplot2`, `patchwork`; `pdflatex` for the report step
(skipped with a message if absent).

```
Rscript run_all.R      # figures, tables, CSV, report.pdf
Rscript verify.R       # self-checks
```
