# Fast cigarette-tax runs

A fast reimplementation of the cigarette-tax analysis in `R/`. It produces the same outputs from the same inputs:

- **`run_tax_grid.R`** writes the TCP tool files for the price-per-pack × added-tax grid for every state and the US. On a 12-core M3 Pro this takes about 25 s, or 15 s when every scenario is already cached. The original pipeline needed `R/03_main_analysis.R`, then `processing.ipynb`, then `update_us_results.py`, an estimated 19 hours on the same machine.
- **`run_tax_figures.R`** runs the no-policy baseline and the +$1 to +$5 tax runs, then draws the tax figures from `R/04_visualization_policy_specific.R` for the US and every state, as SVGs with an `index.html` page showing them all. Each state's prevalence figures also show the US in grey. This takes about 20 s.
  - The US is all states together, each with its own tax increase. Deaths and life-years are summed over states, and prevalence is weighted by each state's census population.

Only the tax policy is covered. Smoke-free air laws, tobacco control spending, warnings and the other scripts in `R/` are unchanged and still run as before.

## Setup

1. Copy `config.example.yml` to `config.yml` and set the input and output paths. `config.yml` is git-ignored because its paths are specific to your machine.
   - Relative paths in the file are resolved from the file's own folder.
   - Any setting can be overridden for a single run with a flag, writing `_` as `-`, for example `--final-dir=/tmp/out --states=AK,CA --cache=off`.
   - `--config=<file>` uses a different config file.
2. R packages install themselves: each script installs any it needs from CRAN the first time it runs (`lib/packages.R`). The grid uses `Rcpp`, `readr`, `vroom`, `usmap` and `yaml`; the figures also use `dplyr`, `ggplot2` and `scales`.
   - The compiled kernel needs a C++ compiler (Xcode Command Line Tools on macOS).
   - Without a compiler, the scripts fall back to pure R. The results are the same; only the speed changes.

## Run

Run from the repository root:

```bash
Rscript R-fast/run_tax_grid.R      # TCP tool files -> grid.final_dir
Rscript R-fast/run_tax_figures.R   # .rds lists, SVG figures and index.html -> figures.out_dir
```

## Files

| File | Purpose |
|---|---|
| `run_tax_grid.R` | TCP tool files for the scenario grid (states and US) |
| `run_tax_figures.R` | `baseline_results.rds`, `tax_out_list.rds`, `figs/tax/*.svg` (and/or `.jpg`) and `index.html` |
| `lib/tcp_fast.R` | The model: same functions as `R/02_model_functions.R`, fast |
| `lib/cohort_kernel.cpp` | Compiled cohort model and smoking-attributable-death sum |
| `lib/us_results.R` | US results: population-weighted average of the state results |
| `lib/config.R` | Reads `config.yml` and command-line overrides |
| `lib/gallery.R` | `index.html`: every figure on one page, by state |
| `lib/packages.R` | Installs missing R packages from CRAN on first run |
| `config.example.yml` | Template listing every setting |
| `tools/check_equivalence.R` | Test: fast functions give `identical()` results to `R/02_model_functions.R` |
| `tools/compare_trees.R` | Test: compares two output trees file by file and cell by cell |

The tax figure functions live in `R/tax_plots.R`. Both `R/04_visualization_policy_specific.R` and `run_tax_figures.R` load that file, so the figures are defined in one place.

## How it differs from the original

The model arithmetic is unchanged. `lib/tcp_fast.R` performs every operation in the same order as `R/02_model_functions.R`, and `tools/check_equivalence.R` checks that every object it returns is `identical()` to the original's. The speed comes from doing less repeated work.

| Original (`R/`) | Fast (`R-fast/`) |
|---|---|
| Runs every grid tag: 51 states × 90 tags × 7 policy years = 32,130 `runstates()` calls | Runs one scenario per distinct total tax increase, `max(0, target − state price) + added tax`: 843 scenarios. Tags with the same total produce identical files. |
| Each `runstates()` call reloads the state's inputs and reruns the no-policy baseline | The inputs, mortality tables and baseline are computed once per state |
| Loops over cohort, age and year in R: 5–8 s per `runstates()` call | Vectorised R, or a compiled C++ kernel with the same arithmetic: about 6 ms per call |
| `generate_TCPoutput()` rebuilds every table for each scenario | `tcp_tables()` computes the row layout and baseline columns once per state and policy year |
| Writes per-scenario CSVs. `processing.ipynb` (pandas) then reshapes them into the TCP tool tree, and `update_us_results.py` builds the US results. | Writes the TCP tool tree and the US files directly, in a single step |
| A full copy of every duplicate file | Duplicates are symlinks to one file (`link: symlink`). `aws s3 cp` / `sync` uploads them as real files. |
| Recomputes everything on each run | The scenario cache reuses results when the inputs, code and settings haven't changed |

Output compared with the original pipeline:

- **Compared with the earlier version of this code**, which emulated pandas to reproduce `processing.ipynb` byte for byte: every state file is byte-identical. Of the 98,280 files, only 7 US deaths and lyg files differ. Each has one near-zero cell that is exactly 0 here and pandas float noise (~1e-29) in that version.
- **`tcp-data` 2.0.4:** both this code and the original R code differ from the published 2.0.4 only in near-zero cells of deaths and lyg, such as `-3.6e-15` against `0`. All differences are at most 1e-8 in absolute size, or 2.3e-12 relative to the column's largest value, and every results file is identical.
- **`run_tax_figures.R` lists:** they keep only `df_mort.outputs` and `df_CSprevs.by.state` (plus `initprice` and `tax_raise`), which is everything the figure scripts read. The original also stored `l_prev_out`, about 50 MB per scenario.

## Checking changes

```bash
Rscript R-fast/tools/check_equivalence.R --states=AK,CA --policy-years=2026,2040   # against HEAD's R/02
Rscript R-fast/tools/compare_trees.R --ref=<previous output> --new=<new output>
```
