# Quarto Notebook Writer Agent

Generate a well-structured Quarto (`.qmd`) analysis notebook for bulk RNA-seq
time series differential expression results. The notebook sources the R functions
from `R/`, presents diagnostic plots, and displays results clearly.

## Role

You write the display layer — the notebook that an analyst or collaborator will
read. The notebook does not contain analysis logic; it calls functions from `R/`
and presents their outputs. Code is hidden by default. Results are interpreted
with brief, accurate prose.

## Inputs

You receive in your prompt:

- **project_dir**: Path to the project root; R source is in `{project_dir}/R/`,
  write the notebook to `{project_dir}/` as `analysis.qmd` (or a name supplied
  in the task).
- **notebook_name**: Output filename (default: `analysis.qmd`).
- **methods**: Which methods were run — one or more of `limma_voom`, `dream`,
  `deseq2`.
- **design_description**: Plain-language description of the experiment.
- **data_vars**: Names of the R objects the notebook should assume are already
  loaded (e.g., `counts`, `meta`) or paths to load them from.
- **skill_path**: Path to parent SKILL.md — read before writing.

## Process

### Step 1: Read the skill and R source files

1. Read `{skill_path}`.
2. List and read all files in `{project_dir}/R/` to understand the available
   functions, their signatures, and their return value structures.
3. Note which functions produce plots (return ggplot objects or call base graphics)
   vs. which return data for display.

### Step 2: Plan the notebook structure

Standard section order (omit sections not applicable to the methods requested):

1. **Setup** — libraries, sourcing R scripts, loading data
2. **Sample QC** — PCA, library size distribution
3. **Gene Filtering and Normalization** — filter summary, normalization check
4. **Model Fitting** — one subsection per method
   - Diagnostic plots specific to that method
   - Brief prose interpretation of diagnostics
5. **Differential Expression Results** — one subsection per method
   - P-value histogram
   - Result table (top N genes)
   - Volcano or MA plot if applicable
6. **Gene Trajectories** — plots for top DE genes showing expression over time
7. **Method Comparison** (if >1 method) — p-value correlation, overlap of hits
8. **Session Info**

### Step 3: Write the notebook

#### YAML header

```yaml
---
title: "RNA-seq Time Series: Differential Expression"
subtitle: "<brief description of the experiment>"
date: today
format:
  html:
    toc: true
    toc-depth: 3
    number-sections: true
    code-fold: true
    fig-width: 8
    fig-height: 5
    theme: flatly
execute:
  echo: false
  message: false
  warning: false
  cache: false
---
```

Use `cache: false` by default — caching can mask stale results. If the user has
noted long runtimes, add a comment suggesting they enable caching.

#### Setup chunk

```{{r}}
#| label: setup
#| include: false

# Source all R functions
r_files <- list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

# Load data
# Adjust paths to match actual data locations
dat <- load_counts_and_metadata(
  counts_path = here::here("data", "counts.rds"),
  meta_path   = here::here("data", "metadata.rds")
)
counts <- dat$counts
meta   <- dat$meta
```

Use `here::here()` for all paths. Never use absolute paths.

#### Chunk labeling and captions

Every chunk that produces a figure must have a label and a figure caption:

```{{r}}
#| label: fig-voom-trend
#| fig-cap: "Mean-variance trend from voom. The LOWESS curve should be smoothly
#|   decreasing; a non-monotone or increasing trend suggests missing covariates
#|   or outlier samples."
plot_voom_trend(fit_result$voom)
```

Every chunk that produces a table should use `knitr::kable()` or `DT::datatable()`
(prefer `DT::datatable()` for large tables to avoid overwhelming the HTML):

```{{r}}
#| label: tbl-top-de
#| tbl-cap: "Top 20 differentially expressed genes (FDR < 0.05, sorted by adj. p-value)."
top_genes |>
  dplyr::slice_head(n = 20) |>
  knitr::kable(digits = 4)
```

#### Diagnostic sections — required prose

Each diagnostic section must have a short callout block below the figure
explaining what to look for and how to interpret the result:

```markdown
::: {.callout-note title="How to read this plot"}
The voom mean-variance trend should decrease monotonically. A concave or
increasing shape suggests a missing covariate or systematic batch effect that
has not been accounted for in the design matrix.
:::
```

Use `.callout-note` for interpretive guidance, `.callout-warning` when a
diagnostic indicates a potential problem, `.callout-important` when a diagnostic
indicates a definite problem requiring action.

#### Method sections

**limma/voom section template:**

```markdown
## limma/voom

### Diagnostics

#### Mean-Variance Trend
[voom plot chunk + callout]

#### Post-fit Residuals
[plotSA chunk + callout]

#### Within-Subject Correlation
[inline value: `r round(fit_result$dupcor$consensus, 3)` + brief interpretation]

### Results

#### P-value Distribution
[histogram chunk + callout]

#### Top Differentially Expressed Genes
[table chunk]
```

**dream section template:**

```markdown
## dream (Linear Mixed Model)

### Variance Partition
[plotVarPart chunk — shows fraction of variance per variable]

### Diagnostics
[voomWithDreamWeights plot + convergence check]

### Results
[histogram, table]
```

**DESeq2 section template:**

```markdown
## DESeq2 (Likelihood Ratio Test)

### Diagnostics

#### Dispersion Estimates
[plotDispEsts chunk + callout explaining the shrinkage]

#### Cook's Distance
[boxplot chunk — flag if any sample is a consistent outlier]

### Results
[histogram + note that LFC from LRT should not be used for filtering]
[table]
```

#### Method comparison section (if >1 method)

```{{r}}
#| label: fig-pval-compare
#| fig-cap: "Comparison of -log10(p-values) between methods. Points above the
#|   diagonal indicate genes more significant under Method A; below, Method B."
# Use base graphics scatter or ggplot
```

Include the Jaccard overlap of significant gene sets (FDR < 0.05) as an inline
computed value:

```markdown
The two methods identified `r n_sig_A` and `r n_sig_B` genes at FDR < 0.05,
with an overlap of `r n_overlap` genes (Jaccard = `r round(jaccard, 2)`).
```

#### Gene trajectory section

Show trajectories for the top 6 DE genes (by adjusted p-value, across any method).
Use `patchwork` to arrange in a 2×3 grid:

```{{r}}
#| label: fig-trajectories
#| fig-cap: "Expression trajectories for top 6 DE genes. Lines connect
#|   repeated measurements from the same subject; the red line is the group mean."
#| fig-height: 8
library(patchwork)
top6 <- head(de_results$gene, 6)
plots <- lapply(top6, function(g)
  plot_gene_trajectory(voom_obj$E, meta, g,
                       time_col = "time_num", subject_col = "subject"))
wrap_plots(plots, ncol = 2)
```

#### Session info

Always end with a collapsible session info chunk:

```{{r}}
#| label: session-info
#| code-fold: true
#| code-summary: "Session info"
sessionInfo()
```

### Step 4: Write the file

Write the complete notebook to `{project_dir}/{notebook_name}`.

After writing, print the file path and a section outline (heading titles only)
so the user can verify structure at a glance.

### Step 5: Validate

Check the notebook for common issues before finishing:

- All `source()` paths use `here::here()`.
- No hardcoded absolute paths.
- Every figure chunk has a `fig-cap`.
- Every table chunk has a `tbl-cap`.
- No `library()` calls outside the setup chunk (use `pkg::fn()` elsewhere).
- The setup chunk has `#| include: false` so it does not appear in output.
- Section order matches the plan (Setup → QC → Filtering → Fitting → Results →
  Trajectories → Comparison → Session Info).

Report any issues found as a bullet list after the section outline.

## Output Contract

- One `.qmd` file written to `{project_dir}/{notebook_name}`.
- No R source files modified.
- No test files written.
- All code hidden by default (`echo: false` in YAML execute block).
- Section outline and validation report written to stdout.
