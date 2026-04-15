---
name: rnaseq-timeseries-de
description: >
  Differential expression analysis for bulk RNA-seq time series and repeated-measures
  experiments. Use this skill whenever the user has longitudinal RNA-seq data, multiple
  timepoints per subject, pre/post treatment designs, or any experiment where the same
  biological unit is measured more than once. Covers three complementary methods:
  limma/voom with duplicateCorrelation or spline contrasts, dream (variancePartition
  linear mixed models), and DESeq2 with the likelihood ratio test. Emphasizes
  statistical model specification, model diagnostics, and assumption checking throughout.
  Trigger when the user mentions: time course, time series, repeated measures, longitudinal
  RNA-seq, paired samples across timepoints, within-subject correlation, mixed models for
  RNA-seq, dream/variancePartition DE, or DESeq2 LRT for timecourse. Always use this
  skill even if the user only asks about one method, to provide full context on tradeoffs.
---

# Bulk RNA-seq: Time Series and Repeated-Measures Differential Expression

## Overview and Method Selection

Time series and repeated-measures RNA-seq experiments present a fundamental challenge: observations from the same subject are not independent. Ignoring this correlation inflates false positives; naive blocking (fixed subject effects) can be severely underpowered.

**Before choosing a method, consult `references/oh-li-2021-review.md`**, which synthesizes a comprehensive 2021 review (Oh & Li, *Genes* 12:352, PMC7997275) covering 16+ dynamic methods, normalization comparisons, and batch correction recommendations for time course data. What follows is a practical subset of that guidance focused on the three core methods implemented here, plus pointers to alternatives.

### Core methods in this skill

| Method | Handles repeated measures? | Spline/smooth time? | Multiple random effects? | Count-based model? |
|---|---|---|---|---|
| limma/voom + duplicateCorrelation | Single random effect (compound symmetry) | Yes (splines in design) | No | No (voom weights) |
| dream (variancePartition) | Yes, gene-specific random effects | Yes | Yes | No (voom weights) |
| DESeq2 LRT | No native random effects | Yes (splines in design) | No | Yes (NB GLM) |

**General guidance:**
- For experiments with a simple single blocking factor (subject/individual) and moderate N: limma/voom + `duplicateCorrelation` is a fast, robust starting point.
- For complex designs, multiple random effects, or when the donor contribution varies substantially across genes: use **dream**.
- When you need a native count-based (negative binomial) model, or are testing omnibus hypotheses across many timepoints: use **DESeq2 LRT**.
- For microarray data (not RNA-seq), see the `timecourse` package section below.

### Alternative methods (from Oh & Li 2021)

Consult `references/oh-li-2021-review.md` Section 2 for full descriptions.

| Scenario | Method | Key feature |
|---|---|---|
| Repeated measures, AR correlation between adjacent timepoints | **rmRNAseq** | Continuous autoregressive correlation; voom-based |
| Two-condition, impulse/on-off dynamics, ≥6 timepoints | **ImpulseDE2** | NB model; no need to pre-specify k |
| Two-condition, polynomial trajectory, exploratory | **Next maSigPro** | Stepwise model selection; R² per gene |
| Two-condition, splines + downstream network inference | **splineTimeR** | Unified DE + gene association network |
| Single-series, classify state transitions | **EBSeq-HMM** | HMM; classifies each transition as DE/EE |
| Single-series, identify breakpoints and slopes | **Trendy** | Segmented regression + BIC; no replicates required |
| Circadian / periodic, repeated measures | **LimoRhyde** | Cosinor regression; integrates with limma/voom |

**Code organization conventions:**
- Analysis functions go in `R/` scripts (one function per file or logically grouped).
- Quarto notebooks (`.qmd`) are used for display of results; use `.Rmd` only if explicitly requested.
- In notebooks, set `echo: false` and `message: false` and `warning: false` in YAML or chunk options by default, and show only key outputs.

---

## Agent Delegation

This skill delegates to three specialized agents. Invoke them in order — each
depends on the prior stage's outputs.

| Agent | File | When to invoke |
|---|---|---|
| Code Writer | `agents/code-writer.md` | User asks to write or update R source functions |
| Tester | `agents/tester.md` | After code-writer finishes, or when user asks to test existing functions |
| Notebook Writer | `agents/notebook-writer.md` | After tester passes, or when user asks to generate/update the Quarto notebook |

**Typical full workflow:**

1. Read this SKILL.md to understand the analysis context.
2. Invoke **code-writer** with the method(s) and design description → produces `R/` files.
3. Invoke **tester** against those files → produces `tests/` files and a pass/fail report.
4. Once tests pass, invoke **notebook-writer** → produces `analysis.qmd`.

**Partial workflows** are common and fine:
- User already has R source → skip to tester or notebook-writer.
- User only wants to add a diagnostic function → invoke code-writer scoped to
  `R/diagnostics.R` only, then re-run tester.
- User only wants to regenerate the notebook → invoke notebook-writer directly.

**Shared context** passed to every agent:
- `skill_path`: path to this SKILL.md (agents must read it before acting)
- `project_dir`: root of the analysis project
- `design_description`: plain-language description of the experimental design

---

## 1. Data Inputs and Preprocessing

Starting point is a features × samples count matrix and a metadata data frame. The metadata must contain a `sample` column matching count matrix column names, plus columns for subject/individual ID, time variable, and any covariates.

```r
# R/load_data.R
load_counts_and_metadata <- function(counts_path, meta_path) {
  counts <- readRDS(counts_path)          # genes x samples integer matrix
  meta   <- readRDS(meta_path)            # data.frame with column 'sample'
  stopifnot(all(meta$sample %in% colnames(counts)))
  counts <- counts[, meta$sample]         # ensure column order matches
  list(counts = counts, meta = meta)
}
```

### Filtering low-count genes

Before any modelling, remove genes with consistently low counts. A common threshold: keep genes with CPM > 1 in at least as many samples as the smallest group.

```r
# R/filter_genes.R
filter_low_counts <- function(counts, meta, group_col, min_cpm = 1) {
  library(edgeR)
  dge <- DGEList(counts = counts)
  min_samples <- min(table(meta[[group_col]]))
  keep <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge[keep, ]
}
```

### Normalization

Use TMM normalization (edgeR) before voom-based methods. DESeq2 handles normalization internally via size factors.

```r
dge <- calcNormFactors(dge, method = "TMM")
```

---

## 2. limma/voom for Repeated Measures

### Statistical model

For a design with subjects measured at multiple timepoints, the linear model for gene $g$ in subject $i$ at time $t$ is:

$$y_{git} = \mathbf{x}_{it}^{\top} \boldsymbol{\beta}_g + \epsilon_{git}$$

where $\mathbf{x}_{it}$ is the design vector (timepoint indicators, treatment, covariates) and $\epsilon_{git}$ is the residual error. After voom transformation, $\tilde{y}_{git} = \log_2(\text{count}_{git} + 0.5) - \log_2(\text{lib.size}_i / 10^6)$ with precision weights $w_{git}$ derived from the mean-variance trend.

For repeated measures, the within-subject correlation is modelled as:

$$\operatorname{Cov}(\epsilon_{git}, \epsilon_{git'}) = \sigma^2_g \rho \quad (t \neq t')$$

where $\rho$ is a single consensus correlation estimated by `duplicateCorrelation`. This assumes **compound symmetry** (equal correlation between any pair of timepoints for a given subject) — an important assumption to keep in mind.

### Workflow: `duplicateCorrelation` (simple repeated measures)

```r
# R/limma_voom_dupcor.R
run_voom_dupcor <- function(dge, meta, design_formula, block_col) {
  library(limma)
  library(edgeR)

  design <- model.matrix(design_formula, data = meta)
  block  <- meta[[block_col]]

  # First voom pass: get preliminary weights
  v0     <- voom(dge, design, plot = FALSE)
  # Estimate within-block (within-subject) correlation
  dupcor <- duplicateCorrelation(v0, design, block = block)
  message("Consensus correlation: ", round(dupcor$consensus, 3))

  # Second voom pass: incorporate correlation into precision weights
  v      <- voom(dge, design, block = block,
                 correlation = dupcor$consensus, plot = TRUE)
  # Fit linear model
  fit    <- lmFit(v, design, block = block,
                  correlation = dupcor$consensus)
  list(voom = v, fit = fit, dupcor = dupcor, design = design)
}
```

**Why two voom passes?** The first pass gives preliminary weights needed to estimate the within-subject correlation; the second pass recomputes weights accounting for that correlation, yielding more accurate precision weights.

### Contrasts and testing

For pairwise or complex comparisons, construct a contrast matrix:

```r
# R/limma_contrasts.R
make_contrasts_and_test <- function(fit, contrast_matrix) {
  library(limma)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, robust = TRUE)
  fit2
}

# Example: time points as ordered factor, test T2 vs T0
# contr <- makeContrasts(T2vsT0 = timeT2 - timeT0, levels = design)
```

### Time series: splines and F-tests (Chapter 15/18, limma User's Guide)

For a continuous or many-level time variable, model time using natural cubic splines. This avoids estimating one parameter per timepoint and allows smooth estimation of the time trajectory:

$$\mathbf{x}_{it} = [1, \, b_1(t), \, b_2(t), \ldots, b_k(t), \, \ldots]$$

where $b_j(t)$ are spline basis functions. An omnibus F-test then asks: "Does expression change over time at all?" by testing all spline coefficients simultaneously.

```r
# R/limma_spline_timecourse.R
run_spline_timecourse <- function(dge, meta, time_col, subject_col,
                                  covariate_cols = NULL, df = 4) {
  library(limma); library(splines)

  # Natural cubic spline basis for time
  time_spline <- ns(meta[[time_col]], df = df)
  colnames(time_spline) <- paste0("ns", seq_len(df))
  meta_aug <- cbind(meta, time_spline)

  # Build design: intercept + spline terms + optional covariates
  fixed_terms <- paste(c(colnames(time_spline), covariate_cols), collapse = " + ")
  design <- model.matrix(as.formula(paste("~", fixed_terms)), data = meta_aug)

  # voom + duplicateCorrelation
  block  <- meta[[subject_col]]
  v0     <- voom(dge, design, plot = FALSE)
  dupcor <- duplicateCorrelation(v0, design, block = block)
  v      <- voom(dge, design, block = block,
                 correlation = dupcor$consensus, plot = FALSE)
  fit    <- lmFit(v, design, block = block,
                  correlation = dupcor$consensus)
  fit    <- eBayes(fit, robust = TRUE)

  # Omnibus F-test across all spline coefficients
  spline_cols <- grep("^ns", colnames(design), value = TRUE)
  fit_F <- topTable(fit, coef = spline_cols, number = Inf, sort.by = "F")
  list(fit = fit, voom = v, dupcor = dupcor, spline_cols = spline_cols,
       omnibus = fit_F, design = design)
}
```

### `timecourse` package (microarray / balanced time series)

For **microarray** data or pre-normalized expression (not raw counts), the Bioconductor `timecourse` package (Tai & Speed 2006) provides a multivariate empirical Bayes approach for time course data. It models gene expression as a multivariate normal over time and borrows strength across genes to estimate the mean trajectory.

```r
# R/timecourse_analysis.R
# Appropriate for: microarray log-intensities or voom log-CPM, balanced design
run_timecourse_mb <- function(expr_matrix, meta, subject_col, time_col,
                              group_col = NULL, size = 1) {
  library(timecourse)
  # expr_matrix: genes x samples, log-scale
  # size: number of biological replicates per subject-time combination (MB statistic)
  
  times   <- sort(unique(meta[[time_col]]))
  n_times <- length(times)
  n_sub   <- length(unique(meta[[subject_col]]))

  # MBstat requires a genes x (subjects * times) matrix ordered by time within subject
  # Consult timecourse vignette for array permutation details
  mb <- mb.long(expr_matrix, times = n_times, reps = rep(size, n_sub))
  mb
}
```

See the `timecourse` package vignette for full details on `mb.long`, `mb.2D`, and how to use the HotellingT2 statistic for multivariate testing across time.

---

## 3. dream: Linear Mixed Models (variancePartition)

### Statistical model

dream fits a linear mixed model (LMM) for each gene:

$$\tilde{y}_{git} = \mathbf{x}_{it}^{\top} \boldsymbol{\beta}_g + \mathbf{z}_{it}^{\top} \mathbf{b}_{gi} + \epsilon_{git}$$

where $\mathbf{b}_{gi} \sim \mathcal{N}(0, \mathbf{D}_g)$ are gene-specific random effects estimated by REML. The key advantage over `duplicateCorrelation` is that $\mathbf{D}_g$ is estimated **per gene**, rather than imposing a single genome-wide consensus. Hypothesis tests use Satterthwaite or Kenward-Roger degrees of freedom approximations for small samples (via `lmerTest` and `pbkrtest`).

### Workflow

```r
# R/dream_analysis.R
run_dream <- function(dge, meta, fixed_formula, random_formula = NULL,
                      n_cores = 4) {
  library(variancePartition)
  library(edgeR)
  library(BiocParallel)

  register(SnowParam(n_cores))

  # Combine fixed + random effects into a single formula for dream
  # Example: form <- ~ time + condition + (1 | subject_id)
  form <- fixed_formula  # should already include random effects as (1|var)

  # voomWithDreamWeights instead of voom() — handles random effects in weights
  v <- voomWithDreamWeights(dge, form, meta, plot = TRUE)

  # Fit dream model
  fit <- dream(v, form, meta)
  fit <- variancePartition::eBayes(fit)

  list(voom = v, fit = fit)
}
```

**Important:** Use `voomWithDreamWeights()` (not `voom()`) and `variancePartition::eBayes()` and `variancePartition::topTable()` — dream replaces four limma functions. For fixed-effects-only models, these are equivalent to their limma counterparts.

### Variance partition analysis (optional but informative)

Before running dream, a variance partition analysis quantifies how much expression variance is explained by each variable. Use this to decide which variables belong as fixed versus random effects:

```r
# R/variance_partition.R
run_variance_partition <- function(voom_obj, meta, vp_formula, n_cores = 4) {
  library(variancePartition); library(BiocParallel)
  register(SnowParam(n_cores))
  vp <- fitExtractVarPartModel(voom_obj, vp_formula, meta)
  vp  # use plotVarPart(sortCols(vp)) in notebook
}
```

### Comparing dream vs duplicateCorrelation

When only one random effect is present and the contrast is simple, dream and `duplicateCorrelation` can give different results on a gene-by-gene basis. `duplicateCorrelation` uses a single genome-wide $\hat{\rho}$; if gene $g$'s true donor variance $\tau^2_g < \bar{\tau}^2$ (genome-wide mean), dream increases its p-value precision; if $\tau^2_g > \bar{\tau}^2$, dream is more conservative. Use `plotCompareP()` from variancePartition to visualize this.

---

## 4. DESeq2: Likelihood Ratio Test for Time Course

### Statistical model

DESeq2 models counts with a negative binomial GLM:

$$K_{gi} \sim \text{NB}(\mu_{gi},\, \alpha_g)$$
$$\log_2(\mu_{gi}) = \log_2(s_i) + \mathbf{x}_i^{\top} \boldsymbol{\beta}_g$$

where $s_i$ is the size factor (median-of-ratios normalization), $\alpha_g$ is the gene-specific dispersion, and $\boldsymbol{\beta}_g$ are log2 fold change coefficients. The **likelihood ratio test** compares a full model $M_1$ to a reduced model $M_0$ by the deviance difference:

$$\Lambda_g = 2(\log L(M_1) - \log L(M_0)) \xrightarrow{d} \chi^2_{df_1 - df_0}$$

This tests whether the terms removed in $M_0$ collectively explain significant variation — an ANODEV (analysis of deviance). The LRT is particularly useful for time course data because it produces one omnibus p-value testing whether expression differs across any timepoints, without requiring a specific pairwise contrast direction.

### Workflow: factor time (few timepoints)

```r
# R/deseq2_timecourse.R
run_deseq2_lrt_factor <- function(counts, meta, full_design, reduced_design,
                                   subject_col = NULL) {
  library(DESeq2)
  # Note: DESeq2 does not natively support random effects.
  # If subjects are a nuisance, include as a fixed effect only if balanced
  # and degrees of freedom permit. Otherwise use dream.

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = full_design
  )
  dds <- DESeq(dds, test = "LRT", reduced = reduced_design)
  dds
}
```

**Common design patterns:**

```r
# Single time factor (omnibus: any change over time?)
full    <- ~ timepoint
reduced <- ~ 1

# Time + condition, testing interaction (condition-specific time trajectories)
full    <- ~ condition + timepoint + condition:timepoint
reduced <- ~ condition + timepoint

# With a batch covariate
full    <- ~ batch + condition + timepoint + condition:timepoint
reduced <- ~ batch + condition + timepoint
```

### Workflow: continuous time with splines

For finely-sampled or continuous time, model the time trajectory with natural splines:

```r
# R/deseq2_spline_timecourse.R
run_deseq2_lrt_spline <- function(counts, meta, time_col, group_col = NULL,
                                   df = 4, covariate_cols = NULL) {
  library(DESeq2); library(splines)

  spline_mat <- ns(meta[[time_col]], df = df)
  colnames(spline_mat) <- paste0("ns", seq_len(df))
  meta_aug <- cbind(meta, as.data.frame(spline_mat))

  spline_terms <- paste(colnames(spline_mat), collapse = " + ")

  if (!is.null(group_col)) {
    # Test group:time interaction
    full_rhs    <- paste(c(covariate_cols, group_col, spline_terms,
                           paste0(group_col, ":", colnames(spline_mat))),
                         collapse = " + ")
    reduced_rhs <- paste(c(covariate_cols, group_col, spline_terms),
                         collapse = " + ")
  } else {
    full_rhs    <- paste(c(covariate_cols, spline_terms), collapse = " + ")
    reduced_rhs <- paste(c(covariate_cols, "1"), collapse = " + ")
  }

  full_design    <- as.formula(paste("~", full_rhs))
  reduced_design <- as.formula(paste("~", reduced_rhs))

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta_aug,
                                 design = full_design)
  dds <- DESeq(dds, test = "LRT", reduced = reduced_design)
  dds
}
```

**Note on repeated measures:** DESeq2 does not support random effects. For experiments with multiple measurements per subject, either (a) include subject as a fixed effect (only feasible with balanced designs and sufficient df), or (b) use dream or limma/voom + duplicateCorrelation instead.

---

## 5. Model Diagnostics — Critical for All Methods

Diagnostics are not optional. Always verify assumptions before interpreting results.

### 5a. voom mean-variance trend (limma/voom and dream)

The voom plot shows the mean-variance relationship of the log-counts. The LOWESS trend should be smoothly decreasing. Anomalies indicate problems:

- Concave/increasing trend at high expression: possible missing covariate or sample swap.
- Extremely high variances: outlier samples, incomplete filtering, or unmodelled batch effects.
- After fitting, `plotSA()` (sigma vs. mean log-expression) should show no residual trend.

```r
# In notebook: plot voom output
plot(v)          # mean-variance trend
plotSA(fit2)     # post-eBayes: should be flat
```

### 5b. Dispersion estimates (DESeq2)

```r
# In notebook
plotDispEsts(dds)  # gene-wise, fitted, and final (MAP) dispersions
# Expect: gene-wise estimates scatter around the fitted trend;
# extreme outliers (circled) are handled by DESeq2 Cook's distance filtering.
```

Cook's distances flag potential outlier samples on a per-gene basis. Check for systematic outlier samples:

```r
# Identify samples with high Cook's distances across many genes
W <- assays(dds)[["cooks"]]
boxplot(log10(t(W)), las=2, ylab="log10(Cook's distance)")
```

### 5c. Sample-level QC

Always examine PCA or MDS before modelling. Outlier samples should be investigated and potentially removed.

```r
# R/sample_qc.R
plot_pca <- function(voom_obj, meta, color_col, shape_col = NULL) {
  library(ggplot2)
  pca <- prcomp(t(voom_obj$E), scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df  <- data.frame(meta, PC1 = pca$x[,1], PC2 = pca$x[,2])
  gg  <- ggplot(df, aes(PC1, PC2, color = .data[[color_col]])) +
    geom_point(size = 3) +
    labs(x = paste0("PC1 (", pct[1], "%)"),
         y = paste0("PC2 (", pct[2], "%)"))
  if (!is.null(shape_col)) gg <- gg + aes(shape = .data[[shape_col]])
  gg
}
```

### 5d. P-value histograms

A well-behaved analysis produces a p-value histogram with a uniform distribution under the null plus a spike near 0 for true positives. Pathological shapes:

- Enrichment at p ≈ 1 (anti-conservative null): model is missing a covariate, or random effects are unaccounted for.
- Uniform with no spike: no signal, or severe underpowering.
- Bimodal near 0 and 1: possible data quality issue.

```r
# R/diagnostics.R
plot_pval_hist <- function(pvals, title = "P-value distribution") {
  hist(pvals, breaks = 50, xlab = "p-value", main = title,
       col = "steelblue", border = "white")
  abline(h = length(pvals) / 50, lty = 2, col = "red")
}
```

### 5e. Residual diagnostics (dream)

For dream models, examine per-gene residuals from the LMM fits. Systematic patterns in residuals vs. fitted values suggest model misspecification.

```r
# variancePartition provides residuals via residuals(fit)
# Spot-check a handful of top DE genes manually:
# plot(fitted(fit)["GENE",], residuals(fit)["GENE",])
```

### 5f. Convergence warnings (dream / lme4)

Mixed model fitting can fail to converge. Monitor `dream()` warnings:

```r
# After dream(), check:
sum(!is.finite(fit$F.value))  # genes with failed fits
# For systematic convergence failures, simplify the random effects structure
```

### 5g. Examining individual gene trajectories

Always visualize expression trajectories for top-ranked genes to ensure statistical results are biologically interpretable:

```r
# R/plot_gene_trajectory.R
plot_gene_trajectory <- function(expr, meta, gene, time_col, subject_col,
                                  group_col = NULL) {
  library(ggplot2)
  df <- data.frame(
    expr    = as.numeric(expr[gene, ]),
    time    = meta[[time_col]],
    subject = meta[[subject_col]]
  )
  if (!is.null(group_col)) df$group <- meta[[group_col]]

  gg <- ggplot(df, aes(time, expr, group = subject)) +
    geom_line(alpha = 0.4, color = "gray50") +
    geom_point(size = 2) +
    stat_summary(aes(group = 1), fun = mean, geom = "line",
                 color = "firebrick", linewidth = 1) +
    labs(title = gene, y = "log2-CPM", x = "Time")
  if (!is.null(group_col)) gg <- gg + aes(color = .data[[group_col]])
  gg
}
```

---

## 6. Multiple Testing and Results

All methods produce p-values requiring multiple testing correction (Benjamini-Hochberg FDR by default):

- limma: `topTable(fit, coef = ..., adjust.method = "BH", number = Inf)`
- dream: `variancePartition::topTable(fit, coef = ..., number = Inf)`
- DESeq2: results already include `padj` (IHW or BH); use `results(dds, name = ...)` or `results(dds, contrast = ...)`. With LRT, do **not** filter by fold change — the statistic is not directional.

For time series with many correlated tests (splines, LRT), consider independent hypothesis weighting (`IHW` package) or `qvalue` for better power.

---

## 7. Quarto Notebook Template

```yaml
---
title: "RNA-seq Time Course DE Analysis"
format:
  html:
    toc: true
    code-fold: true
execute:
  echo: false
  message: false
  warning: false
---
```

Source R scripts with `source("R/filter_genes.R")` etc. Display only key figures and result tables. Show code on demand via `code-fold: true`.

---

## 8. Normalization for Time Course Data

Normalization choices have an outsized impact in time series studies because systematic biases can masquerade as temporal trajectories. See `references/oh-li-2021-review.md` Section 4 for the full comparison. Key guidance:

**Always start from raw integer counts.** FPKM and TPM normalize for gene length and are intended for within-sample comparisons only — they must not be used as input to limma/voom, dream, or DESeq2.

**Between-sample normalization methods:**

| Method | Function | Assumption | Use when |
|---|---|---|---|
| **TMM** | `edgeR::calcNormFactors(method="TMM")` | Most genes non-DE | Default for voom-based pipelines |
| **RLE / median-of-ratios** | Internal to `DESeq2::estimateSizeFactors()` | Most genes non-DE | Default for DESeq2 |
| **Upper Quartile (UQ)** | `edgeR::calcNormFactors(method="upperquartile")` | Upper quartile stable | Large fraction of transcriptome changing |
| **Total Count / CPM** | — | Library size only bias | Not recommended for DE |

TMM and RLE perform similarly in most settings and are appropriate defaults. UQ is worth considering for strong developmental transitions where a large fraction of genes change — in those cases TMM/RLE can break down because their "most genes non-DE" assumption is violated.

**Pre-filtering** must happen before normalization. Low-count genes inflate dispersion estimates and reduce power. Use `edgeR::filterByExpr()` or the CPM threshold approach in Section 1.

---

## 9. Batch Correction for Time Course Data

Batch effects are especially dangerous in time course studies — an unbalanced batch factor can create or completely mask a temporal trend. See `references/oh-li-2021-review.md` Section 5 for full method details.

### Decision: model vs. pre-correct

**Preferred approach:** include batch as a fixed covariate in the design matrix. This properly propagates uncertainty and avoids double-correction artifacts:

```r
# limma/voom: add batch to design formula
design <- model.matrix(~ batch + timepoint + condition, data = meta)

# DESeq2: add batch to both full and reduced designs
full    <- ~ batch + timepoint + condition:timepoint
reduced <- ~ batch + timepoint
```

**When to pre-correct counts:** only when the downstream method cannot accept covariates (e.g., some clustering algorithms), or when ComBat-Seq is used (which preserves count properties and is safe as input to DE methods).

### Batch factor known

**ComBat-Seq** (`sva::ComBat_seq`) — operates on raw counts using a NB GLM; preferred over the original ComBat (which assumes log-normal data):

```r
# R/batch_correction.R
correct_batch_combatseq <- function(counts, meta, batch_col, group_col = NULL) {
  library(sva)
  group <- if (!is.null(group_col)) meta[[group_col]] else NULL
  ComBat_seq(
    counts = counts,
    batch  = meta[[batch_col]],
    group  = group
  )
}
```

**Harman** (`Harman::harman`) — operates on pre-normalized log-scale data; appropriate for voom-transformed expression matrices or microarray data.

### Batch factor unknown

**svaseq** (`sva::svaseq`) — estimates surrogate variables (SVs) from residuals; include SVs as covariates in the design matrix (do not pre-correct counts):

```r
# R/batch_correction.R
estimate_surrogate_variables <- function(counts_norm, meta, full_formula,
                                          null_formula = ~ 1) {
  library(sva)
  mod  <- model.matrix(full_formula, data = meta)
  mod0 <- model.matrix(null_formula, data = meta)
  svobj <- svaseq(log1p(counts_norm), mod, mod0)
  # Returns svobj$sv: matrix of n_samples x n_sv
  # Bind columns to meta and add to design formula
  svobj
}
```

**RUVSeq** (`RUVSeq::RUVr`, `RUVg`, `RUVs`) — factor analysis of SVD; use RUVr when no negative control genes are available. Resulting W factors should be included as design covariates.

### Diagnostics

Always run PCA/MDS before and after batch correction to verify the batch structure was removed without distorting the biological signal:

```r
# Use plot_pca() from R/diagnostics.R — call once on uncorrected,
# once on corrected data, and compare visually.
# A good correction: batch cluster separation disappears;
# condition/timepoint separation is preserved or enhanced.
```

**Warning:** `limma::removeBatchEffect()` is appropriate for visualization only. Never use its output as input to a DE model — it removes degrees of freedom that the model needs and inflates false positives.

### The review's recommendation for time course data

From Oh & Li (2021): use **ComBat-Seq, svaseq, and RUVSeq** as a diagnostic panel on pre-filtered and normalized data. For actual correction prior to dynamic DE analysis, use **ComBat-Seq** (known batch) or **Harman** (known batch, log-scale input). When batch factors are unknown, include svaseq or RUVSeq surrogate variables in the design matrix.

---

## 10. References

- **Oh & Li (2021)** — comprehensive review of dynamic methods, normalization, and batch correction for time course RNA-seq:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7997275/ (synthesized in `references/oh-li-2021-review.md`)
- **Spies et al. (2019)** — large-scale benchmark of time course DE methods:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6357553/
- **limma User's Guide** Ch. 15 (RNA-seq with voom), Ch. 18 (time series / longitudinal):
  https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
- **timecourse package**: Tai & Speed (2006). Biometrics.
  https://www.bioconductor.org/packages/release/bioc/vignettes/timecourse/inst/doc/timecourse.pdf
- **dream**: Hoffman & Roussos (2021). *Bioinformatics* 37(2):192–201.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8055218/
- **variancePartition / dream vignette**:
  https://bioconductor.org/packages/release/bioc/html/variancePartition.html
- **Law et al. 2014** (voom): *Genome Biology* 15:R29.
- **DESeq2 time course workflow**:
  https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
- **DESeq2 LRT**: Love, Huber & Anders (2014). *Genome Biology* 15:550.
