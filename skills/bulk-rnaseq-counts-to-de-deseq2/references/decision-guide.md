# DESeq2 Decision Guide

Decision trees for choosing between DESeq2 methods and approaches.

## Decision 1: Transformation Method

When: after `DESeq()`, before PCA/heatmaps. Question: `vst()` or `rlog()`?

### Option A: VST (Variance Stabilizing Transformation)

Use when: n > 30 samples, need fast computation.

Pros: fast (1000+ samples OK), suitable for large datasets. Cons: less accurate
for small samples (n < 10).

```r
vsd <- vst(dds, blind = FALSE)  # blind=FALSE uses design (recommended)
```

### Option B: rlog (Regularized Log)

Use when: n < 30 samples, want best stabilization.

Pros: better for small samples, better at low counts. Cons: slow for large
datasets (>100 samples).

```r
rld <- rlog(dds, blind = FALSE)
```

### Decision Tree

```
n > 30 samples?
├─ YES → vst() (fast, appropriate for large datasets)
└─ NO → n < 10?
        ├─ YES → rlog() (better stabilization for small samples)
        └─ NO (10-30) → Either works (vst for speed, rlog for accuracy)
```

When to use `blind = TRUE`: exploratory analysis without design, initial QC, want
natural clustering.

## Decision 2: LFC Shrinkage Method

When: after `results()`, before ranking/plotting. Question: which shrinkage method?

### Option A: apeglm (Recommended)

Use when: ranking genes, publication plots, want best performance.

Pros: best shrinkage, preserves large LFC, accurate posteriors. Cons: requires
`coef` (not contrast), needs package, slightly slower.

```r
library(apeglm)
resLFC <- lfcShrink(dds, coef = resultsNames(dds)[2], type = 'apeglm')
```

Note: cannot use contrast specification.

```r
# Works:
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# Doesn't work:
# res <- lfcShrink(dds, contrast = c('condition', 'treated', 'control'), type = 'apeglm')
```

### Option B: ashr

Use when: large datasets (apeglm slow), need contrast specification.

Pros: good performance, fast, works with contrasts. Cons: may over-shrink large
FC, inferior to apeglm for ranking.

```r
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'ashr')
# Also works with contrasts
resLFC <- lfcShrink(dds, contrast = c('condition', 'treated', 'control'), type = 'ashr')
```

### Option C: normal (Legacy)

Use when: backward compatibility, no apeglm/ashr.

Pros: fast, simple, no extra packages. Cons: inferior shrinkage, not recommended
for new analyses.

```r
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'normal')
```

### Shrinkage vs Unshrunk

Use SHRUNK for:

- ✓ Ranking genes by effect size
- ✓ Visualization (volcano, MA plots, heatmaps)
- ✓ GSEA
- ✓ Selecting top genes for validation

Use UNSHRUNK for:

- ✓ Hypothesis testing (p-values)
- ✓ Reporting fold changes with CI
- ✓ Testing specific FC thresholds

```r
res <- results(dds)  # Unshrunk: hypothesis testing
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')  # Shrunk: ranking/viz
```

## Decision 3: Design Formula

When: before creating DESeqDataSet. Question: include covariates in design?

### Option A: Simple (`~ condition`)

Use when: no batch effects, single factor, exploratory.

Pros: straightforward, easy to explain, maximum power for simple comparisons.
Cons: lower power if batch effects present, can't control confounders.

```r
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
```

Check first:

```r
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")
# Samples cluster by condition? → Use simple design
# Samples cluster by batch? → Use multi-factor design
```

### Option B: Multi-Factor (`~ batch + condition`)

Use when: known batch effects, sequencing runs, PCA shows batch clustering.

Pros: controls confounders, higher power with batch, more accurate FC. Cons:
requires balanced design, uses degrees of freedom, complex interpretation.

```r
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ batch + condition)
res <- results(dds, name = 'condition_treated_vs_control')
```

Critical: design must not be confounded.

```r
# ✅ Works (balanced):
# batch 1: control, control, treated, treated
# batch 2: control, control, treated, treated

# ❌ Fails (confounded):
# batch 1: control, control
# batch 2: treated, treated
```

### Option C: Paired (`~ individual + condition`)

Use when: same individuals before/after, matched pairs (tumor/normal).

Pros: controls individual effects, higher power, reduces noise. Cons: requires
paired structure, can't test between-individual effects.

```r
coldata <- data.frame(
  individual = factor(rep(1:5, each = 2)),
  condition = factor(rep(c('before', 'after'), 5))
)
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ individual + condition)
```

### Option D: Interaction (`~ genotype * treatment`)

Use when: test if treatment effect differs by genotype, gene × environment.

Pros: tests interactions, identifies genotype-specific responses. Cons: requires
more samples (n ≥ 4 per group), complex interpretation.

```r
design = ~ genotype + treatment + genotype:treatment
dds <- DESeq(dds)

# Main effects
res_genotype <- results(dds, name = 'genotype_KO_vs_WT')
res_treatment <- results(dds, name = 'treatment_drug_vs_vehicle')

# Interaction: Does treatment effect differ by genotype?
res_interaction <- results(dds, name = 'genotypeKO.treatmentdrug')
```

### Decision Tree

```
Known batch effects or technical covariates?
├─ YES → Balanced design (batches in both conditions)?
│        ├─ YES → ~ batch + condition
│        └─ NO → Cannot adjust (confounded); consider batch correction
│
└─ NO → Paired samples (before/after, tumor/normal)?
        ├─ YES → ~ individual + condition
        └─ NO → Testing interactions?
                ├─ YES → ~ genotype * treatment
                └─ NO → ~ condition
```

### Checking Design

```r
# Check if full rank (not confounded)
design_matrix <- model.matrix(~ batch + condition, coldata)
Matrix::rankMatrix(design_matrix) == ncol(design_matrix)  # Should be TRUE

# Check PCA to see if batch correction helps
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("condition", "batch"))
```

## Decision 4: Pre-Filtering Strategy

When: after creating DESeqDataSet, before `DESeq()`.

### Option A: Minimum Total Counts (Standard)

```r
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

Use: standard approach for most datasets.

### Option B: Minimum Samples with Counts

```r
keep <- rowSums(counts(dds) >= 10) >= 3  # ≥3 samples with 10+ counts
dds <- dds[keep,]
```

Use: small sample size, want genes expressed in multiple samples.

### Option C: Mean Expression

```r
keep <- rowMeans(counts(dds)) >= 10
dds <- dds[keep,]
```

Use: large sample size, filter by average expression.

Recommendation: use Option A for most analyses — simple and effective.

## Decision 5: Significance Thresholds

When: filtering results.

### Standard Thresholds

| Threshold | Use Case | Stringency |
|-----------|----------|------------|
| `padj < 0.1` | DESeq2 default, discovery | Relaxed |
| `padj < 0.05` | Standard, most studies | Moderate |
| `padj < 0.01` | High confidence | Stringent |

### Fold Change Thresholds

| Threshold | FC | Use Case |
|-----------|-----|----------|
| `|log2FC| > 0.5` | 1.4× | Small effect |
| `|log2FC| > 1` | 2× | Standard |
| `|log2FC| > 2` | 4× | Large effect |

### Examples

```r
# Standard
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Discovery
sig <- subset(res, padj < 0.1)

# High confidence
sig <- subset(res, padj < 0.01 & abs(log2FoldChange) > 2)
```

## Quick Reference Table

| Decision | Use This | When |
|----------|----------|------|
| Transformation | `vst()` | n > 30 samples |
| | `rlog()` | n < 30 samples |
| LFC Shrinkage | `apeglm` | Ranking/visualization |
| | `ashr` | Need contrasts/speed |
| Design | `~ condition` | Simple, no batch |
| | `~ batch + condition` | Known batch effects |
| | `~ individual + condition` | Paired samples |
| | `~ genotype * treatment` | Interactions |
| Pre-filtering | `rowSums >= 10` | Standard |
| Significance | `padj < 0.05 & |log2FC| > 1` | Standard |

## Best Practices

- Always check PCA first before finalizing design formula.
- Use shrinkage for visualization, unshrunk for testing.
- Document decisions in analysis code.
- Report both unshrunk and shrunk results.
- Validate top hits with independent method.
