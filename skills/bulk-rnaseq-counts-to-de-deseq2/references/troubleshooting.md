# DESeq2 Troubleshooting Guide

Solutions to common DESeq2 errors and issues.

## Error Messages

### "the model matrix is not full rank"

Error:

```
Error in checkFullRank(modelMatrix): the model matrix is not full rank
```

Cause: confounded variables (one variable perfectly predicts another).

Common example:

```r
# ❌ All treated in batch B, all control in batch A
coldata <- data.frame(
  condition = c('control', 'control', 'treated', 'treated'),
  batch = c('A', 'A', 'B', 'B')
)
design = ~ batch + condition  # Confounded!
```

Solutions:

Remove confounded variable:

```r
design = ~ condition  # Remove batch
```

Combine into single factor:

```r
coldata$group <- factor(paste(coldata$batch, coldata$condition, sep = '_'))
design = ~ 0 + group
```

Check design rank:

```r
design_matrix <- model.matrix(~ batch + condition, coldata)
Matrix::rankMatrix(design_matrix) == ncol(design_matrix)  # Should be TRUE
```

### "counts matrix should contain integer values"

Error:

```
Error: some values in assay are not integers
```

Cause: using normalized data (TPM/FPKM/RPKM) or tximport without proper function.

Solutions:

Use raw counts — DESeq2 requires integers, not normalized data.

For tximport data:

```r
library(tximport)
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
# NOT DESeqDataSetFromMatrix
```

Round if needed (only if counts have minor decimals from artifacts):

```r
counts <- round(counts)
dds <- DESeqDataSetFromMatrix(counts, coldata, design)
```

### "every gene contains at least one zero"

Error:

```
Error in estimateDispersionsGeneEst: every gene contains at least one zero
```

Cause: usually transposed matrix (samples as rows instead of columns).

Solutions:

Check orientation:

```r
dim(counts)  # Should be genes × samples
head(counts)

# If samples are rows, transpose:
counts <- t(counts)
```

Verify data loaded correctly:

```r
sum(counts)  # Should be > 0
colSums(counts)  # Check sample totals
rowSums(counts)  # Check gene totals
```

### "factor levels not in colData" / "subscript out of bounds"

Cause: sample name mismatch or typo in design formula.

Solutions:

Check names match:

```r
colnames(counts)  # Column names in counts
rownames(coldata)  # Row names in colData
all(colnames(counts) == rownames(coldata))  # Should be TRUE
```

Fix mismatch:

```r
# Reorder coldata to match counts
coldata <- coldata[colnames(counts), ]
```

Check design formula:

```r
colnames(colData(dds))  # See available columns
design = ~ condition  # Check spelling matches
```

## Results Issues

### Too many NA values in padj

Observation: many genes have `padj = NA`.

Cause: normal DESeq2 behavior — independent filtering and Cook's distance outlier
detection.

Solutions:

This is expected — DESeq2 sets NA for:

- Low count genes (independent filtering)
- Outlier samples (Cook's distance)
- Genes with extreme outliers

Adjust filtering:

```r
res <- results(dds, alpha = 0.1)  # More lenient
# Or disable
res <- results(dds, independentFiltering = FALSE)
```

Check specific genes:

```r
mcols(res)$filterThreshold
which(is.na(res$padj))
```

### No significant genes found

Observation: `summary(res)` shows 0 genes with `padj < 0.1`.

Possible causes: no real DE, high variability, small sample size, batch effects,
wrong reference.

Solutions:

Check PCA:

```r
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = 'condition')
# Samples cluster by condition? If not, may not have DE
```

Check batch effects:

```r
plotPCA(vsd, intgroup = c('condition', 'batch'))
# If clustered by batch, add to design:
design = ~ batch + condition
```

Verify reference level:

```r
levels(dds$condition)  # First is reference
dds$condition <- relevel(dds$condition, ref = 'control')
```

Check dispersion:

```r
plotDispEsts(dds)
# High dispersion = high variability = low power
```

Relax threshold:

```r
sig <- subset(res, padj < 0.1 & abs(log2FoldChange) > 0.5)
```

Increase sample size — n ≥ 4 per group recommended.

### Extremely large fold changes (log2FC > 10)

Observation: some genes show >1000-fold change.

Cause: low count genes (dividing by near-zero creates large FC).

Solutions:

Use LFC shrinkage:

```r
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
```

Filter low counts:

```r
keep <- rowMeans(counts(dds)) >= 10
dds <- dds[keep,]
```

Check specific genes:

```r
plotCounts(dds, gene = 'gene_with_high_fc', intgroup = 'condition')
# If counts are 0 vs 1, not reliable
```

### Unexpected up/downregulation direction

Observation: expected treated > control, but `log2FC` is negative.

Cause: wrong reference level or contrast direction reversed.

Solutions:

Check reference:

```r
levels(dds$condition)  # First is reference (denominator)
dds$condition <- relevel(dds$condition, ref = 'control')
dds <- DESeq(dds)
```

Understand direction:

```r
# log2FC = log2(numerator / denominator)
# For condition_treated_vs_control:
# log2FC > 0: treated > control (upregulated in treated)
# log2FC < 0: treated < control (downregulated in treated)
```

Flip contrast if needed:

```r
res <- results(dds, contrast = c('condition', 'control', 'treated'))
```

## Performance Issues

### DESeq() is very slow

Cause: large dataset (>50K genes, >100 samples).

Solutions:

Pre-filter more aggressively:

```r
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
```

Use parallel processing:

```r
library(BiocParallel)
register(MulticoreParam(4))  # Use 4 cores
dds <- DESeq(dds, parallel = TRUE)
```

Simplify design:

```r
design = ~ condition  # Instead of complex interaction
```

### rlog() takes forever

Cause: `rlog()` is slow for large datasets.

Solution: use `vst()` instead.

```r
vsd <- vst(dds, blind = FALSE)  # Much faster, similar results
```

## Interpretation Issues

### High dispersion estimates

Observation: `plotDispEsts()` shows points scattered far from fitted line.

Cause: high biological variability, batch effects, or outlier samples.

Solutions:

Check for outliers:

```r
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = 'condition')
# Remove clear outliers if present
```

Add batch to design:

```r
design = ~ batch + condition
```

Accept high variability — some systems naturally have high variance; results still
valid, just lower power.

### Cook's distance outliers

Observation: many genes filtered due to Cook's distance.

Cause: outlier samples with extreme counts.

Solutions:

Inspect outliers:

```r
W <- res$stat
outliers <- apply(W, 1, function(x) any(abs(x) > 3))
```

Remove problem samples:

```r
dds <- dds[, !colData(dds)$sample == 'outlier_sample']
dds <- DESeq(dds)
```

Disable filtering (use with caution):

```r
res <- results(dds, cooksCutoff = FALSE)
```

## Best Practices to Avoid Issues

Check data structure:

- Verify counts are integers
- Check sample/gene names match
- Confirm orientation (genes × samples)

Run QC first:

- PCA to check clustering
- Identify batch effects
- Find outlier samples

Set reference explicitly:

- Use `relevel()` before `DESeq()`
- Don't rely on alphabetical order

Document design:

- Comment why you chose design formula
- Record batch corrections

Start simple:

- Begin with `~ condition`
- Add complexity only if needed

Report versions:

```r
packageVersion('DESeq2')
sessionInfo()
```

## Common Error Quick Reference

| Error | Likely Cause | Quick Fix |
|-------|--------------|-----------|
| "not full rank" | Confounded design | Remove confounded variable |
| "should be integers" | Normalized data | Use raw counts |
| "every gene zero" | Transposed matrix | `counts <- t(counts)` |
| "subscript out of bounds" | Name mismatch | `coldata <- coldata[colnames(counts), ]` |
| Many NA in padj | Normal filtering | Expected behavior |
| No significant genes | No DE / wrong ref | Check PCA, verify reference level |
| Huge fold changes | Low counts | Use LFC shrinkage |
| DESeq() slow | Large dataset | Pre-filter, use parallel |

## Getting Help

Check first:

- DESeq2 vignette: `browseVignettes('DESeq2')`
- Bioconductor support: https://support.bioconductor.org/

When asking for help, provide:

- Full error message
- `sessionInfo()` output
- Code that produces error
- Sample size and experimental design
