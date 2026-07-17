# bulk-rnaseq-counts-to-de-deseq2

Differential expression analysis on bulk RNA-seq count data with
[DESeq2](https://bioconductor.org/packages/DESeq2/) (R/Bioconductor).

- **Domain:** transcriptomics
- **Skill dir:** [`skills/bulk-rnaseq-counts-to-de-deseq2`](../skills/bulk-rnaseq-counts-to-de-deseq2)
- **Language:** R (Bioconductor)

## What it does

Takes a raw count matrix and sample metadata and runs the standard DESeq2
workflow: build a `DESeqDataSet`, pre-filter low-count genes, set the reference
level, fit the negative-binomial GLM, extract results by coefficient or contrast,
apply log-fold-change shrinkage, and export significant genes and QC plots.

Handles simple two-group comparisons through multi-factor, paired, interaction,
and likelihood-ratio-test designs.

## Requirements

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'apeglm'))
# QC plots additionally use: ggplot2, ggprism, ggrepel (and optionally svglite)
```

## Inputs

| Input | Format | Description |
|-------|--------|-------------|
| Count matrix | Integer matrix | Genes (rows) × samples (columns); raw counts, not normalized |
| Sample metadata | data.frame | Row names match count column names |
| Design formula | R formula | Variables from metadata, e.g. `~ condition` or `~ batch + condition` |

Alternative inputs are supported: `SummarizedExperiment`, `tximport`
(Salmon/Kallisto), and `featureCounts`.

## Outputs

Depending on the scripts invoked:

- Results tables (all genes, significant genes, top-N) as CSV
- Size-factor normalized counts and VST/rlog-transformed counts as CSV
- The saved `DESeqDataSet` object (`.rds`)
- QC figures (dispersion, PCA, MA, volcano) as PNG + SVG

## Design formulas

| Formula | Use |
|---------|-----|
| `~ condition` | Single factor, no batch effects |
| `~ batch + condition` | Known batch effects (must be balanced, not confounded) |
| `~ individual + condition` | Paired samples (before/after, tumor/normal) |
| `~ genotype * treatment` | Test whether treatment effect differs by genotype |
| `~ sex + age_group + treatment` | Multiple confounders |

## How to use

Describe the analysis in natural language:

```
> Run DESeq2 on my count matrix with treated vs control comparison
> Analyze differential expression controlling for batch effects
> Get significantly differentially expressed genes with padj < 0.05
```

Or run the bundled scripts directly in R (paths relative to the skill directory):

```r
source("scripts/load_example_data.R")   # example data + input validation
source("scripts/basic_workflow.R")      # end-to-end walkthrough
source("scripts/extract_results.R")     # contrasts, shrinkage, filtering
source("scripts/qc_plots.R")            # dispersion / PCA / MA / volcano
source("scripts/export_results.R")      # write tables and transformed counts
```

## Package layout

```
skills/bulk-rnaseq-counts-to-de-deseq2/
├── SKILL.md                        # Agent instructions + DESeq2 reference
├── references/
│   ├── comprehensive-reference.md  # Full DESeq2 code patterns
│   ├── decision-guide.md           # Method-selection decision trees
│   ├── troubleshooting.md          # Common errors and fixes
│   └── usage-guide.md              # Quick prompts + input requirements
└── scripts/
    ├── basic_workflow.R            # Complete 8-step workflow
    ├── batch_correction.R          # ~ batch + condition example
    ├── multi_condition.R           # Three-group comparisons
    ├── extract_results.R           # Coefficient/contrast extraction, filtering
    ├── export_results.R            # Export tables, counts, objects
    ├── transformations.R           # VST / rlog helpers
    ├── qc_plots.R                  # Dispersion, PCA, MA, volcano
    └── load_example_data.R         # pasilla/airway loaders + validation
```

## Key guidance

- Use **raw** integer counts — never TPM/FPKM/RPKM.
- Pre-filter low-count genes before `DESeq()`.
- Set the reference level explicitly with `relevel()`.
- Use LFC shrinkage (`apeglm`) for ranking and visualization; use **unshrunk**
  results for hypothesis testing.
- Use `vst()` for large datasets (>30 samples), `rlog()` for small.
- Always use `padj` (adjusted p-value), never raw `pvalue`.
- Check PCA before finalizing the design formula, to spot batch effects.

See `references/decision-guide.md` for method selection and
`references/troubleshooting.md` for common errors.

## Tests

No automated tests yet — R scripts require an R + Bioconductor runner
(`testthat`), separate from the repository's Python test suite.
