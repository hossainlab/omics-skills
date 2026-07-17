# omics-dataset-retrieval

Systematically retrieve, deduplicate, classify, and relevance-audit publicly
available omics datasets for a disease, phenotype, gene, or biological process.

- **Domain:** data retrieval
- **Skill dir:** [`skills/omics-dataset-retrieval`](../skills/omics-dataset-retrieval)
- **Boundary:** catalogs metadata only — does **not** download raw data files or
  run downstream analysis.

## What it does

Given a topic (e.g. a disease or gene) plus synonyms, the skill queries the
broadest reachable set of public repositories, normalizes and deduplicates hits
across them, tags each dataset by omics type, scores its relevance to the topic,
and writes a catalog plus a summary report.

Covered omics types: transcriptomics (bulk, single-cell, spatial), proteomics,
metabolomics/lipidomics, epigenomics, genomics, and multi-omics.

## Repositories queried

Worked in priority tiers — all tiers are attempted, since rare topics often only
have data in Tier 2/3.

- **Tier 1 (programmatic APIs, always):** GEO, SRA, ArrayExpress/BioStudies,
  PRIDE/ProteomeXchange, OmicsDI, CZ CELLxGENE, GDC/TCGA.
- **Tier 2 (specialized, when relevant):** ENCODE, Expression Atlas, HCA,
  Metabolomics Workbench, MetaboLights, MassIVE/GNPS, jPOST, iProX, cBioPortal,
  EpiRR/IHEC, HPA, ENA.
- **Tier 3 (web search):** Zenodo, Figshare, Dryad, OSF, Harvard Dataverse,
  Synapse, ICGC, CPTAC, AWS Open Data, dbGaP/EGA/JGA (controlled), UK Biobank,
  FinnGen.

## Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `disease_or_topic` | string | Disease, phenotype, gene, or process |
| `synonyms` | list[str] | Alternative names / abbreviations / gene symbols (more = better recall) |
| `omics_types` | list[str] or `"all"` | Restrict to omics types, or all (default) |
| `organism` | string | `"all"` (default), or restrict to a species on request |
| `year_min` | int | Earliest publication year (default: no limit) |
| `output_dir` | path | Where outputs are written (default: `/mnt/results/`) |

## Outputs

| File | Description |
|------|-------------|
| `<topic>_omics_datasets_MASTER.csv` | Full catalog with metadata + relevance labels |
| `<topic>_omics_datasets_VALIDATED.csv` | Filtered to `CORE/DIRECT` + `ADJACENT` |
| `<topic>_omics_summary.md` | Counts by omics type / repository, top datasets, limitations |
| `<topic>_omics_landscape.png` | Optional: donut + repository bar + timeline |

## Relevance tiers

Each dataset is scored into one of four tiers by keyword analysis of title +
summary:

| Tier | Meaning |
|------|---------|
| `CORE/DIRECT` | Real patient samples or validated disease model |
| `ADJACENT` | Mechanistically relevant, no patient data |
| `WEAK` | Epidemiologically linked or tangential |
| `REMOVE` | False positive (always manually reviewed before exclusion) |

## How to use

Describe the retrieval task in natural language — the skill activates on a match:

```
> Find all publicly available omics datasets for Alzheimer's disease
> What RNA-seq data is available for BCL11A in erythroid cells?
> Survey all cancer proteomics datasets for pancreatic ductal adenocarcinoma
```

The skill first asks a short set of clarification questions (synonyms, omics
types, organism, year range, controlled-access inclusion, output goal, tissue
focus) before running any queries — answers materially change which repositories
are searched and how relevance is tuned.

## Package layout

```
skills/omics-dataset-retrieval/
├── SKILL.md                       # Agent instructions (17-step workflow)
├── references/
│   ├── repositories.md            # Tier 1/2/3 map + API limitations
│   ├── query-code.md              # Per-repository query snippets
│   ├── classification.md          # Omics classifier + relevance audit + assembly
│   └── examples.md                # Worked examples
└── scripts/
    └── omics_retrieval.py         # Pure helpers (classify/audit/dedup)
```

## Tests

`tests/omics-dataset-retrieval/` (pytest) covers the deterministic core:
omics-type classification, the 4-tier relevance audit, and cross-repository
accession deduplication. Network query code is not unit-tested.

```bash
pytest tests/omics-dataset-retrieval/
```

## Key caveats

- GEO's `gdstype` field is unreliable — omics type is inferred from title +
  summary instead.
- The same study often appears under multiple accessions (e.g. `GSE41575` and its
  ArrayExpress mirror `E-GEOD-41575`); these are normalized before deduplication.
- Relevance scoring is keyword-based, so it produces false positives and
  negatives — `REMOVE` candidates should be reviewed manually.
- Controlled-access datasets (dbGaP/EGA/JGA) are catalogued and labelled but
  require institutional data-access agreements to download.
