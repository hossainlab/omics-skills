---
name: omics-dataset-retrieval
description: >-
  Systematically retrieve, deduplicate, classify, and relevance-audit publicly
  available omics datasets for a user-specified disease, phenotype, gene, or
  biological process. Covers all major omics types (transcriptomics, proteomics,
  metabolomics, epigenomics, genomics, single-cell, spatial, lipidomics,
  multi-omics) across the broadest set of public repositories (GEO, SRA,
  ArrayExpress, PRIDE, OmicsDI, CELLxGENE, GDC/TCGA, ENCODE, and 20+ more). Use
  when the user asks to find, survey, catalog, or map available omics/sequencing
  datasets for a topic. Does NOT download raw data or run downstream analysis.
---

# Omics Dataset Retrieval

Build a deduplicated, relevance-audited catalog of public omics datasets for a
disease, phenotype, gene, or biological process. Query every reachable public
repository, classify each hit by omics type, score its relevance, and emit a CSV
catalog + Markdown summary (+ optional landscape figure).

**Boundary:** catalog metadata only. This skill does not download raw data files
and does not perform downstream analysis.

Worked examples (Alzheimer's, BCL11A/erythroid, PDAC) are in
`references/examples.md`.

## Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `disease_or_topic` | string | Disease, phenotype, gene, or process (e.g. `"sickle cell disease"`, `"BCL11A"`) |
| `synonyms` | list[str] | Alternative names, abbreviations, gene symbols (e.g. `["SCD","SCA","HbSS"]`) |
| `omics_types` | list[str] or `"all"` | Restrict to omics types, or `"all"` (default) |
| `organism` | string | `"all"` (default) — human, mouse, and all others. Restrict to `"human"`/`"mouse"` only if explicitly requested |
| `year_min` | int | Earliest publication year (default: no limit) |
| `output_dir` | path | Where to write outputs (default: `/mnt/results/`) |

## Outputs

| File | Description |
|------|-------------|
| `<disease>_omics_datasets_MASTER.csv` | Full catalog with all metadata + relevance labels |
| `<disease>_omics_datasets_VALIDATED.csv` | Filtered to `CORE/DIRECT` + `ADJACENT` only |
| `<disease>_omics_summary.md` | Counts by omics type / repository, top datasets, limitations |
| `<disease>_omics_landscape.png` | Optional overview: donut + repository bar + timeline |

## Workflow

### Step 1 — Ask clarification questions (MANDATORY before any search)

Not optional. Before querying, use `AskUserQuestion` to collect the items below
in a **single call**. Skip only questions already answered in the user's initial
message; always ask the rest. Answers materially change which repositories run,
result counts, and relevance tuning.

- **Q1 Disease/topic** — name; synonyms/abbreviations (more = better recall);
  key genes/pathways (drive extra targeted GEO queries).
- **Q2 Omics types** — all (default) / transcriptomics only / epigenomics only /
  proteomics+metabolomics only / custom.
- **Q3 Organism** — all (default) / human only / specific species.
- **Q4 Year range** — no restriction (default) / recent only / custom.
- **Q5 Controlled-access repos** (dbGaP, EGA, JGA) — include+flag (default) /
  open-access only.
- **Q6 Goal** — browse & select (default) / landscape overview (figure) /
  download & analyze (add links) / all.
- **Q7 Tissue/cell-type focus** (optional, big recall boost) — e.g.
  `"whole blood"`, `"brain cortex"`, `"CD34+ HSCs"`. Added as targeted terms.

Proceed to Step 2 only once answers collected.

### Step 2 — Build search term matrix

Combine: primary name + synonyms × omics-type terms × tissue/cell terms ×
clinical-context terms (treatment, crisis, biomarker, pediatric, longitudinal) ×
key gene/pathway terms.

GEO syntax tips: `"<disease>"[Title]` for precision; `[DataSet Type]` for omics
filtering; `retmax=100` per query; run 20–40 queries; dedupe by GEO UID.

### Step 3–12 — Query repositories

Work through **all tiers** in priority order — rare diseases often only have data
in Tier 2/3. Full API base URLs, per-repository query code, and quirks are in the
reference files:

- **`references/repositories.md`** — Tier 1/2/3 coverage map + known API
  limitations + workarounds.
- **`references/query-code.md`** — runnable query snippets per repository (GEO,
  SRA, ArrayExpress/BioStudies, PRIDE, jPOST, OmicsDI, Metabolomics Workbench,
  CELLxGENE, GDC/TCGA, ENCODE, Expression Atlas, Tier-3 web search).

Tier 1 (always query): GEO, SRA, ArrayExpress/BioStudies, PRIDE/ProteomeXchange,
OmicsDI, CZ CELLxGENE, GDC/TCGA. Tier 2 (when relevant): ENCODE, Expression Atlas,
HCA, Metabolomics Workbench, MetaboLights, MassIVE/GNPS, jPOST, iProX, cBioPortal,
EpiRR/IHEC, HPA, ENA. Tier 3 (WebSearch tool): Zenodo, Figshare, Dryad, OSF,
Harvard Dataverse, Synapse, ICGC, CPTAC, AWS Open Data, dbGaP/EGA/JGA (controlled),
UK Biobank, FinnGen.

GEO is the primary source — run 20–40 targeted queries. Use OmicsDI as the single
entry point for metabolomics (aggregates MetaboLights, MassIVE, GNPS, Metabolomics
Workbench). Only query GDC/TCGA/CPTAC/cBioPortal for cancer topics; only ENCODE
when the topic has epigenomic/regulatory components.

### Step 13 — Classify omics type

Apply rule-based `classify_omics()` to each record. The GEO `gdstype` field is
**unreliable** — always use Title + Summary as the primary signal (check
single-cell/spatial before bulk RNA-seq; ArrayExpress labels arrays as
"transcription profiling by array"). Function in
`references/classification.md`.

### Step 14 — Relevance audit (4 tiers)

Score each record with `audit_relevance()` using disease/weak/exclude keyword
lists tuned per topic:

| Tier | Meaning | Action |
|------|---------|--------|
| `CORE/DIRECT` | Real patient samples or validated disease model | Primary catalog |
| `ADJACENT` | Mechanistically relevant, no patient data | Include, flagged |
| `WEAK` | Epidemiologically linked / tangential | Include, flagged; user decides |
| `REMOVE` | False positive, no disease connection | Exclude — **always manually review before deleting** |

Function + example keyword lists in `references/classification.md`.

### Step 15 — Assemble master catalog

Concat all sources (list `df_geo` **before** `df_biostudies` so richer native GEO
records win dedup). Deduplicate across repositories by normalizing accessions:
`E-GEOD-NNN → GSENNN`, so GEO/ArrayExpress mirrors collapse while genuine
ArrayExpress studies (`E-MTAB-`, `E-MEXP-`, `E-TABM-`) survive. **Do NOT
pre-filter `E-GEOD-` accessions at query time** — that silently drops studies GEO
missed. Backfill authoritative organism for `CORE/ADJACENT` ArrayExpress rows from
`/studies/{accession}`. Add `Access` column (`Controlled` for dbGaP/EGA/JGA/UK
Biobank/FinnGen). Sort by relevance tier then year desc. Assembly code in
`references/classification.md`.

### Step 16 — Markdown summary report

Write `<disease>_omics_summary.md`: overview (totals, repos searched, date);
counts by repository; counts by omics type (CORE/DIRECT only); relevance
breakdown; top datasets (largest N, most recent); controlled-access list with
access instructions; coverage gaps and limitations (which APIs failed / need keys
/ are data-sparse for this topic).

### Step 17 — Overview figure (if requested)

3-panel `<disease>_omics_landscape.png` — Panel A donut (CORE/DIRECT by omics
type), Panel B horizontal bar (all by repository), Panel C stacked-bar timeline
(datasets/year by omics type). After saving, verify render with
`Read(file_path=..., mode="media_output_check")`. Plot boilerplate in
`references/query-code.md`.

## Key scientific caveats

- **GEO `gdstype` is unreliable** — reclassify from Title + Summary. ChIP/ATAC/
  Ribo-seq are often mislabeled; some RNA-seq mislabeled "Hi-C".
- **Cross-repository duplicates are common** — same study in GEO + ArrayExpress +
  SRA + OmicsDI. Normalize `E-GEOD-` → `GSE` before dedup (Step 15).
- **Relevance audit is keyword-based** — produces false positives and negatives.
  Always manually review `REMOVE` candidates before excluding.
- **Controlled-access** (dbGaP/EGA/JGA) — catalog and label `Controlled`; download
  needs institutional data-access agreements.
- **Metabolomics is data-sparse** for rare diseases — check Metabolomics Workbench
  + HMDB disease pages manually if OmicsDI returns nothing.
- **CELLxGENE uses MONDO ontology** — look up the correct MONDO ID before
  filtering cell metadata.
- **NCBI rate limits** — 3 req/sec without key, 10 with. `time.sleep(0.34)`
  between calls; set `api_key` param to raise the limit.
