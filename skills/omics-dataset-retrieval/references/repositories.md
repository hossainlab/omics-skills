# Repository Coverage Map

Work through repositories in priority order. Tier 1 have programmatic APIs and are
highest-yield. Tier 2 have APIs but are specialized. Tier 3 need web search or
manual curation. **Always attempt all tiers** — rare diseases often have data only
in Tier 2/3.

## Tier 1 — High-yield, programmatic APIs (always query)

| Repository | Omics Focus | API Base URL | Notes |
|------------|-------------|--------------|-------|
| GEO (NCBI) | All transcriptomics, epigenomics | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/` | Largest source; run 20–40 targeted queries |
| SRA (NCBI) | Raw sequencing reads | `.../eutils/` (`db=sra`) | Complements GEO; finds studies not deposited as GEO series |
| ArrayExpress / BioStudies (EBI) | Transcriptomics, functional genomics | `https://www.ebi.ac.uk/biostudies/api/v1/search` | European mirror of GEO; many unique studies |
| PRIDE / ProteomeXchange | Proteomics (MS) | `https://www.ebi.ac.uk/pride/ws/archive/v2/projects/search` | Primary proteomics repository |
| OmicsDI (aggregator) | Multi-omics aggregator | `https://www.omicsdi.org/ws/dataset/search` | Covers MetaboLights, MassIVE, GNPS, Metabolomics Workbench, ArrayExpress, GEO, PRIDE — catch metabolomics here, avoid re-querying individual repos |
| CZ CELLxGENE | scRNA-seq, spatial | `https://api.cellxgene.cziscience.com/curation/v1/collections` | 33M+ cells; query by disease via `cellxgene_census` package |
| GDC / TCGA (NCI) | Cancer multi-omics | `https://api.gdc.cancer.gov/` | Best for cancer; covers TCGA, TARGET, CGCI, CPTAC-GDC |

## Tier 2 — Specialized APIs (query when relevant to disease/omics type)

| Repository | Omics Focus | API / Access | When to use |
|------------|-------------|--------------|-------------|
| ENCODE | Epigenomics (ChIP/ATAC/RNA-seq) | `https://www.encodeproject.org/search/?format=json` | Regulatory genomics; TF binding; chromatin accessibility |
| Expression Atlas (EBI) | Bulk + single-cell RNA-seq | `https://www.ebi.ac.uk/gxa/json/experiments` | Curated baseline + differential expression |
| Human Cell Atlas (HCA) | scRNA-seq, spatial | `https://service.azul.data.humancellatlas.org/index/projects` | Reference atlases; healthy tissue baselines |
| Metabolomics Workbench (NIH) | Metabolomics, lipidomics | `https://www.metabolomicsworkbench.org/rest/study/study_id/ST/named_json` | NIH-funded metabolomics; REST API |
| MetaboLights (EBI) | Metabolomics | `https://www.ebi.ac.uk/metabolights/ws/studies/` | European metabolomics repository |
| MassIVE / GNPS | Metabolomics, lipidomics | `https://massive.ucsd.edu/ProteoSAFe/datasets.jsp` | MS-based metabolomics; search via OmicsDI |
| jPOST | Proteomics (MS) | `https://repository.jpostdb.org/search` | Japanese proteomics; ProteomeXchange member |
| iProX | Proteomics (MS) | `https://www.iprox.cn/page/project.html` | Chinese proteomics; ProteomeXchange member |
| cBioPortal | Cancer multi-omics | `https://www.cbioportal.org/api/` | Mutation, CNA, expression, methylation |
| EpiRR / IHEC | Epigenomics reference | `https://www.ebi.ac.uk/epirr/api/` | IHEC reference epigenomes; healthy baselines |
| Human Protein Atlas (HPA) | Proteomics, RNA-seq | `https://www.proteinatlas.org/api/` | Tissue/cell-type protein and RNA expression |
| ENA (EBI) | Raw sequencing | `https://www.ebi.ac.uk/ena/portal/api/search` | European mirror of SRA; raw reads |

## Tier 3 — Web search + manual curation (always attempt via WebSearch tool)

| Repository | Omics Focus | Search Strategy |
|------------|-------------|-----------------|
| Zenodo | All | `WebSearch: "{disease}" omics dataset site:zenodo.org` |
| Figshare | All | `WebSearch: "{disease}" RNA-seq proteomics site:figshare.com` |
| Dryad | All | `WebSearch: "{disease}" omics data site:datadryad.org` |
| OSF | All | `WebSearch: "{disease}" omics dataset site:osf.io` |
| Harvard Dataverse | All | `WebSearch: "{disease}" omics dataset site:dataverse.harvard.edu` |
| Synapse (Sage) | Neuro, cancer | `WebSearch: "{disease}" omics dataset site:synapse.org` |
| ICGC Data Portal | Cancer genomics | `WebSearch: "{disease}" ICGC site:dcc.icgc.org` |
| CPTAC (NCI) | Cancer proteomics | `WebSearch: "{disease}" CPTAC proteomics site:proteomics.cancer.gov` |
| AWS Open Data Registry | All | `WebSearch: "{disease}" omics site:registry.opendata.aws` |
| dbGaP (controlled) | Genomics, WGS | `WebSearch: "{disease}" WGS genomics site:ncbi.nlm.nih.gov/gap` |
| EGA (controlled) | Genomics, WGS | `WebSearch: "{disease}" genome sequencing site:ega-archive.org` |
| JGA (controlled) | Genomics | `WebSearch: "{disease}" genomics site:ddbj.nig.ac.jp/jga` |
| UK Biobank | Population genomics | `WebSearch: "{disease}" UK Biobank omics` |
| FinnGen | Population genomics | `WebSearch: "{disease}" FinnGen GWAS` |

Extract accession IDs, titles, descriptions from results. Flag dbGaP/EGA/JGA as
`Access = "Controlled"`.

## Known API Limitations and Workarounds

| Repository | Issue | Workaround |
|------------|-------|------------|
| GEO | `retmax=100` per query cap | Run 20–40 targeted queries; use `[Title]` field tags |
| PRIDE | Response is list or dict by result count | `data if isinstance(data, list) else data.get("_embedded", {}).get("compactprojects", [])` |
| MetaboLights | Direct search API returns 404 | Use OmicsDI aggregator instead |
| Metabolomics Workbench | REST requires exact study-title match | Keyword search via OmicsDI; use REST for known study IDs |
| ArrayExpress/BioStudies | Search hits expose only `release_date` + free-text `content` (no releaseDate/description/organism/numberOfSamples) | Read `release_date`; mine `content` for summary/organism/assay. Organism mining ~64%; Step 15 backfills from `/studies/{accession}` for CORE/ADJACENT → ~100% where it matters |
| ArrayExpress/BioStudies | ~60% of hits are `E-GEOD-*` GEO mirrors | Do NOT drop at query time; dedupe against GEO at assembly via `E-GEOD-NNN → GSENNN` |
| DISGENET | API v7 needs paid key; returns HTML without | Skip; use GEO + PubMed searches |
| GWAS Catalog | `findByDiseaseTrait` returns 0 for many names | Search by EFO ontology ID; browse `https://www.ebi.ac.uk/gwas/` |
| CellxGene | No disease-specific collection tagging for rare diseases | Filter by `disease` field in Census cell metadata |
| ENCODE | Rate limit 10 req/sec | `time.sleep(0.1)` between requests |
| dbGaP / EGA / JGA | No open programmatic search API | Web search + manual curation; flag Controlled |
| Zenodo / Figshare / Dryad | No disease-specific API | WebSearch with `site:` operator |
| Synapse | Login required for some datasets | Web search public projects; note login requirement |
| GDC/TCGA | Only relevant for cancer | Skip for non-cancer diseases |
| jPOST / iProX | API may be unstable | Fall back to ProteomeCentral `http://central.proteomexchange.org` |

## NCBI API key

Set the `api_key` param on E-utilities calls to raise the rate limit from 3 to 10
req/sec (register free at NCBI). Without a key, keep `time.sleep(0.34)` between
requests.
