# Example Runs

## "Find all publicly available omics datasets for Alzheimer's disease"

- **Synonyms:** `["AD", "Alzheimer", "LOAD", "EOAD", "dementia", "amyloid", "tau"]`
- **Tissue terms:** brain, cortex, hippocampus, CSF, blood, iPSC neurons, microglia, astrocytes
- **Gene terms:** APOE, APP, PSEN1, PSEN2, TREM2, BIN1, CLU, ABCA7
- Run 30+ GEO queries, ArrayExpress, PRIDE, OmicsDI, CELLxGENE, Expression Atlas,
  Zenodo, Figshare. Classify omics types, audit relevance.
- **Output:** `alzheimers_omics_datasets_MASTER.csv`,
  `alzheimers_omics_datasets_VALIDATED.csv`, `alzheimers_omics_summary.md`

## "What RNA-seq data is available for BCL11A in erythroid cells?"

- **Topic:** BCL11A erythroid
- **Synonyms:** `["BCL11A", "CTIP1", "fetal hemoglobin", "HbF", "globin switching"]`
- **Tissue terms:** erythroid, HUDEP, CD34, erythroblast, reticulocyte
- Run targeted GEO + ArrayExpress queries; skip PRIDE/metabolomics (RNA-seq focus).
- **Output:** `bcl11a_erythroid_omics_datasets_MASTER.csv`

## "Survey all cancer proteomics datasets for pancreatic ductal adenocarcinoma"

- **Synonyms:** `["PDAC", "pancreatic cancer", "pancreatic adenocarcinoma"]`
- Query PRIDE, jPOST, iProX, GDC/TCGA (PAAD project), CPTAC, cBioPortal.
- Skip metabolomics-only repos; focus on MS proteomics.
- **Output:** `pdac_omics_datasets_MASTER.csv`
