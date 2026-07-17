# Per-Repository Query Code

Runnable snippets for each repository. Shared variables assumed in scope:
`disease` (str), `synonyms` (list[str]), `output_dir` (path), `disease_slug` (str).

Common preamble:

```python
import requests, time, json, re
import pandas as pd

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
```

**Organism policy for every source:** keep all organisms by default. Apply an
organism filter only if the user explicitly requested a specific species.

---

## Step 3 — GEO (NCBI E-utilities) — PRIMARY SOURCE

```python
# Build 20-40 queries combining disease + omics + tissue + gene terms
search_queries = {
    "title_rnaseq":    f'"{disease}"[Title] AND "RNA-seq"',
    "title_scrna":     f'"{disease}"[Title] AND "single cell"',
    "title_array":     f'"{disease}"[Title] AND "expression profiling by array"[DataSet Type]',
    "title_chip":      f'"{disease}" AND "ChIP-seq"',
    "title_atac":      f'"{disease}" AND "ATAC-seq"',
    "title_methyl":    f'"{disease}" AND "methylation"',
    "title_wgs":       f'"{disease}" AND "whole genome sequencing"',
    "title_proteom":   f'"{disease}" AND "proteomics"',
    "title_metabolom": f'"{disease}" AND "metabolomics"',
    # + f'"{synonym}"[Title] AND "RNA-seq"' for each synonym
    # + f'"{disease}" AND "{tissue_term}"' for key tissues
    # + f'"{disease}" AND "{key_gene}"' for key disease genes
}

all_ids = set()
for label, query in search_queries.items():
    r = requests.get(f"{BASE}esearch.fcgi",
                     params={"db": "gds", "term": query, "retmax": 100, "retmode": "json"})
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    all_ids.update(ids)
    time.sleep(0.34)  # NCBI rate limit (3 req/sec without key)

records = []
for i in range(0, len(list(all_ids)), 20):
    batch = list(all_ids)[i:i+20]
    r = requests.get(f"{BASE}esummary.fcgi",
                     params={"db": "gds", "id": ",".join(batch), "retmode": "json"})
    result = r.json().get("result", {})
    for uid in result.get("uids", []):
        item = result[uid]
        records.append({
            "Accession": item.get("accession", ""),
            "Title": item.get("title", ""),
            "GEO_Type": item.get("gdstype", ""),
            "Organism": item.get("taxon", ""),
            "N_Samples": item.get("n_samples", ""),
            "Date": item.get("pdat", ""),
            "Summary": item.get("summary", "")[:500],
            "Repository": "GEO",
        })
    time.sleep(0.34)

df_geo = pd.DataFrame(records)
# Keep GSE/GDS only (exclude GPL platform records). Keep ALL organisms by default;
# organism filter only if user requested a species, e.g.:
#   df_geo = df_geo[df_geo["Organism"].str.contains("Homo sapiens", na=False)]
df_geo = df_geo[df_geo["Accession"].str.startswith(("GSE", "GDS"))]
```

---

## Step 4 — SRA (NCBI) — raw sequencing studies not in GEO

```python
sra_records = []
for query in [disease] + synonyms[:3]:
    r = requests.get(f"{BASE}esearch.fcgi",
                     params={"db": "sra", "term": f'"{query}"',
                             "retmax": 100, "retmode": "json"})
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if ids:
        r2 = requests.get(f"{BASE}esummary.fcgi",
                          params={"db": "sra", "id": ",".join(ids), "retmode": "json"})
        result = r2.json().get("result", {})
        for uid in result.get("uids", []):
            item = result[uid]
            exp_acc = item.get("experiment", {}).get("acc", "")
            sra_records.append({
                "Accession": item.get("runs", {}).get("Run", [{}])[0].get("acc", exp_acc),
                "Title": item.get("title", ""),
                "GEO_Type": item.get("exptype", ""),
                "Organism": item.get("organism", {}).get("ScientificName", "unknown"),
                "N_Samples": item.get("runs", {}).get("@total", ""),
                "Date": item.get("createdate", ""),
                "Summary": item.get("summary", "")[:500],
                "Repository": "SRA",
            })
    time.sleep(0.34)
df_sra = pd.DataFrame(sra_records)
```

---

## Step 5 — ArrayExpress / BioStudies (EBI)

Search hits expose `accession`, `title`, `type` (always `"study"`),
`release_date`, `files`, `links`, and a free-text `content` blob (assay type, EFO
terms, organism, description). There is NO releaseDate/description/organism/
numberOfSamples field on the hit — those live at `/studies/{accession}` (too
expensive per row). Read `release_date`; mine `content`.

```python
biostudies_url = "https://www.ebi.ac.uk/biostudies/api/v1/search"
biostudies_records = []

ORGANISM_PATTERNS = ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio",
                     "Drosophila melanogaster", "Saccharomyces cerevisiae",
                     "Caenorhabditis elegans", "Gallus gallus", "Sus scrofa",
                     "Macaca mulatta", "Pan troglodytes", "Arabidopsis thaliana"]

for keyword in [disease] + synonyms[:2]:
    r = requests.get(biostudies_url,
                     params={"query": keyword, "collection": "arrayexpress",
                             "pageSize": 100, "page": 1},
                     headers={"Accept": "application/json"}, timeout=60)
    if r.status_code == 200:
        for hit in r.json().get("hits", []):
            content = hit.get("content", "") or ""
            organism = "; ".join([o for o in ORGANISM_PATTERNS if o in content])
            biostudies_records.append({
                "Accession": hit.get("accession", ""),
                "Title": hit.get("title", ""),
                "GEO_Type": "",  # real assay type is in `content`; classify_omics recovers it
                "Organism": organism,
                "N_Samples": "",
                "Date": hit.get("release_date", ""),
                "Summary": content[:500],
                "Repository": "ArrayExpress/BioStudies",
            })
    time.sleep(1)

df_biostudies = pd.DataFrame(biostudies_records).drop_duplicates(subset="Accession")
# Do NOT drop E-GEOD-* accessions here — many are studies GEO missed. They dedupe
# against native GEO at assembly (Step 15) via E-GEOD-NNN <-> GSENNN normalization.
```

---

## Step 6 — PRIDE / ProteomeXchange (proteomics)

```python
pride_url = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/search"
pride_records = []
for keyword in [disease] + synonyms:
    r = requests.get(pride_url,
                     params={"keyword": keyword, "pageSize": 100, "page": 0},
                     headers={"Accept": "application/json"}, timeout=60)
    if r.status_code == 200:
        data = r.json()
        projects = data if isinstance(data, list) else data.get("_embedded", {}).get("compactprojects", [])
        for p in projects:
            pride_records.append({
                "Accession": p.get("accession", ""),
                "Title": p.get("title", ""),
                "GEO_Type": "Proteomics (MS)",
                "Organism": "; ".join(p.get("organisms", [])),
                "N_Samples": p.get("numberOfSamples", ""),
                "Date": p.get("submissionDate", ""),
                "Summary": p.get("projectDescription", "")[:500],
                "Repository": "PRIDE",
            })
    time.sleep(1)
df_pride = pd.DataFrame(pride_records).drop_duplicates(subset="Accession")
```

jPOST and iProX (ProteomeXchange members not always in PRIDE search):

```python
jpost_records = []
for keyword in [disease] + synonyms[:2]:
    r = requests.get("https://repository.jpostdb.org/search",
                     params={"keyword": keyword, "format": "json"}, timeout=30)
    if r.status_code == 200:
        for p in r.json().get("projects", []):
            jpost_records.append({
                "Accession": p.get("projectId", ""),
                "Title": p.get("title", ""),
                "GEO_Type": "Proteomics (MS)",
                "Organism": p.get("organism", ""),
                "N_Samples": "",
                "Date": p.get("releaseDate", ""),
                "Summary": p.get("description", "")[:500],
                "Repository": "jPOST",
            })
    time.sleep(1)
df_jpost = pd.DataFrame(jpost_records)
```

---

## Step 7 — OmicsDI (metabolomics aggregator)

Best single entry point for metabolomics. Filter by source to avoid GEO dupes.

```python
omicsdi_url = "https://www.omicsdi.org/ws/dataset/search"
omicsdi_records = []
METABOLOMICS_REPOS = {"MetaboLights", "MassIVE", "GNPS", "Metabolon",
                      "Metabolomics Workbench", "HMDB", "Lipidmaps"}

for keyword in [disease] + synonyms[:2]:
    r = requests.get(omicsdi_url,
                     params={"query": keyword, "size": 100, "start": 0},
                     headers={"Accept": "application/json"}, timeout=60)
    if r.status_code == 200:
        for ds in r.json().get("datasets", []):
            repo = ds.get("source", "")
            if repo in METABOLOMICS_REPOS:
                omicsdi_records.append({
                    "Accession": ds.get("id", ""),
                    "Title": ds.get("name", ""),
                    "GEO_Type": "Metabolomics",
                    "Organism": "; ".join(ds.get("organisms", {}).get("name", [])),
                    "N_Samples": "",
                    "Date": ds.get("publicationDate", ""),
                    "Summary": ds.get("description", "")[:500],
                    "Repository": repo,
                })
    time.sleep(1)

# Metabolomics Workbench REST — NIH-funded studies
mw_url = "https://www.metabolomicsworkbench.org/rest/study/study_title"
for keyword in [disease] + synonyms[:2]:
    r = requests.get(f"{mw_url}/{requests.utils.quote(keyword)}/summary/json", timeout=30)
    if r.status_code == 200 and r.text.strip():
        for sid, meta in r.json().items():
            omicsdi_records.append({
                "Accession": sid,
                "Title": meta.get("study_title", ""),
                "GEO_Type": "Metabolomics",
                "Organism": meta.get("subject_species", ""),
                "N_Samples": meta.get("subject_count", ""),
                "Date": meta.get("submit_date", ""),
                "Summary": meta.get("study_summary", "")[:500],
                "Repository": "Metabolomics Workbench",
            })

df_metabolomics = pd.DataFrame(omicsdi_records).drop_duplicates(subset="Accession")
```

---

## Step 8 — CZ CELLxGENE (single-cell)

```python
cxg_url = "https://api.cellxgene.cziscience.com/curation/v1/collections"
r = requests.get(cxg_url, headers={"Accept": "application/json"}, timeout=60)
cxg_records = []
if r.status_code == 200:
    for col in r.json():
        title = col.get("name", "")
        desc = col.get("description", "")
        text = (title + " " + desc).lower()
        if any(kw.lower() in text for kw in [disease] + synonyms):
            cxg_records.append({
                "Accession": col.get("collection_id", ""),
                "Title": title,
                "GEO_Type": "scRNA-seq",
                "Organism": "; ".join(sorted({d.get("organism","") for d in col.get("datasets",[]) if d.get("organism")})) or "unknown",
                "N_Samples": col.get("cell_count", ""),
                "Date": col.get("published_at", ""),
                "Summary": desc[:500],
                "Repository": "CZ CELLxGENE",
            })
df_cxg = pd.DataFrame(cxg_records)

# Cell-level metadata via Census (uses MONDO ontology terms):
# import cellxgene_census
# census = cellxgene_census.open_soma()
# obs = census["census_data"]["homo_sapiens"].obs.read(
#     value_filter=f'disease == "{disease_ontology_term}"').concat().to_pandas()
```

---

## Step 9 — GDC / TCGA (cancer only — skip for non-cancer)

```python
gdc_url = "https://api.gdc.cancer.gov/projects"
r = requests.get(gdc_url,
                 params={"filters": json.dumps({"op": "in", "content": {
                             "field": "disease_type", "value": [disease] + synonyms}}),
                         "fields": "project_id,name,disease_type,primary_site,summary",
                         "format": "json", "size": 100}, timeout=60)
# Parse @graph/hits into records with Repository = "GDC/TCGA"
```

---

## Step 10 — ENCODE (epigenomics — only when regulatory component)

```python
encode_url = "https://www.encodeproject.org/search/"
encode_records = []
r = requests.get(encode_url,
                 params={"searchTerm": disease, "type": "Experiment",
                         "status": "released", "format": "json", "limit": 100},
                 headers={"Accept": "application/json"}, timeout=60)
if r.status_code == 200:
    for exp in r.json().get("@graph", []):
        encode_records.append({
            "Accession": exp.get("accession", ""),
            "Title": exp.get("description", ""),
            "GEO_Type": exp.get("assay_title", ""),
            "Organism": exp.get("organism", {}).get("scientific_name", ""),
            "N_Samples": len(exp.get("replicates", [])),
            "Date": exp.get("date_released", ""),
            "Summary": exp.get("biosample_summary", "")[:500],
            "Repository": "ENCODE",
        })
        time.sleep(0.1)  # 10 req/sec limit
df_encode = pd.DataFrame(encode_records)
```

---

## Step 11 — Expression Atlas (EBI)

```python
atlas_url = "https://www.ebi.ac.uk/gxa/json/experiments"
r = requests.get(atlas_url, timeout=60)  # omit species= to get all organisms
atlas_records = []
if r.status_code == 200:
    for exp in r.json().get("experiments", []):
        text = (exp.get("experimentDescription","") + " " +
                " ".join(exp.get("factors",[])) + " " +
                " ".join(exp.get("experimentalFactors",[]))).lower()
        if any(kw.lower() in text for kw in [disease] + synonyms):
            atlas_records.append({
                "Accession": exp.get("experimentAccession", ""),
                "Title": exp.get("experimentDescription", ""),
                "GEO_Type": exp.get("experimentType", ""),
                "Organism": "; ".join(exp.get("species", [])) or "unknown",
                "N_Samples": exp.get("numberOfAssays", ""),
                "Date": exp.get("lastUpdate", ""),
                "Summary": "",
                "Repository": "Expression Atlas",
            })
df_atlas = pd.DataFrame(atlas_records)
```

---

## Step 12 — Web search for Tier 3 repositories

Use the `WebSearch` tool. Run all of these; extract accession IDs, titles,
descriptions; add rows with the appropriate `Repository` label. Flag
dbGaP/EGA/JGA as `Access = "Controlled"`.

```
# Open
WebSearch: "{disease}" omics dataset site:zenodo.org
WebSearch: "{disease}" RNA-seq proteomics site:figshare.com
WebSearch: "{disease}" omics data site:datadryad.org
WebSearch: "{disease}" omics dataset site:osf.io
WebSearch: "{disease}" omics dataset site:dataverse.harvard.edu
WebSearch: "{disease}" omics dataset site:synapse.org
# Controlled
WebSearch: "{disease}" WGS genomics site:ncbi.nlm.nih.gov/gap
WebSearch: "{disease}" genome sequencing site:ega-archive.org
WebSearch: "{disease}" genomics site:ddbj.nig.ac.jp/jga
# Consortium / specialized
WebSearch: "{disease}" ICGC genomics dataset
WebSearch: "{disease}" CPTAC proteomics dataset
WebSearch: "{disease}" UK Biobank omics
WebSearch: "{disease}" AWS open data omics
```

---

## Step 17 — Overview figure

```python
import matplotlib.pyplot as plt, matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
# Panel A: donut — CORE/DIRECT datasets by omics type
# Panel B: horizontal bar — all datasets by repository
# Panel C: stacked-bar timeline — datasets per year, colored by omics type
fig.savefig(f"{output_dir}/{disease_slug}_omics_landscape.png", dpi=150, bbox_inches="tight")
```

After saving, verify: `Read(file_path=..., mode="media_output_check")`.
