# Classification, Relevance Audit, and Catalog Assembly

## Step 13 — Classify omics type

The GEO `gdstype` field is unreliable — always use Title + Summary as the primary
signal. Check single-cell/spatial **before** bulk RNA-seq (order matters).

```python
def classify_omics(row):
    combined = (str(row.get("GEO_Type","")) + " " +
                str(row.get("Title","")) + " " +
                str(row.get("Summary",""))).lower()

    if any(x in combined for x in ["single cell", "scrna", "10x chromium", "dropseq",
                                    "smart-seq", "single-nucleus", "snrna"]):
        return "scRNA-seq"
    elif any(x in combined for x in ["spatial transcriptom", "visium", "slide-seq",
                                      "merfish", "seqfish", "stereo-seq"]):
        return "Spatial Transcriptomics"
    elif "atac" in combined:
        return "ATAC-seq"
    elif any(x in combined for x in ["cut&run", "cut and run", "cutana", "cut&tag"]):
        return "CUT&RUN/CUT&TAG"
    elif any(x in combined for x in ["chip-seq", "chip seq", "binding/occupancy",
                                      "chromatin immunoprecipitation"]):
        return "ChIP-seq"
    elif any(x in combined for x in ["hi-c", "hic", "3d genome", "chromatin conformation",
                                      "capture-c", "4c-seq", "5c"]):
        return "Hi-C / 3D Genome"
    elif any(x in combined for x in ["methylat", "bisulfite", "wgbs", "rrbs",
                                      "epic array", "850k", "450k", "dnmt"]):
        return "DNA Methylation"
    elif any(x in combined for x in ["mirna", "microrna", "ncrna", "lncrna",
                                      "small rna", "pirna", "circrna"]):
        return "miRNA/ncRNA"
    elif any(x in combined for x in ["ribo-seq", "ribosome profiling", "rip-seq",
                                      "clip-seq", "iclip"]):
        return "Other (Ribo/RIP/CLIP-seq)"
    elif any(x in combined for x in ["high throughput sequencing", "rna-seq", "rnaseq",
                                      "mrna-seq", "total rna", "poly-a"]):
        return "Bulk RNA-seq"
    elif any(x in combined for x in ["expression profiling by array",
                                      "transcription profiling by array", "microarray",
                                      "affymetrix", "illumina beadchip", "agilent"]):
        return "Microarray"
    elif any(x in combined for x in ["wgs", "whole genome sequencing", "whole-genome",
                                      "snp array", "genotyping", "gwas", "exome",
                                      "comparative genomic hybridization"]):
        return "Genomics (WGS/SNP/Exome)"
    elif any(x in combined for x in ["proteom", "mass spectrometry", "lc-ms", "tmtpro",
                                      "label-free", "dia-ms", "dda-ms", "phosphoproteom"]):
        return "Proteomics (MS)"
    elif any(x in combined for x in ["metabolom", "metabolite", "nmr", "gcms", "lcms",
                                      "lipidom", "lipidome", "untargeted ms"]):
        return "Metabolomics/Lipidomics"
    elif any(x in combined for x in ["multi-omics", "multiomics", "multi omics",
                                      "integrat", "joint profil"]):
        return "Multi-omics"
    else:
        return "Other/Mixed"
```

## Step 14 — Relevance audit (4-tier)

Tune keyword lists per disease. The example below is for sickle cell disease (SCD).

```python
def audit_relevance(row, disease_keywords, weak_keywords, exclude_keywords):
    """
    disease_keywords: confirm direct relevance (patient data or validated model)
    weak_keywords:    adjacent/tangential relevance
    exclude_keywords: NOT about the disease (false positives)
    """
    text = (str(row.get("Title","")) + " " + str(row.get("Summary",""))).lower()

    if any(kw.lower() in text for kw in exclude_keywords):   # hard exclusion first
        return "REMOVE"

    strong_hits = sum(1 for kw in disease_keywords if kw.lower() in text)
    weak_hits   = sum(1 for kw in weak_keywords    if kw.lower() in text)

    if strong_hits >= 2:                       return "CORE/DIRECT"
    elif strong_hits == 1 and weak_hits >= 1:  return "CORE/DIRECT"
    elif strong_hits == 1 and weak_hits == 0:  return "ADJACENT"
    elif weak_hits >= 1 and strong_hits == 0:  return "WEAK"
    else:                                      return "WEAK"
```

Tier meanings:

| Tier | Meaning | Action |
|------|---------|--------|
| `CORE/DIRECT` | Actual patient samples or validated model (e.g. HbSS mice) | Primary catalog |
| `ADJACENT` | Mechanistically relevant (e.g. HbF reactivation, globin switching), no patient data | Include, flagged |
| `WEAK` | Epidemiologically linked (e.g. sickle trait, RMC) or tangential | Include, flagged; user decides |
| `REMOVE` | False positive | Exclude — **always manually review first** |

Example keyword lists (SCD — retune per run):

```python
disease_keywords = ["sickle cell", "sickle-cell", "HbSS", "hemoglobin S",
                    "vaso-occlus", "sickling", "SCA patient", "SCD patient",
                    "sickle cell anemia", "sickle cell disease"]
weak_keywords    = ["sickle cell trait", "HbAS", "renal medullary carcinoma",
                    "RMC", "beta-thalassemia", "thalassemia only"]
exclude_keywords = ["lung cancer", "gastric cancer", "colon cancer",
                    "BACH1 lung", "unrelated disease"]  # tune per run
```

Keyword matching produces false negatives — always manually review `REMOVE`
candidates before deleting.

## Step 15 — Assemble master catalog

```python
import pandas as pd, re, requests, time

# Order matters: df_geo BEFORE df_biostudies so the richer native GEO record wins
# dedup when a study exists as both GSE41575 and its E-GEOD-41575 mirror.
all_dfs = [df for df in [df_geo, df_sra, df_biostudies, df_pride, df_jpost,
                          df_metabolomics, df_cxg, df_encode, df_atlas,
                          df_zenodo, df_figshare, df_controlled]
           if df is not None and len(df) > 0]
df_all = pd.concat(all_dfs, ignore_index=True)

# Dedup across repositories. Normalize E-GEOD-NNN -> GSENNN so GEO/ArrayExpress
# mirrors collapse, while unique ArrayExpress studies (E-MTAB-, E-MEXP-, E-TABM-)
# are retained.
def dedup_key(acc):
    acc = str(acc)
    m = re.match(r"^E-GEOD-(\d+)$", acc)
    return f"GSE{m.group(1)}" if m else acc
df_all["_dedup"] = df_all["Accession"].apply(dedup_key)
df_all = df_all.drop_duplicates(subset="_dedup").drop(columns=["_dedup"])

df_all["Omics_Type"] = df_all.apply(classify_omics, axis=1)
df_all["Relevance"]  = df_all.apply(
    lambda row: audit_relevance(row, disease_keywords, weak_keywords, exclude_keywords), axis=1)

# --- Organism accuracy pass for ArrayExpress/BioStudies ---
# Scoped to CORE/DIRECT + ADJACENT so this stays cheap (a few calls, not one per row).
def fetch_biostudies_organism(accession):
    """'; '-joined organism(s) from the per-study record, walking nested subsections.
    '' on any failure."""
    try:
        r = requests.get(f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}",
                         headers={"Accept": "application/json"}, timeout=30)
        if r.status_code != 200:
            return ""
        organisms, stack = [], [r.json().get("section", {})]
        while stack:
            node = stack.pop()
            if isinstance(node, dict):
                for a in node.get("attributes", []) or []:
                    if str(a.get("name", "")).strip().lower() == "organism" and a.get("value"):
                        organisms.append(a["value"].strip())
                if node.get("subsections"):
                    stack.append(node["subsections"])
            elif isinstance(node, list):
                stack.extend(node)
        seen, out = set(), []
        for o in organisms:
            if o.lower() not in seen:
                seen.add(o.lower()); out.append(o)
        return "; ".join(out)
    except Exception:
        return ""

needs_organism = df_all[
    (df_all["Repository"] == "ArrayExpress/BioStudies")
    & (df_all["Relevance"].isin(["CORE/DIRECT", "ADJACENT"]))
    & (df_all["Organism"].fillna("").str.strip() == "")
]
for idx, row in needs_organism.iterrows():
    org = fetch_biostudies_organism(row["Accession"])
    if org:
        df_all.at[idx, "Organism"] = org
    time.sleep(0.2)

# Access level
CONTROLLED_REPOS = {"dbGaP", "EGA", "JGA", "UK Biobank", "FinnGen"}
df_all["Access"] = df_all["Repository"].apply(
    lambda r: "Controlled" if any(c in r for c in CONTROLLED_REPOS) else "Open")

# Sort: relevance tier, then year descending
tier_order = {"CORE/DIRECT": 0, "ADJACENT": 1, "WEAK": 2, "REMOVE": 3}
df_all["_tier"] = df_all["Relevance"].map(tier_order)
df_all = df_all.sort_values(["_tier", "Date"], ascending=[True, False]).drop(columns=["_tier"])

# Save
df_all.to_csv(f"{output_dir}/{disease_slug}_omics_datasets_MASTER.csv", index=False)
df_validated = df_all[df_all["Relevance"].isin(["CORE/DIRECT", "ADJACENT"])]
df_validated.to_csv(f"{output_dir}/{disease_slug}_omics_datasets_VALIDATED.csv", index=False)

print(f"Total datasets: {len(df_all)}")
print(f"CORE/DIRECT: {(df_all.Relevance=='CORE/DIRECT').sum()}")
print(f"ADJACENT: {(df_all.Relevance=='ADJACENT').sum()}")
print(f"WEAK: {(df_all.Relevance=='WEAK').sum()}")
print(f"REMOVE: {(df_all.Relevance=='REMOVE').sum()}")
```
