"""Pure, offline-testable helpers for the omics-dataset-retrieval skill.

These are the deterministic building blocks used by the workflow in ../SKILL.md:
omics-type classification, relevance auditing, and cross-repository accession
deduplication. Network query code stays in references/query-code.md — only the
side-effect-free logic lives here so it can be unit tested.
"""

from __future__ import annotations

import re
from typing import Iterable, Mapping

__all__ = ["classify_omics", "audit_relevance", "dedup_key"]


def _text(row: Mapping[str, object], *fields: str) -> str:
    """Lowercased concatenation of the named row fields (missing -> '')."""
    return " ".join(str(row.get(f, "") or "") for f in fields).lower()


def classify_omics(row: Mapping[str, object]) -> str:
    """Classify a dataset record into an omics type from Title + Summary + GEO_Type.

    The GEO ``gdstype`` field is unreliable, so Title + Summary are the primary
    signal. Order matters: single-cell / spatial are checked before bulk RNA-seq,
    ATAC before ChIP, so the most specific label wins.
    """
    combined = _text(row, "GEO_Type", "Title", "Summary")

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


def audit_relevance(row: Mapping[str, object],
                    disease_keywords: Iterable[str],
                    weak_keywords: Iterable[str],
                    exclude_keywords: Iterable[str]) -> str:
    """Score a record into one of four relevance tiers.

    ``disease_keywords`` confirm direct relevance (patient data / validated model),
    ``weak_keywords`` mark adjacent/tangential relevance, and ``exclude_keywords``
    force a REMOVE. Exclusion is checked first so a false-positive hit can never be
    promoted to CORE. Returns one of: CORE/DIRECT, ADJACENT, WEAK, REMOVE.
    """
    text = _text(row, "Title", "Summary")

    if any(kw.lower() in text for kw in exclude_keywords):
        return "REMOVE"

    strong_hits = sum(1 for kw in disease_keywords if kw.lower() in text)
    weak_hits = sum(1 for kw in weak_keywords if kw.lower() in text)

    if strong_hits >= 2:
        return "CORE/DIRECT"
    elif strong_hits == 1 and weak_hits >= 1:
        return "CORE/DIRECT"
    elif strong_hits == 1 and weak_hits == 0:
        return "ADJACENT"
    else:  # weak_hits >= 1 and strong_hits == 0, or no hits at all
        return "WEAK"


def dedup_key(accession: object) -> str:
    """Normalize an accession so GEO/ArrayExpress mirrors collapse on dedup.

    ArrayExpress mirrors a GEO series ``GSE41575`` as ``E-GEOD-41575``; both refer
    to one study. Map ``E-GEOD-NNN -> GSENNN`` so they dedupe together, while native
    ArrayExpress studies (E-MTAB-, E-MEXP-, E-TABM-, …) keep their own key.
    """
    acc = str(accession)
    m = re.match(r"^E-GEOD-(\d+)$", acc)
    return f"GSE{m.group(1)}" if m else acc
