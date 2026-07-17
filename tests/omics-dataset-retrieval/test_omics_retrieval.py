"""Unit tests for the omics-dataset-retrieval skill's pure helpers.

Covers the deterministic logic that drives catalog quality: omics-type
classification, 4-tier relevance auditing, and cross-repository accession
deduplication. Network query code is not exercised here (see references/).
"""

import importlib.util
from pathlib import Path

import pytest

# Load the skill's helper module by path (skill dir name has a hyphen, so it is
# not importable as a normal package).
SKILL_DIR = Path(__file__).resolve().parents[2] / "skills" / "omics-dataset-retrieval"
MODULE_PATH = SKILL_DIR / "scripts" / "omics_retrieval.py"

_spec = importlib.util.spec_from_file_location("omics_retrieval", MODULE_PATH)
omics_retrieval = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(omics_retrieval)

classify_omics = omics_retrieval.classify_omics
audit_relevance = omics_retrieval.audit_relevance
dedup_key = omics_retrieval.dedup_key


# --------------------------------------------------------------------------- #
# classify_omics
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("title, summary, expected", [
    ("Single cell RNA-seq of cortex", "10x chromium", "scRNA-seq"),
    ("Spatial map of tumor", "Visium slides", "Spatial Transcriptomics"),
    ("Chromatin accessibility", "ATAC-seq of HSCs", "ATAC-seq"),
    ("TF occupancy", "ChIP-seq for GATA1", "ChIP-seq"),
    ("Genome folding", "Hi-C in erythroid cells", "Hi-C / 3D Genome"),
    ("Methylome", "whole genome bisulfite sequencing (WGBS)", "DNA Methylation"),
    ("Small RNA study", "microRNA profiling", "miRNA/ncRNA"),
    ("Translation", "ribosome profiling (Ribo-seq)", "Other (Ribo/RIP/CLIP-seq)"),
    ("Expression", "bulk RNA-seq of blood", "Bulk RNA-seq"),
    ("Array study", "expression profiling by array, Affymetrix", "Microarray"),
    ("Variant discovery", "whole genome sequencing WGS", "Genomics (WGS/SNP/Exome)"),
    ("Protein levels", "mass spectrometry proteomics", "Proteomics (MS)"),
    ("Small molecules", "untargeted metabolomics by NMR", "Metabolomics/Lipidomics"),
    ("Integrated study", "multi-omics integration", "Multi-omics"),
    ("Mystery dataset", "no recognizable assay terms", "Other/Mixed"),
])
def test_classify_omics_basic(title, summary, expected):
    assert classify_omics({"Title": title, "Summary": summary}) == expected


def test_classify_prefers_single_cell_over_bulk():
    """A record mentioning both single-cell and RNA-seq must resolve to scRNA-seq."""
    row = {"Title": "single cell RNA-seq atlas", "Summary": "rna-seq"}
    assert classify_omics(row) == "scRNA-seq"


def test_classify_atac_before_chip():
    row = {"Title": "ATAC-seq and ChIP-seq combined", "Summary": ""}
    assert classify_omics(row) == "ATAC-seq"


def test_classify_lcms_collides_to_proteomics():
    """Known ordering caveat: 'lc-ms' is a proteomics keyword checked before
    metabolomics, so a metabolomics study described only as 'LC-MS' mislabels as
    Proteomics. Reclassify or add context terms (NMR, GC-MS, 'metabolite') to
    disambiguate."""
    row = {"Title": "untargeted metabolomics", "Summary": "LC-MS platform"}
    assert classify_omics(row) == "Proteomics (MS)"


def test_classify_uses_geo_type_field():
    row = {"GEO_Type": "Expression profiling by array", "Title": "", "Summary": ""}
    assert classify_omics(row) == "Microarray"


def test_classify_handles_missing_and_none_fields():
    assert classify_omics({}) == "Other/Mixed"
    assert classify_omics({"Title": None, "Summary": None}) == "Other/Mixed"


# --------------------------------------------------------------------------- #
# audit_relevance
# --------------------------------------------------------------------------- #

DISEASE_KW = ["sickle cell", "hbss", "sickle cell disease"]
WEAK_KW = ["sickle cell trait", "beta-thalassemia"]
EXCLUDE_KW = ["lung cancer", "unrelated disease"]


def _audit(title, summary=""):
    return audit_relevance({"Title": title, "Summary": summary},
                           DISEASE_KW, WEAK_KW, EXCLUDE_KW)


def test_audit_two_strong_hits_is_core():
    assert _audit("sickle cell disease patient HbSS samples") == "CORE/DIRECT"


def test_audit_one_strong_plus_weak_is_core():
    assert _audit("sickle cell study", "beta-thalassemia comparison") == "CORE/DIRECT"


def test_audit_one_strong_only_is_adjacent():
    assert _audit("sickle cell mechanism in cell line") == "ADJACENT"


def test_audit_weak_only_is_weak():
    assert _audit("beta-thalassemia patients only") == "WEAK"


def test_audit_no_hits_is_weak():
    assert _audit("generic erythroid study") == "WEAK"


def test_audit_substring_collision_promotes_trait_to_core():
    """Known keyword caveat: 'sickle cell trait' contains the strong keyword
    'sickle cell', so it scores 1 strong + 1 weak = CORE/DIRECT even though a
    trait-carrier study is really WEAK. This is why REMOVE/CORE candidates need
    manual review. Use word-boundary or exclusion terms to harden per run."""
    assert _audit("sickle cell trait carriers") == "CORE/DIRECT"


def test_audit_exclusion_wins_over_strong_hits():
    """Exclusion is checked first: a strong disease hit cannot rescue an excluded row."""
    row = "sickle cell disease HbSS but in a lung cancer cohort"
    assert _audit(row) == "REMOVE"


def test_audit_is_case_insensitive():
    assert _audit("SICKLE CELL DISEASE HBSS") == "CORE/DIRECT"


# --------------------------------------------------------------------------- #
# dedup_key
# --------------------------------------------------------------------------- #

def test_dedup_maps_arrayexpress_geo_mirror_to_gse():
    assert dedup_key("E-GEOD-41575") == "GSE41575"


def test_dedup_collapses_mirror_and_native_gse():
    assert dedup_key("E-GEOD-41575") == dedup_key("GSE41575")


def test_dedup_preserves_native_arrayexpress_accessions():
    for acc in ["E-MTAB-1234", "E-MEXP-999", "E-TABM-42"]:
        assert dedup_key(acc) == acc


def test_dedup_preserves_unrelated_accessions():
    for acc in ["GSE100", "PXD000123", "SRP0456", "MTBLS77"]:
        assert dedup_key(acc) == acc


def test_dedup_accepts_non_string_input():
    assert dedup_key(12345) == "12345"
