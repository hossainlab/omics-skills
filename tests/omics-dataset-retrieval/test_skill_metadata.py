"""Structural checks on the omics-dataset-retrieval skill package.

Guards the contract that makes the skill discoverable: valid SKILL.md frontmatter
and the presence of the reference files the workflow points at.
"""

from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[2] / "skills" / "omics-dataset-retrieval"
SKILL_MD = SKILL_DIR / "SKILL.md"


def _parse_frontmatter(text):
    """Minimal YAML-frontmatter extractor (name/description only, no PyYAML dep)."""
    assert text.startswith("---\n"), "SKILL.md must open with a --- frontmatter fence"
    _, fm, _body = text.split("---\n", 2)
    meta, key = {}, None
    for line in fm.splitlines():
        if not line.strip():
            continue
        if line[0] not in " \t" and ":" in line:          # top-level key
            key, _, val = line.partition(":")
            key = key.strip()
            meta[key] = val.strip().lstrip(">|").strip()
        elif key:                                          # folded continuation
            meta[key] = (meta[key] + " " + line.strip()).strip()
    return meta


@pytest.fixture(scope="module")
def frontmatter():
    return _parse_frontmatter(SKILL_MD.read_text(encoding="utf-8"))


def test_skill_md_exists():
    assert SKILL_MD.is_file()


def test_name_matches_directory(frontmatter):
    assert frontmatter.get("name") == "omics-dataset-retrieval"


def test_description_present_and_substantial(frontmatter):
    desc = frontmatter.get("description", "")
    assert len(desc) >= 50, "description should meaningfully say when to use the skill"


def test_description_states_boundary(frontmatter):
    """Description must flag that the skill does not download or analyze data."""
    desc = frontmatter.get("description", "").lower()
    assert "not" in desc and ("download" in desc or "analysis" in desc)


@pytest.mark.parametrize("ref", [
    "references/repositories.md",
    "references/query-code.md",
    "references/classification.md",
    "references/examples.md",
    "scripts/omics_retrieval.py",
])
def test_referenced_files_exist(ref):
    assert (SKILL_DIR / ref).is_file(), f"missing {ref}"
