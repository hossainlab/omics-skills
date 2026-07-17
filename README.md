# omics-skills

> A curated collection of [Claude](https://claude.com/claude-code) skills for omics data analysis.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Skills](https://img.shields.io/badge/skills-1%2B-blue.svg)](#available-skills)
[![Status](https://img.shields.io/badge/status-active%20development-green.svg)](#roadmap)

## Overview

**omics-skills** packages reusable, domain-specific workflows for high-throughput biology as [Agent Skills](https://docs.claude.com/en/docs/claude-code/skills) — modular capabilities that Claude loads on demand to perform reproducible omics analyses.

Each skill bundles the instructions, scripts, and references needed to carry a task end to end: from raw reads and count matrices to quality control, statistics, and publication-ready figures. Skills are added incrementally, one domain at a time, and are designed to compose with one another and with the broader scientific-computing ecosystem.

Covered (and planned) domains include:

- **Genomics** — variant calling, annotation, alignment QC
- **Transcriptomics** — bulk and single-cell RNA-seq, differential expression
- **Proteomics** — mass-spec quantification, protein identification
- **Metabolomics** — feature detection, pathway enrichment
- **Epigenomics** — ChIP-seq, ATAC-seq, methylation
- **Multi-omics** — data integration and joint analysis

## Available skills

| Skill | Domain | Description | Status |
|-------|--------|-------------|--------|
| [omics-dataset-retrieval](skills/omics-dataset-retrieval) | Data retrieval | Retrieve, deduplicate, classify, and relevance-audit public omics datasets for a disease, gene, or process across 25+ repositories (GEO, SRA, ArrayExpress, PRIDE, OmicsDI, CELLxGENE, GDC/TCGA, ENCODE, …) | ✅ |

> This table is updated as each skill lands. See the [Roadmap](#roadmap) for what's next.

## What is a skill?

A skill is a self-contained directory that Claude discovers and loads when a task matches its purpose. The minimal layout:

```
skill-name/
├── SKILL.md          # Metadata (name, description) + instructions
├── scripts/          # Optional executable helpers
├── references/       # Optional docs loaded on demand
└── assets/           # Optional templates, schemas, example data
```

The `SKILL.md` frontmatter tells Claude *when* to use the skill; the body tells it *how*. Claude reads only what it needs, keeping context lean.

## Installation

Clone the repository into your Claude skills directory.

**Project-level** (available in one project):

```bash
git clone https://github.com/hossainlab/omics-skills.git .claude/skills/omics-skills
```

**User-level** (available everywhere):

```bash
# macOS / Linux
git clone https://github.com/hossainlab/omics-skills.git ~/.claude/skills/omics-skills

# Windows (PowerShell)
git clone https://github.com/hossainlab/omics-skills.git $env:USERPROFILE\.claude\skills\omics-skills
```

Claude Code auto-discovers skills placed under a `skills/` directory. Restart your session or run `/skills` to confirm they are loaded.

## Usage

Skills activate automatically when your request matches a skill's description. Just describe the task in natural language:

```
> Run differential expression on these RNA-seq counts and give me a volcano plot
```

You can also invoke a skill by name:

```
> /omics:<skill-name> <arguments>
```

## Requirements

- [Claude Code](https://claude.com/claude-code) (or another Agent Skills-compatible client)
- Python 3.10+ for skills that ship analysis scripts
- Domain tools are documented per skill (e.g. `samtools`, `scanpy`, `DESeq2`)

## Roadmap

- [x] Omics dataset retrieval — catalog public datasets across 25+ repositories
- [ ] Transcriptomics — bulk RNA-seq differential expression
- [ ] Single-cell RNA-seq QC and clustering
- [ ] Variant annotation and filtering
- [ ] Pathway and enrichment analysis
- [ ] Multi-omics integration

## Contributing

Contributions are welcome. To propose a new skill:

1. Fork the repository and create a branch.
2. Add a directory with a well-formed `SKILL.md` (clear `name` and `description` frontmatter).
3. Keep scripts deterministic, documented, and dependency-light.
4. Open a pull request describing the workflow and its inputs/outputs.

## License

Released under the [MIT License](LICENSE). © 2026 Jubayer Hossain.
