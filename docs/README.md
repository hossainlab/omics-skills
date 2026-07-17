# Documentation

Human-facing documentation for the skills in this repository. Each page explains
what a skill does, its inputs and outputs, how to invoke it, and worked examples.

For the agent-facing instructions Claude actually loads, see each skill's
`SKILL.md` under `skills/<name>/`.

## Skills

| Skill | Domain | Documentation |
|-------|--------|---------------|
| omics-dataset-retrieval | Data retrieval | [docs](omics-dataset-retrieval.md) |
| bulk-rnaseq-counts-to-de-deseq2 | Transcriptomics | [docs](bulk-rnaseq-counts-to-de-deseq2.md) |

## Conventions

- **SKILL.md** — agent instructions (frontmatter `name` + `description`, loaded on
  demand by Claude).
- **references/** — deep-dive material loaded only when a task needs it.
- **scripts/** — runnable helpers.
- **docs/** — this folder; prose documentation for humans.
