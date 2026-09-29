# Contributing to Phylo3D-Trait

Bug reports, feature requests, documentation improvements, tests, and pull requests are welcome.

## Bug reports

Use the GitHub bug-report template. A useful rendering/parsing report should include:
- Python version and operating system;
- browser/version when the Four-Layer WebGL2 renderer is involved;
- exact command;
- renderer and relevant options;
- minimal tree file;
- minimal node-value CSV/TSV;
- observed behavior and expected behavior;
- screenshot or generated HTML when it materially helps diagnosis.

Remove unpublished or sensitive biological data from reproducing examples whenever possible.

## Feature requests

Describe the scientific or visualization problem first. If an existing CLI option, renderer, parser, or data model can solve it, prefer extending that implementation instead of creating a parallel subsystem.

## Pull requests

Before coding:
1. read `README.md`, `AGENTS.md`, and `docs/PHYLO3D_TRAIT_USAGE_GUIDE.md`;
2. search the repository for existing implementations/helpers/tests;
3. prefer reuse > extension > refactor > new implementation;
4. do not hard-code dataset-specific taxa, node IDs, values, or topology.

Before opening a PR:

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

Include:
- what problem the PR solves;
- why the chosen implementation fits the existing architecture;
- tests added/updated;
- any user-visible CLI or documentation changes.

Small, reviewable PRs are preferred. Scientific-coordinate invariants and existing input semantics must not be changed merely to make a figure look better.

## Scope

Phylo3D-Trait visualizes supplied phylogenies and continuous node values. Tree inference, tree dating, ASR, sequence analysis, and evolutionary model fitting are outside package scope unless the project scope is explicitly changed.
