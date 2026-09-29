# AGENTS.md

Operating guidelines and project rules for AI agents working in `phylo3d-trait`.

## 1. Mandatory documentation preflight
Before running, testing, or modifying Phylo3D-Trait, read:
- `README.md`
- `docs/PHYLO3D_TRAIT_USAGE_GUIDE.md`

For interface details, the current CLI and source code are authoritative. If documentation disagrees with `python -m phylo3d_trait.cli --help` or the implementation, report the discrepancy instead of guessing.

## 2. CLI-first principle
For a new tree or trait dataset, use the existing CLI. New data do not imply code changes.

## 3. Deterministic clade IDs
Generate internal ancestral node IDs with `template-values`. Never guess or hand-code internal-node IDs.

## 4. Time and branch-length verification
Do not interpret the Z axis as time before present unless branch lengths represent evolutionary time and the tree is appropriate for that interpretation.

## 5. No built-in ancestral-state reconstruction
Phylo3D-Trait visualizes supplied node values. It does not infer ancestral states. All required tip, internal-node, and root values must be supplied explicitly.

## 6. Renderer choice
- `--renderer four-layer`: recommended for publication-oriented output and transparent curtain meshes; WebGL2 depth peeling, opacity 0.5–1.0.
- `--renderer plotly`: legacy/general exploration backend.
Do not create another renderer unless the existing architecture cannot support a documented requirement.

## 7. Standard execution flow
`inspect input -> verify branch-length meaning -> template-values -> map traits -> validate -> plot -> verify HTML`

## 8. Code modification invariant
Before changing source code:
`Understand -> Search -> Reuse > Extend > Refactor > Create -> Test -> Re-search for duplicate logic`

Do not hard-code taxa, datasets, node numbers, or project-specific scientific values into package code.
