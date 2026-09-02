# AGENTS.md

Operating guidelines and project rules for AI Agents working in `phylo3d-trait`.

---

## 1. Mandatory Documentation Preflight
Before running, testing, or modifying Phylo3D-Trait, completely read:
- `README.md`
- `docs/PHYLO3D_TRAIT_USAGE_GUIDE.md`

## 2. CLI-First Principle (New Data $\neq$ Code Changes)
For any new tree or trait dataset, run via CLI commands (`python -m phylo3d_trait.cli`). Do not modify package source code for new datasets.

## 3. Deterministic Clade IDs
Internal ancestral node IDs must be generated via `template-values` (derived from sorted descendant tip hashes). Never guess or hand-code internal node IDs.

## 4. Time & Branch Length Verification
Before interpreting the $Z$ axis as "Time before present (Ma)", verify that the tree branch lengths strictly represent evolutionary time.

## 5. No Built-in ASR
Phylo3D-Trait does not reconstruct ancestral states. All tip and ancestral node trait values must be supplied explicitly.

## 6. Standard 6-Step Execution Flow
$$\text{inspect input} \to \text{template-values} \to \text{map traits} \to \text{validate} \to \text{plot} \to \text{verify HTML}$$

## 7. Code Modification Invariant
If code modification is ever required, strictly follow:
$$\text{Understand first} \to \text{Search before coding} \to \text{Reuse} > \text{Extend} > \text{Refactor} > \text{Create}$$
