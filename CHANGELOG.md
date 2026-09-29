# Changelog

All notable user-visible changes are documented here.

## [0.3.2] - 2026-09-30

### Archival and distribution
- Metadata/distribution release intended to trigger Zenodo archival and DOI registration.
- Version metadata synchronized across Python package, citation metadata, README, AI-readable metadata, and CI package smoke tests.
- No scientific computation, CLI semantics, or rendering behavior changed from v0.3.1.

## [0.3.1] - 2026-09-30

### Packaging and distribution
- Added PyPI Trusted Publishing through GitHub Actions using OIDC; no long-lived PyPI token is stored in the repository.
- Added package-build validation for source distribution and wheel artifacts.
- Added install/import/CLI smoke tests for built distributions.
- Improved package, citation, documentation, and contribution metadata for public distribution.

## [0.3.0] - 2026-09-29

### Added
- Four-Layer native WebGL2 renderer with fixed four-layer fragment depth peeling.
- Order-independent transparency for retained peel layers.
- Adaptive front-facing scientific axes.
- Camera-aware outward terminal labels.
- Interactive picking for terminal and internal nodes.
- PNG and hybrid SVG export.
- Public Four-Layer HTML API exports.

### Retained
- Plotly renderer for legacy/general exploration workflows.
- Deterministic clade IDs generated from descendant tip sets.
- CLI workflow for generating node-value templates and rendering supplied continuous traits.

### Scientific scope
Phylo3D-Trait visualizes supplied node values and does not perform ancestral-state reconstruction.
