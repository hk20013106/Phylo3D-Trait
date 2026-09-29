# Changelog

All notable user-visible changes are documented here.

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
