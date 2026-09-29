# Phylo3D-Trait

**Interactive 3D visualization of continuous trait evolution on phylogenetic trees.**

Phylo3D-Trait is a Python CLI/library for visualizing supplied continuous trait values and ancestral-state estimates across phylogenies. For dated trees, evolutionary depth can be displayed as time before present. The recommended renderer produces standalone WebGL2 HTML using fixed four-layer depth peeling for transparent curtain meshes.

![Phylo3D-Trait example](assets/preview_eulipotyphla.png)

## What it accepts
- Newick or Nexus phylogenetic trees
- CSV/TSV node values for tips, internal nodes, and root
- continuous phenotypic, physiological, morphological, genomic, or other quantitative traits

## What it produces
- interactive 3D rectangular phylograms
- continuous trait height/color mapping
- tip and ancestral-node inspection
- standalone offline-viewable HTML
- publication-oriented PNG and hybrid SVG export

## What it does not do
Phylo3D-Trait does not infer phylogenies, date trees, fit evolutionary models, or perform ancestral-state reconstruction.

## Quick start

```bash
git clone https://github.com/hk20013106/Phylo3D-Trait.git
cd Phylo3D-Trait
pip install -e .

python -m phylo3d_trait.cli template-values -t tree.nwk -o values_template.csv

python -m phylo3d_trait.cli plot \
  -t tree.nwk \
  -v values.csv \
  -o tree3d.html \
  --renderer four-layer \
  --opacity 0.85
```

## Documentation
- [Repository README](https://github.com/hk20013106/Phylo3D-Trait)
- [User & AI Agent Guide](PHYLO3D_TRAIT_USAGE_GUIDE.md)
- [Contributing](https://github.com/hk20013106/Phylo3D-Trait/blob/main/CONTRIBUTING.md)
- [Issues](https://github.com/hk20013106/Phylo3D-Trait/issues)
- [Citation metadata](https://github.com/hk20013106/Phylo3D-Trait/blob/main/CITATION.cff)
