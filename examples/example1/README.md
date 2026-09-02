# Example 1: 4-Taxon Ultrametric Dated Phylogeny

A simple 4-taxon dated tree demonstrating standard 3D rectangular phylogram rendering with continuous trait height and color mapping.

## Tree
- **Number of tips**: 4 (`A`, `B`, `C`, `D`)
- **Topology**: `((A:10,B:10):20,(C:15,D:15):15);`
- **Branch-length meaning**: Evolutionary divergence time (Ma before present)
- **Rooted**: Yes (ultrametric dated tree, root age = 30.0, all tips at $Z = 0.0$)

## Trait
- **Supplied node Trait range**: 1.0 – 5.0 (Root = 1.0, Tip D = 5.0)
- **Visual baseline (`baseline_y`)**: Default (`trait_min` = 1.0)
- **Visualization color range**: 1.0 – 5.0

## Run Command
```bash
python -m phylo3d_trait.cli plot \
  --tree examples/example1/tree.nwk \
  --values examples/example1/node_values.csv \
  --output examples/example1/tree3d.html
```

## Output
- `examples/example1/tree3d.html`
