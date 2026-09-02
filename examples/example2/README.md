# Example 2: 6-Taxon Multi-Level Phylogeny with Baseline Y = 0

A 6-taxon nested phylogeny demonstrating rectangular 3D phylogram rendering with a custom visual baseline `baseline_y = 0.0` ensuring positive curtain elevation at all ancestral nodes.

## Tree
- **Number of tips**: 6 (`A`, `B`, `C`, `D`, `E`, `F`)
- **Topology**: `(((A:10,B:10):15,(C:15,D:15):10):25,(E:20,F:20):30);`
- **Branch-length meaning**: Evolutionary divergence time (Ma before present)
- **Rooted**: Yes (ultrametric dated tree, root age = 50.0, all tips at $Z = 0.0$)

## Trait
- **Supplied node Trait range**: 1.0 – 5.0 (Root = 1.0, Tip D = 5.0)
- **Visual baseline (`baseline_y`)**: `0.0` (via `--baseline-y 0`)
- **Visualization color range**: 0.0 – 5.0 (Mesh color range strictly spans $[0.0, 5.0]$)
- **Root curtain height**: Top $Y = 1.0$, Bottom $Y = 0.0 \implies \text{height} = 1.0$ (preserves visible立面, does not collapse into a line)

## Run Command
```bash
python -m phylo3d_trait.cli plot \
  --tree examples/example2/tree.nwk \
  --values examples/example2/node_values.csv \
  --output results/example2/tree3d.html \
  --baseline-y 0
```

## Output
- `results/example2/tree3d.html`
