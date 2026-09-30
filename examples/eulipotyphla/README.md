# Eulipotyphla real-data example

This directory contains the **starting data** for the 38-species Eulipotyphla hemoglobin-buffering analysis used to validate Phylo3D-Trait.

## Files committed here

- `tree_eulipotyphla.nwk` — dated, ultrametric 38-tip tree; root age = 70.02701 in the supplied branch-length units.
- `tip_traits.csv` — extant/tip hemoglobin buffering trait values.
- `README.md` — provenance and the required analysis path.

Taxon labels are preserved exactly because stable ancestral-node IDs depend on the literal descendant tip names.

### Public-example provenance

The original project files live outside this public repository. The two public example inputs here were recovered from the validated Eulipotyphla rendering artifact used during Phylo3D-Trait development, then checked for internal consistency:

```text
tips = 38
internal/root nodes = 37
all root-to-tip distances = 70.02701
root stable ID = clade:6747b5f19c9e
extant trait range = 5.52690 .. 9.82470
```

They reproduce the validated public example at the stored numerical precision, but they are not presented as a substitute for the upstream raw-analysis archive. If the original source files are later published, compare them directly before replacing these public-example copies.

## Data flow

These two files are the true starting point:

```text
tree_eulipotyphla.nwk
+
tip_traits.csv
        |
        |  phylo3d-trait template-values
        v
node_values_template.csv
        |
        |  external ancestral reconstruction / biochemical calculation
        v
ancestral traits keyed by clade:<hash>
        |
        |  merge tip + ancestral values
        v
node_values_sequence_based.csv
        |
        |  phylo3d-trait plot
        v
interactive HTML
```

Generate all required node IDs from the exact tree:

```bash
phylo3d-trait template-values \
  --tree examples/eulipotyphla/tree_eulipotyphla.nwk \
  --output examples/eulipotyphla/node_values_template.csv
```

Expected checks:

```text
tips = 38
internal/root nodes = 37
total nodes = 75
root stable ID = clade:6747b5f19c9e
```

## Ancestral Hb4 values: sequence-based workflow

The historical Eulipotyphla visualization used character-based Brownian/fastAnc values for internal nodes. Those historical internal values are **not** bundled here as the scientific example.

For the current analysis, ancestral Hb4 values must be derived from ancestral globin sequences:

```text
IQ-TREE 2 ML ancestral HBA_T1 sequence
+
IQ-TREE 2 ML ancestral HBB_T1 sequence
        |
        | same sequence-to-buffer calculation used for extant species
        v
beta_HBA_T1 + beta_HBB_T1
        |
        v
beta_Hb4 = 2 * (beta_HBA_T1 + beta_HBB_T1)
        |
        | map by the same stable clade:<hash> IDs
        v
37 ancestral/root Hb4 values
```

The extant values in `tip_traits.csv` remain unchanged.

The exact 37-node sequence-derived final table is not fabricated here. It should be added as `node_values_sequence_based.csv` only after the upstream sequence-derived values are available and verified. Do **not** substitute the historical Brownian/fastAnc internal-node table.

Once the verified final table exists, the intended rendering command is:

```bash
phylo3d-trait plot \
  --tree examples/eulipotyphla/tree_eulipotyphla.nwk \
  --values examples/eulipotyphla/node_values_sequence_based.csv \
  --output eulipotyphla_Hb_3D_sequence_based.html \
  --renderer four-layer \
  --camera-preset elife \
  --trait-display-offset 4 \
  --baseline-y 0 \
  --trait-axis-scale 0.5 \
  --tip-label-offset 0.03 \
  --reverse-colorscale \
  --curtain-color-mode branch \
  --opacity 0.8
```

See the main [Tutorial](../../docs/tutorial.md) for the generic end-to-end workflow.
