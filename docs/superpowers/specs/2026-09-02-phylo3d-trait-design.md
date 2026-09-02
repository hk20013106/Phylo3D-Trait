# Phylo3D-Trait v0.1 Design

## Goal

Build a minimal, reusable renderer that converts a phylogenetic tree plus user-supplied trait values for every tip and internal node into an interactive 3D phylogeny. The result must open locally in a browser, support arbitrary mouse rotation/zoom, and preserve the chosen axis scaling during camera rotation.

## Scope

v0.1 is intentionally small. It does not estimate ancestral states, align sequences, infer trees, place mutation labels, add species images, or reproduce every publication-layout detail from Figure 5. It only solves the core reusable problem:

`tree + node/tip trait values -> 3D coordinates -> interactive HTML`

## Inputs

### 1. Tree

Accepted format for v0.1: Newick (`.nwk`, `.newick`, `.tree`).

Requirements:
- rooted tree;
- branch lengths required;
- tip labels unique;
- internal node labels strongly recommended when supplied by the user's upstream ancestral-state workflow.

### 2. Trait table

CSV with exactly these required columns:

```text
node,value
```

`node` may identify either a tip or an internal node. Every tree node must resolve to exactly one trait value.

Example:

```csv
node,value
A,1.2
B,2.4
C,0.8
D,1.5
ANC_AB,1.8
ANC_CD,1.1
ROOT,1.0
```

## Internal-node identity

Scientific correctness requires node values not to drift onto the wrong ancestor after tree reordering. v0.1 therefore uses two resolution modes:

1. **Named-node mode (preferred):** exact match to the internal node label already present in the Newick tree.
2. **Clade-signature mode (fallback):** for unlabeled internal nodes, derive a stable ID from the sorted descendant tip names. Tree rotation/ladderization does not change this signature.

The program must fail closed if an input row is unmatched, duplicated, ambiguous, or if any tree node lacks a value. It must never silently guess a node mapping.

## Coordinate model

For every node `i`, construct one point:

```text
P_i = (x_i, y_i, z_i)
```

where:
- `x_i` = topology layout coordinate used to separate tips/clades;
- `y_i` = evolutionary time coordinate derived from cumulative branch length;
- `z_i` = supplied trait value.

### X: topology layout

Tips receive evenly spaced x positions in traversal order. Each internal node receives the mean x position of its direct children. This produces a conventional rectangular-phylogram layout embedded in 3D.

### Y: time

Root-to-node path length is computed from branch lengths. For an ultrametric tree, present-day tips share the same time coordinate. Display mode defaults to `time-before-present`, so present = 0 and deeper ancestors have larger positive ages.

For non-ultrametric trees, v0.1 must not invent ages. It may still render branch-length distance, but must label the axis `Root-to-node branch length` unless the user explicitly requests time-before-present and the tree is ultrametric within tolerance.

### Z: trait

Trait value is taken directly from the input table. No ancestral-state estimation occurs in the renderer.

## Edge geometry

For each parent-child edge, v0.1 renders a continuous 3D segment between the parent and child coordinates.

Default interpolation is linear:

```text
z(s) = z_parent + s * (z_child - z_parent), 0 <= s <= 1
```

This interpolation is a visualization convention, not a claim that the biological trait changed at a constant evolutionary rate. The README must state this explicitly.

## Rendering

Technology:
- Python >= 3.10
- Biopython for Newick parsing
- pandas for trait-table parsing
- Plotly for WebGL 3D rendering and standalone HTML output

v0.1 rendering layers:
1. 3D branch lines with color mapped continuously to trait value;
2. optional node markers;
3. tip labels;
4. hover text with node ID, branch/time coordinate, and trait value;
5. fixed user-configurable `scene.aspectmode="manual"` so camera rotation does not distort configured x:y:z scale;
6. optional orthographic camera mode to remove perspective size distortion.

The initial version uses lines, not Figure-5-style filled ribbons/walls. Ribbon/wall geometry is a v0.2 enhancement after the coordinate engine is validated.

## Command-line interface

Primary command:

```bash
python -m phylo3d_trait tree.nwk node_values.csv -o tree3d.html
```

Required behavior:
- parse tree and values;
- validate complete node mapping;
- compute coordinates;
- render standalone HTML;
- print a concise summary of number of tips, internal nodes, edges, and output path.

Useful v0.1 options:

```text
--projection perspective|orthographic
--x-scale FLOAT
--y-scale FLOAT
--z-scale FLOAT
--show-nodes / --hide-nodes
--show-tip-labels / --hide-tip-labels
--title TEXT
```

## Output

Primary output: one self-contained or CDN-backed Plotly HTML file that can be opened locally in a browser and manipulated with mouse drag/zoom.

Default should favor a compact HTML by using Plotly JS from CDN. A `--self-contained` option may embed Plotly JS for offline use.

## Example dataset

Repository must include a small synthetic tree with named internal nodes and a complete trait table. The example should deliberately contain one clade with increasing trait values so the 3D displacement is visually obvious.

Expected quick-start:

```bash
pip install -e .
phylo3d-trait examples/example_tree.nwk examples/example_traits.csv -o example.html
```

## Error handling

Fail with a non-zero exit status for:
- missing input file;
- malformed Newick;
- missing branch length;
- duplicate tip label;
- duplicate trait-table node identifier;
- non-numeric trait value;
- missing node value;
- unmatched node value;
- ambiguous internal-node mapping.

Do not silently drop or impute data.

## Tests

Automated tests must cover at minimum:
- tree parsing and node counting;
- deterministic tip/internal x layout;
- branch-length/time coordinate calculation;
- trait mapping for tips and named ancestors;
- stable clade-signature generation for unlabeled ancestors;
- failure on incomplete/ambiguous mapping;
- edge interpolation endpoints;
- HTML creation containing a Plotly 3D scene.

## Repository layout

```text
Phylo3D-Trait/
├── pyproject.toml
├── README.md
├── src/
│   └── phylo3d_trait/
│       ├── __init__.py
│       ├── cli.py
│       ├── tree.py
│       ├── traits.py
│       ├── layout.py
│       └── render.py
├── examples/
│   ├── example_tree.nwk
│   └── example_traits.csv
├── tests/
│   ├── test_tree.py
│   ├── test_traits.py
│   ├── test_layout.py
│   └── test_render.py
└── docs/
```

Responsibilities:
- `tree.py`: parse/validate tree and generate stable node identities;
- `traits.py`: parse/validate trait table and resolve values to tree nodes;
- `layout.py`: compute x/y/z node coordinates and edge samples;
- `render.py`: convert coordinates to Plotly traces and write HTML;
- `cli.py`: command-line boundary only.

## Acceptance criteria for v0.1

v0.1 is accepted when all of the following are true:

1. A supplied Newick tree and complete node-value CSV generate an HTML file with no manual editing.
2. Browser mouse interaction supports arbitrary rotation and zoom.
3. Axis scaling remains fixed according to configured x:y:z proportions while the camera rotates.
4. Every rendered endpoint corresponds exactly to the supplied node value.
5. Internal-node mapping is deterministic and fails rather than guessing when unsafe.
6. Automated tests pass.
7. The bundled example can be generated with one documented command.

## Deferred to v0.2+

Only after the coordinate/mapping engine is validated:
- Figure-5-style vertical ribbons/walls or tubes;
- continuous surface color gradients on filled branch geometry;
- geologic-period background planes;
- mutation/event labels on branches;
- species-image placement;
- static publication-quality export presets;
- Nexus input;
- richer event-based rather than linear within-branch trait trajectories.
