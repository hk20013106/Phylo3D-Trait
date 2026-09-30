# Phylo3D-Trait User & AI Agent Guide

Version target: **0.3.x**

Phylo3D-Trait is a Python CLI/library for **interactive 3D visualization of continuous trait evolution on phylogenetic trees**. It maps tree layout, supplied node-associated trait values, and—when branch lengths represent time—evolutionary age into a rectangular 3D phylogram.

It is a visualization engine. It does **not** infer trees, date trees, fit evolutionary models, or perform ancestral-state reconstruction (ASR).

## 1. Scientific coordinate system

```text
X = tree layout
Y = continuous trait value
Z = branch-length-derived evolutionary depth / time before present
```

For time-calibrated trees with contemporaneous tips, Z can be interpreted as time before present. If branch lengths are substitutions/site or their meaning is unknown, do not label or interpret Z as Ma.

Top branch geometry follows the supplied trait trajectory. Curtain color can either follow geometric height (`height`) or project the local branch trait vertically to the baseline (`branch`).

## 2. Inputs

### Tree
Supported formats:
- Newick
- Nexus

Before scientific interpretation verify:
- tree is parseable;
- tip names are unique;
- rooting is appropriate;
- branch lengths exist;
- branch-length units are known;
- a time interpretation is justified before calling Z “time before present”.

### Trait table
The minimum columns are:

```csv
node_id,trait
Species_A,1.25
Species_B,2.10
clade:xxxxxxxxxxxx,1.64
```

All terminal tips, internal nodes, and the root require explicit numeric values. Missing required node values must fail loudly; never silently impute them.

Internal-node identifiers are deterministic hashes of descendant tip sets. For each internal node, Phylo3D-Trait trims descendant tip names, removes duplicates, sorts them lexicographically, joins them with a literal comma and no spaces, computes SHA-256 on the UTF-8 string, keeps the first 12 lowercase hexadecimal characters, and prefixes `clade:`.

Exact rule:

```python
clade_id = "clade:" + sha256(
    ",".join(sorted(set(descendant_tips))).encode("utf-8")
).hexdigest()[:12]
```

Example:

```text
descendant tips: Species_B, Species_A
canonical key:   Species_A,Species_B
stable node ID:  clade:9a97b9510492
```

This makes the ID independent of child ordering in the Newick tree. Tip-name spelling, capitalization, punctuation, or whitespace after trimming still matter because they change the canonical key.

Users should normally generate these IDs with `template-values` rather than calculate or guess them manually; the generated template includes each `node_id` together with its descendant-tip list.

## 3. Recommended workflow

### Step 1 — generate node template

```bash
python -m phylo3d_trait.cli template-values \
  --tree tree.nwk \
  --output node_values_template.csv
```

The template contains tip names plus stable `clade:<hash>` identifiers for internal nodes/root.

### Step 2 — supply measured/reconstructed values
Fill the `trait` column using the values produced by the user's scientific workflow. ASR may come from tools such as `phytools::fastAnc()`, `ape::ace()`, or other BM/OU/ML/Bayesian workflows; Phylo3D-Trait itself does not calculate these values.

### Step 3 — render
Recommended Four-Layer renderer:

```bash
python -m phylo3d_trait.cli plot \
  --tree tree.nwk \
  --values node_values.csv \
  --output tree3d.html \
  --renderer four-layer \
  --opacity 0.85 \
  --curtain-color-mode branch \
  --centerline-color trait
```

The HTML is standalone and can be opened locally in a modern WebGL2-capable browser.

## 4. Rendering backends

### Four-Layer WebGL2 (`--renderer four-layer`) — recommended
- native WebGL2;
- fixed four-layer fragment depth peeling;
- order-independent transparency for the retained layers;
- supported opacity: **0.5–1.0**;
- adaptive front-facing scientific axes;
- camera-aware outward species labels;
- tip and internal-node picking;
- PNG and hybrid SVG export;
- self-contained HTML with no external JavaScript runtime.

Four-layer depth peeling is bounded rather than infinite-depth compositing. Very deep overlap can omit fragments beyond the peeled layers; the renderer reports the relevant omission bound.

### Plotly (`--renderer plotly`)
Use for legacy workflows or when Plotly-specific interaction/layout behavior is desired. Transparent Mesh3d surfaces use Plotly/WebGL transparency semantics and may show order-dependent artifacts.

## 5. Important display controls

Run the live CLI help before relying on a copied command:

```bash
python -m phylo3d_trait.cli --help
python -m phylo3d_trait.cli template-values --help
python -m phylo3d_trait.cli plot --help
```

Frequently used options:
- `--renderer {four-layer,plotly}`
- `--opacity`
- `--curtain-color-mode {height,branch}`
- `--colorscale`
- `--reverse-colorscale`
- `--centerline-color`
- `--baseline-y`
- `--baseline-raw-value`
- `--trait-display-offset`
- `--trait-display-range`
- `--trait-axis-scale`
- `--camera-preset {elife,root_front,tips_front}`
- `--tip-label-offset`
- `--no-labels`
- `--no-x-axis`, `--no-y-axis`, `--no-z-axis`
- `--no-tip-hover`, `--no-internal-hover`
- `--no-mesh`, `--no-centerline`
- `--show-node-markers`

`--reverse-colorscale` reverses color lookup only; it must not be used as a substitute for changing trait geometry.

`--trait-axis-scale` is visual scaling only; it does not alter scientific trait values.

## 6. Installation and development

Stable installation from PyPI:

```bash
pip install phylo3d-trait
```

For development from source:

```bash
git clone https://github.com/hk20013106/Phylo3D-Trait.git
cd Phylo3D-Trait
pip install -e .
```

For tests:

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

Python support declared by the package: **3.10+**.

## 7. Python API

The CLI is the recommended user interface. Current public imports include:

```python
from phylo3d_trait import (
    parse_tree,
    build_plot_data,
    build_figure,
    build_four_layer_html,
)
from phylo3d_trait.io import load_trait_values
```

For Four-Layer output:

```python
tree = parse_tree("tree.nwk")
traits = load_trait_values("node_values.csv")
plot_data = build_plot_data(tree, traits)

html = build_four_layer_html(
    plot_data,
    opacity=0.85,
    camera_preset="elife",
    centerline_color="trait",
)
with open("tree3d.html", "w", encoding="utf-8") as fh:
    fh.write(html)
```

Always check the current package API/CLI before generating code; do not infer function names from old documentation.

## 8. Agent guardrails

For AI agents and coding assistants:

1. Read `README.md` and this guide.
2. Inspect the actual tree and trait table.
3. Verify branch-length meaning before scientific interpretation.
4. Run `template-values`; do not invent internal-node IDs.
5. Confirm every required node has a numeric trait value.
6. Prefer the existing CLI and renderer architecture.
7. For a new dataset, do not modify package source code.
8. If code changes are genuinely necessary: understand -> search -> reuse -> extend -> refactor -> create.
9. Run targeted tests plus `pytest tests/ -v` after source changes.
10. Report exact command, input paths, output path, validation status, and warnings.

## 9. Non-goals

Phylo3D-Trait is not:
- phylogenetic inference software;
- ASR software;
- sequence-analysis software;
- a tree-dating program;
- an evolutionary model-fitting package.

It is intended to visualize **already defined phylogenies and already supplied continuous node values**.

## 10. Bug reports and contributions

Use GitHub Issues for reproducible bugs and feature requests. For rendering bugs, include:
- operating system;
- browser and version;
- Python version;
- renderer;
- exact CLI command;
- minimal Newick/Nexus tree;
- minimal node-value table;
- screenshot or exported HTML when useful.

Pull requests are welcome. See `CONTRIBUTING.md`.
