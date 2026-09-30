# Phylo3D-Trait Tutorial

**Goal:** go from a phylogenetic tree and continuous trait values to an interactive 3D HTML visualization.

This tutorial is written for both human users and AI agents. The canonical workflow is:

```text
install
  -> prepare tree
  -> generate node IDs with template-values
  -> fill trait values for every tip + internal node + root
  -> run plot
  -> open the generated HTML
```

> **Important:** Phylo3D-Trait is a visualization tool. It does **not** infer phylogenies, date trees, or reconstruct ancestral trait values. The `template-values` command generates stable ancestral-node IDs only.

---

## 1. Install

Install the stable release from PyPI:

```bash
pip install phylo3d-trait
```

Verify the installation:

```bash
phylo3d-trait --help
```

Expected subcommands include:

```text
plot
template-values
```

If the console command is unavailable in the current environment, the equivalent module form is:

```bash
python -m phylo3d_trait.cli --help
```

Python requirement: **3.10+**.

---

## 2. Prepare the input files

Phylo3D-Trait needs two inputs:

1. a phylogenetic tree in **Newick** or **Nexus** format;
2. a CSV/TSV table containing a numeric trait value for **every tip, every internal node, and the root**.

### 2.1 Tree file

Example `tree.nwk`:

```text
((A:10.0,B:10.0):20.0,(C:15.0,D:15.0):15.0);
```

For a dated/ultrametric tree, branch lengths can be interpreted as evolutionary time and the Z axis can represent time before present.

If branch lengths are substitutions/site or their meaning is unknown, do **not** interpret the Z axis as Ma.

### 2.2 Trait-value file

Do **not** invent ancestral node IDs manually. First generate them from the tree as described in Step 3.

After the IDs are known, the trait table can look like this:

```csv
node_id,trait
A,1.5
B,3.0
C,4.5
D,5.0
clade:b17c8419f544,2.0
clade:6a5756530335,4.0
clade:17f5f129f4c7,1.0
```

The four tip rows are observed/measured tip values in this example. The three `clade:...` rows are values supplied for the two internal ancestors and the root.

These ancestral **trait values are not calculated by Phylo3D-Trait**. They may come from, for example, `phytools::fastAnc()`, `ape::ace()`, Brownian-motion/OU models, Bayesian reconstruction, or another external workflow.

---

## 3. Obtain ancestral-node IDs

Run:

```bash
phylo3d-trait template-values \
  --tree tree.nwk \
  --output node_values_template.csv
```

Equivalent module form:

```bash
python -m phylo3d_trait.cli template-values \
  --tree tree.nwk \
  --output node_values_template.csv
```

For the example tree above, the generated template contains rows equivalent to:

```csv
node_id,label,node_type,descendant_count,descendant_tips,trait
A,A,tip,1,A,
B,B,tip,1,B,
C,C,tip,1,C,
D,D,tip,1,D,
clade:17f5f129f4c7,clade:17f5f129f4c7,root,4,A;B;C;D,
clade:b17c8419f544,clade:b17c8419f544,internal,2,A;B,
clade:6a5756530335,clade:6a5756530335,internal,2,C;D,
```

The `descendant_tips` column tells a human or AI agent exactly which clade each hash refers to:

```text
A;B     -> clade:b17c8419f544
C;D     -> clade:6a5756530335
A;B;C;D -> clade:17f5f129f4c7   (root)
```

### How the ID is generated

For each internal node, Phylo3D-Trait:

1. collects all descendant tip names;
2. strips leading/trailing whitespace;
3. removes duplicates;
4. sorts names lexicographically;
5. joins them with a literal comma and no spaces;
6. computes SHA-256 on the UTF-8 string;
7. keeps the first 12 lowercase hexadecimal characters;
8. prefixes `clade:`.

Exact rule:

```python
clade_id = "clade:" + sha256(
    ",".join(sorted(set(descendant_tips))).encode("utf-8")
).hexdigest()[:12]
```

Example:

```text
descendant tips: B, A
canonical key:   A,B
SHA-256 prefix:  b17c8419f544
node ID:         clade:b17c8419f544
```

This makes the ID independent of child order in the tree: `(A,B)` and `(B,A)` produce the same clade ID.

Tip-name spelling still matters. Changing capitalization, punctuation, or the taxon name changes the hash.

**Recommended practice:** always use `template-values`; do not calculate or guess hashes manually.

### Fill the template

After generating the template, fill the `trait` column for every row:

```csv
node_id,label,node_type,descendant_count,descendant_tips,trait
A,A,tip,1,A,1.5
B,B,tip,1,B,3.0
C,C,tip,1,C,4.5
D,D,tip,1,D,5.0
clade:17f5f129f4c7,Root,root,4,A;B;C;D,1.0
clade:b17c8419f544,Ancestor_AB,internal,2,A;B,2.0
clade:6a5756530335,Ancestor_CD,internal,2,C;D,4.0
```

Save it as, for example, `node_values.csv`.

---

## 4. Run Phylo3D-Trait

Recommended command:

```bash
phylo3d-trait plot \
  --tree tree.nwk \
  --values node_values.csv \
  --output tree3d.html \
  --renderer four-layer \
  --opacity 0.85 \
  --curtain-color-mode branch \
  --centerline-color trait
```

Equivalent module form:

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

A successful run ends with a message similar to:

```text
Successfully generated 3D phylogenetic visualization: tree3d.html
```

Open `tree3d.html` in a modern browser. The Four-Layer renderer produces a standalone HTML file; no web server is required.

### Minimal command

Only three paths are mandatory:

```bash
phylo3d-trait plot \
  --tree tree.nwk \
  --values node_values.csv \
  --output tree3d.html
```

However, the default renderer is currently `plotly`. For publication-oriented transparent curtain rendering, explicitly use:

```text
--renderer four-layer
```

---

## AI execution contract

An AI agent using Phylo3D-Trait should follow these rules exactly:

```text
INPUTS:
- one Newick/Nexus tree
- one numeric continuous trait per tip
- one supplied/reconstructed numeric trait per internal node and root

PROCEDURE:
1. Install phylo3d-trait.
2. Inspect the tree and determine what branch lengths mean.
3. Run template-values on the exact tree.
4. Read node_id + descendant_tips from the generated template.
5. Map externally obtained ancestral trait estimates to those generated node IDs.
6. Confirm every template row has a numeric trait value.
7. Run plot.
8. Report the exact command and output HTML path.

MUST NOT:
- invent clade hashes;
- infer ancestral trait values unless explicitly asked to do so with an external method;
- treat substitution branch lengths as Ma;
- modify package source code merely because the user supplied a new dataset.
```

---

## Built-in example

The repository contains the same four-taxon example used above:

- [Example tree](https://github.com/hk20013106/Phylo3D-Trait/blob/main/examples/example1/tree.nwk)
- [Example node values](https://github.com/hk20013106/Phylo3D-Trait/blob/main/examples/example1/node_values.csv)

For all CLI options, see the [User & AI Agent Guide](PHYLO3D_TRAIT_USAGE_GUIDE.md).
