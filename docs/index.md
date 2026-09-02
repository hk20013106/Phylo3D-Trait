# Phylo3D-Trait Documentation

## Scientific Overview

Phylo3D-Trait is a Python library and CLI tool designed for phylogenetic comparative biology. It bridges the gap between traditional 2D phylogenetic trees and continuous trait evolution by mapping evolutionary trajectories directly into 3D Cartesian coordinates:

- **X (Horizontal layout)**: Positions taxa horizontally based on leaf order and ancestral centroids.
- **Y (Trait / Height)**: Maps the phenotypic trait value directly to height.
- **Z (Evolutionary Time)**: Maps divergence time from root ($Z = \text{root\_age}$) to present day ($Z = 0$).

## Coordinate Specifications

$$\text{node}_i = (X_i, Y_i, Z_i) = (X_i, \text{Trait}_i, \text{Time}_i)$$

1. **Topology & Time Preservation**:
   Branch lengths and divergence times are strictly preserved on the Z axis.
2. **Branch Trait Interpolation**:
   Branches are divided into $N$ linear segments. Along each branch, the trait and height transition smoothly from parent to child:
   $$Y(t) = (1-t) Y_p + t Y_c, \quad t \in [0, 1]$$
3. **Global Colormap**:
   The entire tree is normalized using global $\text{trait}_{\min}$ and $\text{trait}_{\max}$, ensuring that color and height convey identical, calibrated quantitative information across all clades.
4. **Stable Internal Node IDs**:
   Clades are identified by `clade:<hash>` computed from the sorted set of descendant leaves.
