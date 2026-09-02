"""Phylo3D-Trait: Interactive 3D phylogenetic tree tool with continuous trait evolution.

Visualizes phylogenetic tree topology and evolutionary time combined with
continuous trait evolution in 3D WebGL (Plotly).

Coordinate system:
  X = Tree layout / terminal taxon ordering
  Y = Trait value (Height)
  Z = Evolutionary time / divergence age before present
"""

from phylo3d_trait.models import AnnotatedNode, EdgeSegment, PlotData
from phylo3d_trait.tree import annotate_tree, compute_stable_node_id, parse_tree
from phylo3d_trait.renderer import build_figure, build_plot_data
from phylo3d_trait.template import generate_template_csv

__version__ = "0.2.0"

__all__ = [
    "AnnotatedNode",
    "EdgeSegment",
    "PlotData",
    "annotate_tree",
    "compute_stable_node_id",
    "parse_tree",
    "build_figure",
    "build_plot_data",
    "generate_template_csv",
]
