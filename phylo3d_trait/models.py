"""Data models for Phylo3D-Trait visualization.

Defines the core data structures for representing annotated nodes, interpolated
edge segments, and complete 3D plot datasets.

Coordinates:
  X: Phylogeny horizontal layout (tip ordering and internal node centroids)
  Y: Trait value (Height in 3D space)
  Z: Evolutionary time / divergence age (tips at 0, root at maximum age)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class AnnotatedNode:
    """Represents a phylogenetic node (tip or internal) with 3D coordinates and trait values.

    Attributes:
        node_id: Unique identifier (taxon name for tips, 'clade:<hash>' for internal nodes).
        label: Display name / label.
        is_tip: True if this is a leaf/tip taxon, False if internal ancestral node.
        x: Horizontal tree layout coordinate.
        y: Display Trait value, mapped to height (Y axis).
        z: Evolutionary time coordinate (Z axis, tips at 0, root at max age).
        trait: Trait value mapped to Y (equal to display_trait).
        time: Evolutionary time / age before present (identical to z).
        raw_trait: Original unscaled scientific trait value.
        display_trait: Linearly rescaled trait value used for 3D Y coordinate and color.
        parent_id: Node ID of immediate parent, or None for root.
        children_ids: List of immediate children node IDs.
        branch_length: Evolutionary branch length from parent.
        descendant_tips: Sorted list of all leaf taxon names under this node.
    """

    node_id: str
    label: str
    is_tip: bool
    x: float
    y: float
    z: float
    trait: float
    time: float
    raw_trait: Optional[float] = None
    display_trait: Optional[float] = None
    parent_id: Optional[str] = None
    children_ids: List[str] = field(default_factory=list)
    branch_length: Optional[float] = None
    descendant_tips: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.raw_trait is None:
            self.raw_trait = self.trait
        if self.display_trait is None:
            self.display_trait = self.y
        # Enforce consistency: y is display_trait, trait is display_trait, z is time
        self.y = self.display_trait
        self.trait = self.display_trait
        if self.z != self.time:
            self.z = self.time


@dataclass
class EdgeSegment:
    """Represents a subdivided linear segment of a phylogenetic branch in 3D space.

    Attributes:
        parent_id: Node ID of the parent node.
        child_id: Node ID of the child node.
        x0: Start X coordinate.
        y0: Start Y coordinate (interpolated display trait value).
        z0: Start Z coordinate (interpolated evolutionary time).
        x1: End X coordinate.
        y1: End Y coordinate (interpolated display trait value).
        z1: End Z coordinate (interpolated evolutionary time).
        trait0: Start display trait value (equal to y0).
        trait1: End display trait value (equal to y1).
        raw_trait0: Start raw scientific trait value.
        raw_trait1: End raw scientific trait value.
        display_trait0: Start display trait value.
        display_trait1: End display trait value.
        segment_index: Index of this segment along the edge (0 .. total_segments - 1).
        total_segments: Total number of subdivided segments along this branch.
        segment_type: 'connector' (horizontal layout adjust) or 'lineage' (through-time evolution).
    """

    parent_id: str
    child_id: str
    x0: float
    y0: float
    z0: float
    x1: float
    y1: float
    z1: float
    trait0: float
    trait1: float
    raw_trait0: Optional[float] = None
    raw_trait1: Optional[float] = None
    display_trait0: Optional[float] = None
    display_trait1: Optional[float] = None
    segment_index: int = 0
    total_segments: int = 1
    segment_type: str = "lineage"  # 'connector' or 'lineage'

    def __post_init__(self) -> None:
        if self.raw_trait0 is None:
            self.raw_trait0 = self.trait0
        if self.raw_trait1 is None:
            self.raw_trait1 = self.trait1
        if self.display_trait0 is None:
            self.display_trait0 = self.y0
        if self.display_trait1 is None:
            self.display_trait1 = self.y1
        self.y0 = self.display_trait0
        self.y1 = self.display_trait1
        self.trait0 = self.display_trait0
        self.trait1 = self.display_trait1


@dataclass
class PlotData:
    """Encapsulates all 3D geometry, annotations, and global scaling for plotting.

    Attributes:
        nodes: Dictionary mapping node_id to AnnotatedNode.
        segments: List of interpolated EdgeSegments connecting nodes.
        trait_min: Global minimum display trait value across all nodes.
        trait_max: Global maximum display trait value across all nodes.
        time_min: Minimum evolutionary time (0.0 at tips).
        time_max: Maximum evolutionary time (root divergence age).
        x_min: Minimum X layout coordinate.
        x_max: Maximum X layout coordinate.
        raw_trait_min: Global minimum raw scientific trait value across all nodes.
        raw_trait_max: Global maximum raw scientific trait value across all nodes.
        trait_display_range: Optional target display range (start, end) used for linear remapping.
        colorscale: Plotly continuous colorscale name (e.g., 'Turbo').
        title: Title of the visualization.
        baseline_y: Custom or default baseline Y plane height.
    """

    nodes: Dict[str, AnnotatedNode]
    segments: List[EdgeSegment]
    trait_min: float
    trait_max: float
    time_min: float
    time_max: float
    x_min: float
    x_max: float
    raw_trait_min: float = 0.0
    raw_trait_max: float = 0.0
    trait_display_range: Optional[tuple[float, float]] = None
    colorscale: str = "Turbo"
    title: str = "3D Phylogenetic Tree with Continuous Trait Evolution"
    baseline_y: Optional[float] = None

    def __post_init__(self) -> None:
        if self.baseline_y is None:
            self.baseline_y = self.trait_min

    def raw_to_display(self, raw_val: float) -> float:
        """Convert raw scientific trait value to display coordinate Y."""
        if self.trait_display_range is None:
            return raw_val
        t_start, t_end = self.trait_display_range
        if self.raw_trait_max == self.raw_trait_min:
            return t_start
        return t_start + (raw_val - self.raw_trait_min) * (t_end - t_start) / (self.raw_trait_max - self.raw_trait_min)

    def display_to_raw(self, display_val: float) -> float:
        """Convert display coordinate Y back to scientific raw trait value (inverse transform)."""
        if self.trait_display_range is None:
            return display_val
        t_start, t_end = self.trait_display_range
        if t_end == t_start:
            return self.raw_trait_min
        return self.raw_trait_min + (display_val - t_start) * (self.raw_trait_max - self.raw_trait_min) / (t_end - t_start)


