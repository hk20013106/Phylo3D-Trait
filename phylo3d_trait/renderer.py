"""Plotly 3D renderer for Phylo3D-Trait visualization.

Builds interactive 3D WebGL plots mapping:
- X axis: Tree layout
- Y axis: Trait value (Height)
- Z axis: Evolutionary time before present

Visual representations:
- Continuous vertical curtain / ribbon surfaces (Mesh3d) descending from each
  branch's trait height down to a common trait baseline plane.
- Crisp branch top outlines (Scatter3d lines).
- Text labels for terminal taxa without intrusive marker dots.
- Pure white / transparent background with clean axis gridlines.
- eLife-style camera preset with screen-vertical Y (Trait) and +Z foreground (MRCA).
- Global trait normalization (trait_min, trait_max) across all surfaces.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import plotly.graph_objects as go

from phylo3d_trait.models import EdgeSegment, PlotData
from phylo3d_trait.tree import annotate_tree, parse_tree


CAMERA_PRESETS: Dict[str, Dict[str, Any]] = {
    "elife": dict(
        eye=dict(x=1.35, y=0.65, z=2.1),
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        projection=dict(type="orthographic"),
    ),
    "root_front": dict(
        eye=dict(x=0.0, y=0.5, z=2.4),
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        projection=dict(type="orthographic"),
    ),
    "tips_front": dict(
        eye=dict(x=1.35, y=0.65, z=-2.1),
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        projection=dict(type="orthographic"),
    ),
}


def build_plot_data(
    tree_input: Any,
    trait_values: Dict[str, float],
    num_segments: int = 10,
    colorscale: str = "Turbo",
    title: str = "3D Phylogenetic Tree with Continuous Trait Evolution",
    baseline_y: Optional[float] = None,
    trait_display_range: Optional[Tuple[float, float]] = None,
) -> PlotData:
    """Convenience helper to parse tree, annotate traits, and construct PlotData.

    Args:
        tree_input: Newick/Nexus string, path, or Tree object.
        trait_values: Dictionary of node IDs to trait values.
        num_segments: Number of interpolation segments per branch.
        colorscale: Plotly colorscale name.
        title: Visualization title.
        baseline_y: Optional custom baseline Trait height for curtain meshes.
        trait_display_range: Optional custom (start, end) target display range for linear remapping.

    Returns:
        PlotData object.
    """
    from Bio.Phylo.BaseTree import Clade, Tree

    if isinstance(tree_input, (Tree, Clade)) or hasattr(tree_input, "get_terminals"):
        tree = tree_input
    else:
        tree = parse_tree(tree_input)

    data = annotate_tree(
        tree_input=tree,
        trait_values=trait_values,
        num_segments=num_segments,
        colorscale=colorscale,
        title=title,
        trait_display_range=trait_display_range,
    )
    if baseline_y is not None:
        data.baseline_y = baseline_y
    return data


def _build_branch_curtains_geometry(
    plot_data: PlotData,
    baseline_y: float,
) -> Tuple[List[float], List[float], List[float], List[int], List[int], List[int], List[float]]:
    """Construct 3D mesh vertices, triangle indices, and vertex color intensities for branch curtains.

    For every parent -> child edge:
      - Obtains the sequence of sampled vertices P_0 .. P_M along the branch.
      - Constructs Top_k = (x_k, y_k, z_k) at the branch trait height (intensity = y_k).
      - Constructs Bottom_k = (x_k, baseline_y, z_k) on the baseline plane (intensity = baseline_y).
      - Generates 2 triangles for each adjacent step (k, k+1):
          Triangle A: (Top_k, Bottom_k, Top_{k+1})
          Triangle B: (Bottom_k, Bottom_{k+1}, Top_{k+1})
      - Each edge is triangulated independently, strictly preserving topology
        without cross-branch Delaunay triangulation.

    Args:
        plot_data: PlotData containing edge segments and scaling bounds.
        baseline_y: The constant Y height of the baseline plane.

    Returns:
        Tuple of (mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity).
    """
    mesh_x: List[float] = []
    mesh_y: List[float] = []
    mesh_z: List[float] = []
    mesh_i: List[int] = []
    mesh_j: List[int] = []
    mesh_k: List[int] = []
    mesh_intensity: List[float] = []

    # Group segments by parent-child edge
    edge_map: Dict[tuple, List[EdgeSegment]] = {}
    for seg in plot_data.segments:
        key = (seg.parent_id, seg.child_id)
        if key not in edge_map:
            edge_map[key] = []
        edge_map[key].append(seg)

    vertex_offset = 0

    for (p_id, c_id), segs in edge_map.items():
        sorted_segs = sorted(segs, key=lambda s: s.segment_index)
        if not sorted_segs:
            continue

        # Extract sequence of points along this branch
        branch_pts = [(sorted_segs[0].x0, sorted_segs[0].y0, sorted_segs[0].z0)]
        for s in sorted_segs:
            branch_pts.append((s.x1, s.y1, s.z1))

        num_pts = len(branch_pts)

        # Add top and bottom vertices for this branch
        for k in range(num_pts):
            xk, yk, zk = branch_pts[k]

            # Top vertex (at trait height, intensity == yk)
            mesh_x.append(xk)
            mesh_y.append(yk)
            mesh_z.append(zk)
            mesh_intensity.append(yk)

            # Bottom vertex (at baseline_y, intensity == baseline_y)
            mesh_x.append(xk)
            mesh_y.append(baseline_y)
            mesh_z.append(zk)
            mesh_intensity.append(baseline_y)

        # Build 2 triangles per segment quad
        for k in range(num_pts - 1):
            top_k = vertex_offset + 2 * k
            bot_k = vertex_offset + 2 * k + 1
            top_k1 = vertex_offset + 2 * (k + 1)
            bot_k1 = vertex_offset + 2 * (k + 1) + 1

            # Triangle A: (Top_k, Bottom_k, Top_{k+1})
            mesh_i.append(top_k)
            mesh_j.append(bot_k)
            mesh_k.append(top_k1)

            # Triangle B: (Bottom_k, Bottom_{k+1}, Top_{k+1})
            mesh_i.append(bot_k)
            mesh_j.append(bot_k1)
            mesh_k.append(top_k1)

        vertex_offset += 2 * num_pts

    return mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity


def _generate_rescaled_ticks(
    plot_data: PlotData,
    baseline_y: Optional[float] = None,
    num_ticks: int = 5,
) -> Tuple[List[float], List[str]]:
    """Generate (tickvals, ticktext) in display space labeled with raw scientific trait values.

    Args:
        plot_data: PlotData container with raw trait bounds and trait_display_range.
        baseline_y: Custom baseline Y plane height.
        num_ticks: Number of trait tick steps across the display domain (default: 5).

    Returns:
        Tuple of (tickvals, ticktext) for Plotly axis and colorbar.
    """
    if plot_data.trait_display_range is None:
        return [], []

    d_start, d_end = plot_data.trait_display_range
    d_min = min(d_start, d_end)
    d_max = max(d_start, d_end)

    if d_max == d_min:
        disp_ticks = [d_min]
    else:
        step = (d_max - d_min) / float(num_ticks - 1)
        disp_ticks = [d_min + i * step for i in range(num_ticks)]

    tickvals: List[float] = []
    ticktext: List[str] = []

    eff_baseline = baseline_y if baseline_y is not None else plot_data.baseline_y
    # If baseline is explicitly below the display trait domain, label it as "baseline"
    if eff_baseline is not None and eff_baseline < d_min - 1e-4:
        tickvals.append(float(eff_baseline))
        ticktext.append("baseline")

    for d_val in disp_ticks:
        raw_val = plot_data.display_to_raw(d_val)
        tickvals.append(round(d_val, 6))
        if abs(raw_val - round(raw_val)) < 1e-6:
            label = f"{int(round(raw_val))}"
        else:
            label = f"{raw_val:.4f}".rstrip("0").rstrip(".")
        ticktext.append(label)

    return tickvals, ticktext


def build_figure(
    plot_data: PlotData,
    title: Optional[str] = None,
    branch_width: float = 1.0,
    show_tip_labels: bool = True,
    aspect_ratio: Optional[Dict[str, float]] = None,
    show_mesh: bool = True,
    mesh_opacity: float = 1.0,
    show_centerline: bool = True,
    centerline_color: str = "dark",
    baseline_y: Optional[float] = None,
    show_node_markers: bool = False,
    internal_marker_size: float = 4.0,
    background: str = "white",
    camera_preset: str = "elife",
    custom_camera: Optional[Dict[str, Any]] = None,
) -> go.Figure:
    """Construct an interactive Plotly 3D Figure in clean publication style with eLife-style camera.

    Args:
        plot_data: PlotData containing annotated nodes, edge segments, and scaling limits.
        title: Optional title override.
        branch_width: Line width for 3D branch top outline (default: 1.0).
        show_tip_labels: Whether to display text labels for tip taxa (default: True).
        aspect_ratio: Optional custom aspect ratio dictionary {'x': float, 'y': float, 'z': float}.
        show_mesh: Whether to render continuous vertical curtain meshes.
        mesh_opacity: Opacity for curtain meshes (0.0 to 1.0, default 1.0).
        show_centerline: Whether to render top-edge outline along branches.
        centerline_color: Color mode for centerline ('dark', 'trait', or CSS color).
        baseline_y: Custom baseline Y plane height (defaults to plot_data.baseline_y).
        show_node_markers: Whether to render ancestral node markers (default: False).
        internal_marker_size: Marker size if show_node_markers is True.
        background: 'white' (default) or 'transparent'.
        camera_preset: Preset camera angle ('elife', 'root_front', 'tips_front', default: 'elife').
        custom_camera: Optional dict to override camera config completely.

    Returns:
        Plotly go.Figure configured for interactive 3D display.
    """
    plot_title = title if title is not None else plot_data.title

    eff_baseline_y = baseline_y if baseline_y is not None else plot_data.baseline_y
    if eff_baseline_y is None:
        eff_baseline_y = plot_data.trait_min

    # Color range covers all traits and baseline plane
    cmin = min(plot_data.trait_min, eff_baseline_y)
    cmax = max(plot_data.trait_max, eff_baseline_y)
    # Handle single constant trait edge case
    if cmin == cmax:
        cmin -= 0.5
        cmax += 0.5

    fig = go.Figure()

    # 1. Build continuous curtain mesh surfaces (Mesh3d)
    is_transformed = plot_data.trait_display_range is not None
    colorbar_title = "Trait Value"
    y_axis_title = "Trait value"
    rescaled_tickvals, rescaled_ticktext = _generate_rescaled_ticks(
        plot_data=plot_data,
        baseline_y=eff_baseline_y,
        num_ticks=5,
    )

    if show_mesh and plot_data.segments:
        (
            mesh_x,
            mesh_y,
            mesh_z,
            mesh_i,
            mesh_j,
            mesh_k,
            mesh_intensity,
        ) = _build_branch_curtains_geometry(
            plot_data=plot_data,
            baseline_y=eff_baseline_y,
        )

        if mesh_x and mesh_i:
            cb_dict = dict(
                title=dict(text=colorbar_title, side="top", font=dict(size=12, color="#333333")),
                thickness=18,
                len=0.75,
                x=1.02,
            )
            if is_transformed and rescaled_tickvals:
                cb_dict["tickmode"] = "array"
                cb_dict["tickvals"] = rescaled_tickvals
                cb_dict["ticktext"] = rescaled_ticktext

            fig.add_trace(
                go.Mesh3d(
                    x=mesh_x,
                    y=mesh_y,  # Y is TRAIT (Top) and baseline_y (Bottom)
                    z=mesh_z,  # Z is TIME
                    i=mesh_i,
                    j=mesh_j,
                    k=mesh_k,
                    intensity=mesh_intensity,
                    colorscale=plot_data.colorscale,
                    cmin=cmin,
                    cmax=cmax,
                    opacity=mesh_opacity,
                    flatshading=False,
                    lighting=dict(
                        ambient=0.85,
                        diffuse=0.5,
                        specular=0.08,
                        roughness=0.8,
                    ),
                    hoverinfo="none",
                    name="Branch Curtains",
                    showscale=True,
                    colorbar=cb_dict,
                )
            )

    # 2. Build branch top centerline outlines (Scatter3d lines)
    branch_x: List[Optional[float]] = []
    branch_y: List[Optional[float]] = []
    branch_z: List[Optional[float]] = []
    branch_colors: List[Optional[float]] = []

    edge_map: Dict[tuple, List[EdgeSegment]] = {}
    for seg in plot_data.segments:
        key = (seg.parent_id, seg.child_id)
        if key not in edge_map:
            edge_map[key] = []
        edge_map[key].append(seg)

    for (p_id, c_id), segs in edge_map.items():
        sorted_segs = sorted(segs, key=lambda s: s.segment_index)
        if not sorted_segs:
            continue

        # Start vertex
        branch_x.append(sorted_segs[0].x0)
        branch_y.append(sorted_segs[0].y0)
        branch_z.append(sorted_segs[0].z0)
        branch_colors.append(sorted_segs[0].trait0)

        # End vertices
        for s in sorted_segs:
            branch_x.append(s.x1)
            branch_y.append(s.y1)
            branch_z.append(s.z1)
            branch_colors.append(s.trait1)

        # Disconnect line from next branch
        branch_x.append(None)
        branch_y.append(None)
        branch_z.append(None)
        branch_colors.append(sorted_segs[-1].trait1)

    if show_centerline and branch_x:
        if centerline_color == "trait":
            line_cfg = dict(
                color=branch_colors,
                colorscale=plot_data.colorscale,
                cmin=cmin,
                cmax=cmax,
                width=branch_width,
            )
        elif centerline_color == "dark":
            line_cfg = dict(
                color="#2b2b2b",
                width=branch_width,
            )
        else:
            line_cfg = dict(
                color=centerline_color,
                width=branch_width,
            )

        fig.add_trace(
            go.Scatter3d(
                x=branch_x,
                y=branch_y,  # Y is TRAIT
                z=branch_z,  # Z is TIME
                mode="lines",
                line=line_cfg,
                hoverinfo="none",
                name="Branch Centerlines",
                showlegend=False,
            )
        )

    # 3. Optional Internal Nodes trace (Default: False)
    if show_node_markers:
        internal_nodes = [n for n in plot_data.nodes.values() if not n.is_tip]
        if internal_nodes:
            if is_transformed:
                customdata_internal = [
                    [
                        n.node_id,
                        "Ancestral Node",
                        n.raw_trait,
                        n.display_trait,
                        f"Descendants ({len(n.descendant_tips)} tips): {', '.join(n.descendant_tips[:3])}{'...' if len(n.descendant_tips) > 3 else ''}",
                    ]
                    for n in internal_nodes
                ]
                hovertemplate_internal = (
                    "<b>Node: %{customdata[0]}</b><br>"
                    "Type: %{customdata[1]}<br>"
                    "Raw Trait: %{customdata[2]:.4f}<br>"
                    "Display Trait (internal Y): %{customdata[3]:.4f}<br>"
                    "Time before present (Z): %{z:.4f}<br>"
                    "Tree Layout (X): %{x:.2f}<br>"
                    "%{customdata[4]}<extra></extra>"
                )
            else:
                customdata_internal = [
                    [
                        n.node_id,
                        "Ancestral Node",
                        f"Descendants ({len(n.descendant_tips)} tips): {', '.join(n.descendant_tips[:3])}{'...' if len(n.descendant_tips) > 3 else ''}",
                    ]
                    for n in internal_nodes
                ]
                hovertemplate_internal = (
                    "<b>Node: %{customdata[0]}</b><br>"
                    "Type: %{customdata[1]}<br>"
                    "Trait value (Y / Height): %{y:.4f}<br>"
                    "Time before present (Z): %{z:.4f}<br>"
                    "Tree Layout (X): %{x:.2f}<br>"
                    "%{customdata[2]}<extra></extra>"
                )

            fig.add_trace(
                go.Scatter3d(
                    x=[n.x for n in internal_nodes],
                    y=[n.y for n in internal_nodes],  # Y is TRAIT
                    z=[n.z for n in internal_nodes],  # Z is TIME
                    mode="markers",
                    marker=dict(
                        size=internal_marker_size,
                        color=[n.trait for n in internal_nodes],
                        colorscale=plot_data.colorscale,
                        cmin=cmin,
                        cmax=cmax,
                        symbol="diamond",
                        opacity=0.95,
                    ),
                    customdata=customdata_internal,
                    hovertemplate=hovertemplate_internal,
                    name="Internal Nodes",
                    showlegend=False,
                )
            )

    # 4. Build Terminal Tips text trace (aligned on top-front reference line: y = global_trait_max, z = 0.0)
    tip_nodes = [n for n in plot_data.nodes.values() if n.is_tip]
    if tip_nodes and show_tip_labels:
        label_x = [n.x for n in tip_nodes]
        label_y = [plot_data.trait_max for _ in tip_nodes]
        label_z = [0.0 for _ in tip_nodes]
        if is_transformed:
            customdata_tip = [[n.label, n.node_id, n.raw_trait, n.display_trait] for n in tip_nodes]
            hovertemplate_tip = (
                "<b>Taxon: %{customdata[0]}</b><br>"
                "Node ID: %{customdata[1]}<br>"
                "Raw Trait: %{customdata[2]:.4f}<br>"
                "Display Trait (internal Y): %{customdata[3]:.4f}<br>"
                "Time before present: 0.0000<br>"
                "Tree Layout (X): %{x:.2f}<extra></extra>"
            )
        else:
            customdata_tip = [[n.label, n.node_id, n.raw_trait] for n in tip_nodes]
            hovertemplate_tip = (
                "<b>Taxon: %{customdata[0]}</b><br>"
                "Node ID: %{customdata[1]}<br>"
                "Tip Trait value: %{customdata[2]:.4f}<br>"
                "Time before present: 0.0000<br>"
                "Tree Layout (X): %{x:.2f}<extra></extra>"
            )

        fig.add_trace(
            go.Scatter3d(
                x=label_x,
                y=label_y,  # Aligned uniformly at global maximum trait value
                z=label_z,  # Aligned at Time before present = 0.0
                mode="text",
                text=[n.label for n in tip_nodes],
                textposition="top center",
                textfont=dict(size=11, color="#222222"),
                customdata=customdata_tip,
                hovertemplate=hovertemplate_tip,
                name="Terminal Taxa",
                showlegend=False,
            )
        )

    # Calculate default balanced manual aspect ratio
    if aspect_ratio is None:
        ratio_x = 1.4
        ratio_y = 1.0
        ratio_z = 1.2
        ratio_dict = dict(x=ratio_x, y=ratio_y, z=ratio_z)
    else:
        ratio_dict = aspect_ratio

    # Background colors
    is_transparent = background.lower() == "transparent"
    paper_bg = "rgba(0,0,0,0)" if is_transparent else "white"
    plot_bg = "rgba(0,0,0,0)" if is_transparent else "white"

    # Camera configuration
    if custom_camera is not None:
        camera_cfg = custom_camera
    else:
        camera_cfg = CAMERA_PRESETS.get(camera_preset.lower(), CAMERA_PRESETS["elife"])

    yaxis_cfg = dict(
        title=dict(text=y_axis_title, font=dict(size=13, color="#333333")),
        showbackground=False,
        gridcolor="#e5e5e5",
        zerolinecolor="#d0d0d0",
    )
    if is_transformed and rescaled_tickvals:
        yaxis_cfg["tickmode"] = "array"
        yaxis_cfg["tickvals"] = rescaled_tickvals
        yaxis_cfg["ticktext"] = rescaled_ticktext

    # Scene and Camera Configuration
    fig.update_layout(
        title=dict(
            text=plot_title,
            x=0.5,
            xanchor="center",
            font=dict(size=18, color="#222222"),
        ),
        paper_bgcolor=paper_bg,
        plot_bgcolor=plot_bg,
        scene=dict(
            xaxis=dict(
                title=dict(text="Tree layout", font=dict(size=13, color="#333333")),
                showbackground=False,
                gridcolor="#e5e5e5",
                zerolinecolor="#d0d0d0",
            ),
            yaxis=yaxis_cfg,
            zaxis=dict(
                title=dict(text="Time before present", font=dict(size=13, color="#333333")),
                showbackground=False,
                gridcolor="#e5e5e5",
                zerolinecolor="#d0d0d0",
            ),
            aspectmode="manual",
            aspectratio=ratio_dict,
            camera=camera_cfg,
        ),
        margin=dict(l=20, r=20, t=50, b=20),
    )

    return fig
