"""Unit tests for 3D Plotly rendering, clean presentation styling, and eLife camera presets."""

import pytest
import plotly.graph_objects as go

from phylo3d_trait.renderer import (
    DEFAULT_TIP_LABEL_OFFSET_FRACTION,
    build_figure,
    build_plot_data,
)
from phylo3d_trait.tree import compute_stable_node_id, parse_tree


def test_global_color_normalization():
    """Verify global trait_min and trait_max are uniformly applied across Mesh3d and lines."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)

    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    # Trait values ranging from -5.0 to 12.5
    trait_values = {
        "A": -5.0,
        "B": 0.0,
        "C": 8.0,
        "D": 12.5,
        id_ab: -2.0,
        id_cd: 10.0,
        id_root: 3.5,
    }

    plot_data = build_plot_data(tree, trait_values)
    assert plot_data.trait_min == pytest.approx(-5.0)
    assert plot_data.trait_max == pytest.approx(12.5)

    fig = build_figure(plot_data)
    assert isinstance(fig, go.Figure)

    # Check color limits across traces
    for trace in fig.data:
        if isinstance(trace, go.Mesh3d) or getattr(trace, "type", None) == "mesh3d":
            assert trace.cmin == pytest.approx(-5.0)
            assert trace.cmax == pytest.approx(12.5)
        elif hasattr(trace, "mode") and trace.mode:
            if trace.mode == "lines" and trace.line and hasattr(trace.line, "cmin") and trace.line.cmin is not None:
                assert trace.line.cmin == pytest.approx(-5.0)
                assert trace.line.cmax == pytest.approx(12.5)


def test_reverse_colorscale_changes_color_mapping_only():
    """Reverse colorscale must not alter trait heights or raw scientific values."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 5.0, "B": 10.0, id_root: 7.5}

    plot_data = build_plot_data(tree, trait_values)
    original_y = {node_id: node.y for node_id, node in plot_data.nodes.items()}
    original_raw = {node_id: node.raw_trait for node_id, node in plot_data.nodes.items()}

    fig = build_figure(
        plot_data,
        reverse_colorscale=True,
        centerline_color="trait",
        show_node_markers=True,
    )

    mesh = [t for t in fig.data if getattr(t, "type", None) == "mesh3d"][0]
    line = [t for t in fig.data if t.name == "Branch Centerlines"][0]
    markers = [t for t in fig.data if t.name == "Internal Nodes"][0]

    assert mesh.reversescale is True
    assert line.line.reversescale is True
    assert markers.marker.reversescale is True

    assert {node_id: node.y for node_id, node in plot_data.nodes.items()} == original_y
    assert {node_id: node.raw_trait for node_id, node in plot_data.nodes.items()} == original_raw
    assert plot_data.nodes["A"].y == pytest.approx(5.0)
    assert plot_data.nodes["B"].y == pytest.approx(10.0)


def test_default_colorscale_direction_is_unchanged():
    """Existing plots remain backward-compatible unless reversal is requested."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 5.0, "B": 10.0, id_root: 7.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)
    mesh = [t for t in fig.data if getattr(t, "type", None) == "mesh3d"][0]

    assert mesh.reversescale is False


def test_default_hover_markers_unobtrusive():
    """Verify that node markers are invisible (opacity=0.0) by default but present for 3D raycast hover."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0,
        id_ab: 1.5, id_cd: 3.5, id_root: 2.5
    }
    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, show_node_markers=False)

    # Verify internal node markers are opacity=0.0
    internal_traces = [t for t in fig.data if t.name == "Internal Nodes"]
    assert len(internal_traces) == 1
    assert internal_traces[0].marker.opacity == 0.0

    # Verify terminal taxa hover markers are opacity=0.0
    tip_traces = [t for t in fig.data if t.name == "Terminal Taxa"]
    assert len(tip_traces) == 1
    assert tip_traces[0].marker.opacity == 0.0

    # Verify annotations are empty and axis tick labels are hidden
    assert len(fig.layout.scene.annotations) == 0
    xaxis = fig.layout.scene.xaxis
    assert xaxis.showticklabels is False
    assert len(xaxis.tickvals) == 0


def test_optional_internal_node_markers():
    """Verify that ancestral node markers become visible (opacity=0.95) when show_node_markers=True."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, show_node_markers=True)

    internal_traces = [t for t in fig.data if t.name == "Internal Nodes"]
    assert len(internal_traces) == 1
    assert internal_traces[0].mode == "markers"
    assert internal_traces[0].marker.opacity == 0.95


def test_scene_axes_and_clean_background():
    """Verify no gray background walls, paper background is pure white, and 3D spikes enabled."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, background="white")

    assert fig.layout.paper_bgcolor == "white"
    assert fig.layout.plot_bgcolor == "white"

    scene = fig.layout.scene
    assert scene.xaxis.showbackground is False
    assert scene.yaxis.showbackground is False
    assert scene.zaxis.showbackground is False

    # Goal 2: Tree layout title deleted
    assert scene.xaxis.title.text == ""
    assert scene.yaxis.title.text == "Trait value"
    assert scene.zaxis.title.text == "Time before present"
    assert scene.aspectmode == "manual"

    # Goal 4: 3D spikes on all axes
    assert scene.xaxis.showspikes is True
    assert scene.yaxis.showspikes is True
    assert scene.zaxis.showspikes is True


def test_transparent_background_option():
    """Verify transparent background configuration."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, background="transparent")

    assert fig.layout.paper_bgcolor == "rgba(0,0,0,0)"
    assert fig.layout.plot_bgcolor == "rgba(0,0,0,0)"


def test_elife_camera_preset():
    """Verify eLife camera preset maintains vertical Y (Trait), +Z foreground, and orthographic projection."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, camera_preset="elife")

    camera = fig.layout.scene.camera
    assert camera.up.x == 0
    assert camera.up.y == 1  # Trait (Y) is screen-vertical
    assert camera.up.z == 0
    assert camera.eye.z > 0  # +Z side foreground (MRCA)
    assert camera.projection.type == "orthographic"


def test_species_labels_fixed_at_present_side():
    """TEST 1 (Goal 1): All species labels stay fixed on the present side (single common Z <= time_min)."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.2,
        "B": 2.5,
        "C": 3.8,
        "D": 8.0,
        id_ab: 1.8,
        id_cd: 5.0,
        id_root: 3.0,
    }

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    label_traces = [t for t in fig.data if t.name == "Species Labels"]
    assert len(label_traces) == 1
    label_trace = label_traces[0]

    # All labels share one fixed present-side anchor, pushed outward beyond the
    # present plane by the generic fraction of the time span (never drifting).
    assert len(set(label_trace.z)) == 1
    assert label_trace.z[0] <= plot_data.time_min
    span = plot_data.time_max - plot_data.time_min
    assert plot_data.time_min - label_trace.z[0] == pytest.approx(
        DEFAULT_TIP_LABEL_OFFSET_FRACTION * span
    )
    assert plot_data.time_min == pytest.approx(0.0)


def test_species_label_layout_mapping():
    """TEST 2 (Goal 1): Verify species labels match terminal tips in tree layout order (X coordinates)."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.2,
        "B": 2.5,
        "C": 3.8,
        "D": 8.0,
        id_ab: 1.8,
        id_cd: 5.0,
        id_root: 3.0,
    }

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    label_trace = [t for t in fig.data if t.name == "Species Labels"][0]
    tip_nodes = [n for n in plot_data.nodes.values() if n.is_tip]
    sorted_tips = sorted(tip_nodes, key=lambda n: n.x)

    assert len(label_trace.x) == len(sorted_tips)
    assert list(label_trace.text) == [n.label for n in sorted_tips]
    for i, tip in enumerate(sorted_tips):
        assert label_trace.x[i] == pytest.approx(tip.x)


def test_species_label_constant_trait_plane():
    """TEST 3 (Goal 1): All species labels must lie on the constant baseline Y plane."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.2,
        "B": 2.5,
        "C": 3.8,
        "D": 8.0,
        id_ab: 1.8,
        id_cd: 5.0,
        id_root: 3.0,
    }

    # Test default baseline
    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)
    label_trace = [t for t in fig.data if t.name == "Species Labels"][0]
    assert len(set(label_trace.y)) == 1
    assert label_trace.y[0] == pytest.approx(plot_data.baseline_y)

    # Test custom baseline_y
    fig_custom = build_figure(plot_data, baseline_y=4.0)
    label_trace_custom = [t for t in fig_custom.data if t.name == "Species Labels"][0]
    assert len(set(label_trace_custom.y)) == 1
    assert label_trace_custom.y[0] == pytest.approx(4.0)


def test_tree_layout_axis_title_removed():
    """TEST 4 (Goal 2): Scene X axis title text must be completely empty."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    assert fig.layout.scene.xaxis.title.text == ""


def test_tree_layout_numeric_ticks_hidden():
    """TEST 5 (Goal 2): Numeric tick labels on Tree Layout axis must be hidden."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0,
        id_ab: 1.5, id_cd: 3.5, id_root: 2.5
    }

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    xaxis = fig.layout.scene.xaxis
    assert xaxis.showticklabels is False
    assert len(xaxis.tickvals) == 0
    assert len(xaxis.ticktext) == 0
    assert xaxis.ticks == ""


def test_species_font_scale_130_percent():
    """TEST 6 (Goal 3): Species label font size must be 12 (130% of previous size 9)."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    label_trace = [t for t in fig.data if t.name == "Species Labels"][0]
    assert label_trace.textfont.size == 12


def test_label_trace_hover_disabled():
    """TEST 7 (Goal 4): Species label text trace must have hoverinfo='skip' to avoid interfering with 3D hover."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    label_trace = [t for t in fig.data if t.name == "Species Labels"][0]
    assert label_trace.hoverinfo == "skip"


def test_node_hover_restored():
    """TEST 8 (Goal 4): Tree tips and internal nodes must have dedicated hover traces with raw trait."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 7.42,
        "B": 6.85,
        "C": 8.10,
        "D": 5.90,
        id_ab: 7.10,
        id_cd: 7.00,
        id_root: 7.05,
    }

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    tip_traces = [t for t in fig.data if t.name == "Terminal Taxa"]
    internal_traces = [t for t in fig.data if t.name == "Internal Nodes"]

    assert len(tip_traces) == 1
    assert len(internal_traces) == 1

    # Terminal taxa hover template contains raw trait
    tip_trace = tip_traces[0]
    assert "customdata[2]:.4f" in tip_trace.hovertemplate
    # Check first tip customdata raw trait
    assert tip_trace.customdata[0][2] == pytest.approx(7.42)

    # Spikes enabled on scene
    scene = fig.layout.scene
    assert scene.xaxis.showspikes is True
    assert scene.yaxis.showspikes is True
    assert scene.zaxis.showspikes is True


def test_raw_value_in_hover():
    """TEST 9 (Goal 4): Hover tooltip must display unscaled raw trait even when display transform is active."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 5.0, "B": 10.0, id_root: 7.5}

    # Apply linear scale to display range (13.0, 5.0)
    plot_data = build_plot_data(
        tree,
        trait_values,
        trait_display_range=(13.0, 5.0),
    )

    fig = build_figure(plot_data)
    tip_trace = [t for t in fig.data if t.name == "Terminal Taxa"][0]

    # Raw trait should remain 5.0 and 10.0, while display trait is 13.0 and 5.0
    tip_a = [c for c in tip_trace.customdata if c[0] == "A"][0]
    assert tip_a[2] == pytest.approx(5.0)  # raw trait
    assert tip_a[3] == pytest.approx(13.0)  # display trait

    tip_b = [c for c in tip_trace.customdata if c[0] == "B"][0]
    assert tip_b[2] == pytest.approx(10.0)  # raw trait
    assert tip_b[3] == pytest.approx(5.0)  # display trait


def test_geometry_regression():
    """TEST 10: Label refactoring must not alter any 3D tree coordinates, lines, or mesh geometry."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.2, "B": 2.5, "C": 3.8, "D": 8.0,
        id_ab: 1.8, id_cd: 5.0, id_root: 3.0,
    }

    plot_data = build_plot_data(tree, trait_values)
    fig_labels = build_figure(plot_data, show_tip_labels=True)
    fig_nolabels = build_figure(plot_data, show_tip_labels=False)

    mesh_a = [t for t in fig_labels.data if getattr(t, "type", None) == "mesh3d"]
    mesh_b = [t for t in fig_nolabels.data if getattr(t, "type", None) == "mesh3d"]
    assert len(mesh_a) == len(mesh_b) == 6  # one curtain trace per biological edge

    line_a = [t for t in fig_labels.data if t.name == "Branch Centerlines"][0]
    line_b = [t for t in fig_nolabels.data if t.name == "Branch Centerlines"][0]

    # Mesh geometry strictly identical per edge (visible labels never touch geometry)
    for mesh_edge_a, mesh_edge_b in zip(mesh_a, mesh_b):
        assert list(mesh_edge_a.x) == list(mesh_edge_b.x)
        assert list(mesh_edge_a.y) == list(mesh_edge_b.y)
        assert list(mesh_edge_a.z) == list(mesh_edge_b.z)
        assert list(mesh_edge_a.i) == list(mesh_edge_b.i)
        assert list(mesh_edge_a.j) == list(mesh_edge_b.j)
        assert list(mesh_edge_a.k) == list(mesh_edge_b.k)
        assert list(mesh_edge_a.intensity) == list(mesh_edge_b.intensity)

    # Line geometry strictly identical
    assert list(line_a.x) == list(line_b.x)
    assert list(line_a.y) == list(line_b.y)
    assert list(line_a.z) == list(line_b.z)

    # Node coordinates in plot_data strictly intact
    assert plot_data.nodes["A"].x == pytest.approx(0.0)
    assert plot_data.nodes["B"].x == pytest.approx(1.0)
    assert plot_data.nodes["C"].x == pytest.approx(2.0)
    assert plot_data.nodes["D"].x == pytest.approx(3.0)
    assert plot_data.nodes[id_ab].x == pytest.approx(0.5)
    assert plot_data.nodes[id_cd].x == pytest.approx(2.5)
    assert plot_data.nodes[id_root].x == pytest.approx(1.5)


def test_no_labels_option_hides_species_trace():
    """TEST 11: --no-labels (show_tip_labels=False) hides the Species Labels trace."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0,
        id_ab: 1.5, id_cd: 3.5, id_root: 2.5
    }

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, show_tip_labels=False)

    # Species Labels trace must not be present
    label_traces = [t for t in fig.data if t.name == "Species Labels"]
    assert len(label_traces) == 0

    xaxis = fig.layout.scene.xaxis
    assert xaxis.showticklabels is False
    assert len(xaxis.tickvals) == 0
    assert len(xaxis.ticktext) == 0
    assert len(fig.layout.scene.annotations) == 0


def test_previous_pr1_features_no_regression():
    """TEST 12: PR #1 color reversal, branch curtain mode, and trait centerline have no regression."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0,
        id_ab: 1.5, id_cd: 3.5, id_root: 2.5
    }

    plot_data = build_plot_data(tree, trait_values, baseline_y=0.0)
    fig = build_figure(
        plot_data,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        centerline_color="trait",
        show_tip_labels=True,
    )

    mesh = [t for t in fig.data if getattr(t, "type", None) == "mesh3d"][0]
    line = [t for t in fig.data if t.name == "Branch Centerlines"][0]

    assert mesh.reversescale is True
    assert line.line.reversescale is True
    assert line.line.colorscale is not None

    # Species Labels trace present, fixed on the present side (single common Z <= time_min)
    label_traces = [t for t in fig.data if t.name == "Species Labels"]
    assert len(label_traces) == 1
    assert len(set(label_traces[0].z)) == 1
    assert label_traces[0].z[0] <= plot_data.time_min
    assert list(label_traces[0].text) == ["A", "B", "C", "D"]

