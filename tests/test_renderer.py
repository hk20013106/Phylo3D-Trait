"""Unit tests for 3D Plotly rendering, clean presentation styling, and eLife camera presets."""

import pytest
import plotly.graph_objects as go

from phylo3d_trait.renderer import build_figure, build_plot_data
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


def test_no_node_or_tip_markers_by_default():
    """Verify no internal node markers or tip circle markers are drawn by default."""
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

    # Verify no marker traces exist
    for trace in fig.data:
        if hasattr(trace, "mode") and trace.mode:
            assert "markers" not in trace.mode, f"Found unexpected marker mode: {trace.mode}"

    # Verify tip trace mode is strictly 'text'
    tip_traces = [t for t in fig.data if t.name == "Terminal Taxa"]
    assert len(tip_traces) == 1
    assert tip_traces[0].mode == "text"


def test_optional_internal_node_markers():
    """Verify that ancestral node markers are added only when show_node_markers=True."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, show_node_markers=True)

    internal_traces = [t for t in fig.data if t.name == "Internal Nodes"]
    assert len(internal_traces) == 1
    assert internal_traces[0].mode == "markers"


def test_scene_axes_and_clean_background():
    """Verify no gray background walls and paper background is pure white."""
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

    assert scene.xaxis.title.text == "Tree layout"
    assert scene.yaxis.title.text == "Trait value"
    assert scene.zaxis.title.text == "Time before present"
    assert scene.aspectmode == "manual"


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


def test_tip_labels_aligned_to_top_front_reference_line():
    """Verify all tip labels are placed at y=global_trait_max and z=0.0 while retaining tip.x."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.2,
        "B": 2.5,
        "C": 3.8,
        "D": 8.0,  # Global maximum trait
        id_ab: 1.8,
        id_cd: 5.0,
        id_root: 3.0,
    }

    plot_data = build_plot_data(tree, trait_values)
    global_max = plot_data.trait_max
    assert global_max == pytest.approx(8.0)

    fig = build_figure(plot_data)

    tip_traces = [t for t in fig.data if t.name == "Terminal Taxa"]
    assert len(tip_traces) == 1
    trace = tip_traces[0]

    # 1. All tip labels y must equal global_trait_max (8.0)
    assert len(trace.y) == 4
    for y_val in trace.y:
        assert y_val == pytest.approx(8.0)

    # 2. All tip labels z must equal 0.0
    assert len(trace.z) == 4
    for z_val in trace.z:
        assert z_val == pytest.approx(0.0)

    # 3. Each label x must match corresponding tip x
    tip_nodes = [n for n in plot_data.nodes.values() if n.is_tip]
    expected_xs = [n.x for n in tip_nodes]
    assert list(trace.x) == expected_xs

    # 4. Biological tip node coordinates in plot_data remain unchanged
    assert plot_data.nodes["A"].y == pytest.approx(1.2)
    assert plot_data.nodes["B"].y == pytest.approx(2.5)
    assert plot_data.nodes["C"].y == pytest.approx(3.8)
    assert plot_data.nodes["D"].y == pytest.approx(8.0)
