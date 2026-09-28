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

    # Verify tip labels exist on Tree layout axis array ticks and annotations are empty
    annotations = fig.layout.scene.annotations
    assert len(annotations) == 0
    xaxis = fig.layout.scene.xaxis
    assert xaxis.tickmode == "array"
    assert list(xaxis.tickvals) == [0.0, 1.0, 2.0, 3.0]
    assert list(xaxis.ticktext) == ["A", "B", "C", "D"]


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


def test_axis_labels_map_exactly_to_tips():
    """TEST 1: Verify Tree layout axis labels map exactly to terminal tips in layout order."""
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

    xaxis = fig.layout.scene.xaxis
    tip_nodes = [n for n in plot_data.nodes.values() if n.is_tip]
    sorted_tips = sorted(tip_nodes, key=lambda n: n.x)

    assert len(xaxis.tickvals) == len(sorted_tips)
    assert len(xaxis.ticktext) == len(sorted_tips)

    for i, tip in enumerate(sorted_tips):
        assert xaxis.tickvals[i] == pytest.approx(tip.x)
        assert xaxis.ticktext[i] == tip.label


def test_no_internal_nodes_on_axis():
    """TEST 2: Verify internal ancestral nodes are strictly excluded from Tree layout axis ticks."""
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
    internal_ids = {id_ab, id_cd, id_root}

    # No internal IDs in axis ticks
    for label in xaxis.ticktext:
        assert label not in internal_ids
        assert not label.startswith("clade:")

    # Exactly matches terminal tip count
    assert len(xaxis.ticktext) == 4


def test_no_duplicate_terminal_3d_text():
    """TEST 3: Verify default mode has no duplicate floating 3D text traces or annotations."""
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

    # No annotations in scene
    assert len(fig.layout.scene.annotations) == 0

    # No Scatter3d trace with mode="text"
    for trace in fig.data:
        if hasattr(trace, "mode") and trace.mode:
            assert "text" not in trace.mode


def test_geometry_unchanged_by_label_rendering():
    """TEST 4: Label refactoring must not alter any 3D geometry or coordinates."""
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

    mesh_a = [t for t in fig_labels.data if t.name == "Branch Curtains"][0]
    mesh_b = [t for t in fig_nolabels.data if t.name == "Branch Curtains"][0]

    line_a = [t for t in fig_labels.data if t.name == "Branch Centerlines"][0]
    line_b = [t for t in fig_nolabels.data if t.name == "Branch Centerlines"][0]

    # Mesh geometry strictly identical
    assert list(mesh_a.x) == list(mesh_b.x)
    assert list(mesh_a.y) == list(mesh_b.y)
    assert list(mesh_a.z) == list(mesh_b.z)
    assert list(mesh_a.i) == list(mesh_b.i)
    assert list(mesh_a.j) == list(mesh_b.j)
    assert list(mesh_a.k) == list(mesh_b.k)
    assert list(mesh_a.intensity) == list(mesh_b.intensity)

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


def test_no_labels_option_hides_ticks_and_annotations():
    """TEST 5: --no-labels (show_tip_labels=False) hides species ticks and floating annotations."""
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

    xaxis = fig.layout.scene.xaxis
    assert xaxis.showticklabels is False
    assert len(xaxis.tickvals) == 0
    assert len(xaxis.ticktext) == 0
    assert len(fig.layout.scene.annotations) == 0


def test_previous_pr1_features_no_regression():
    """TEST 6: PR #1 color reversal, branch curtain mode, and trait centerline have no regression."""
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

    mesh = [t for t in fig.data if t.name == "Branch Curtains"][0]
    line = [t for t in fig.data if t.name == "Branch Centerlines"][0]

    assert mesh.reversescale is True
    assert line.line.reversescale is True
    assert line.line.colorscale is not None

    # Tree layout axis species labels intact
    xaxis = fig.layout.scene.xaxis
    assert list(xaxis.ticktext) == ["A", "B", "C", "D"]
    assert list(xaxis.tickvals) == [0.0, 1.0, 2.0, 3.0]
    assert len(fig.layout.scene.annotations) == 0
