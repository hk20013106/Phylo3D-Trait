"""Unit and integration tests for linear trait rescaling (reverse/forward display range mapping)."""

import pytest
import plotly.graph_objects as go

from phylo3d_trait.models import PlotData
from phylo3d_trait.renderer import build_figure, build_plot_data
from phylo3d_trait.tree import annotate_tree, compute_stable_node_id, parse_tree


def test_default_no_transform_display_trait_equals_raw_trait():
    """Verify default behavior when trait_display_range is None."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.5, "B": 3.0, "C": 4.5, "D": 8.0,
        id_ab: 2.0, id_cd: 6.0, id_root: 3.5
    }

    plot_data = build_plot_data(tree, trait_values)
    assert plot_data.trait_display_range is None
    assert plot_data.raw_trait_min == pytest.approx(1.5)
    assert plot_data.raw_trait_max == pytest.approx(8.0)
    assert plot_data.trait_min == pytest.approx(1.5)
    assert plot_data.trait_max == pytest.approx(8.0)

    for nid, node in plot_data.nodes.items():
        assert node.raw_trait == pytest.approx(trait_values[nid])
        assert node.display_trait == pytest.approx(trait_values[nid])
        assert node.y == pytest.approx(trait_values[nid])
        assert node.trait == pytest.approx(trait_values[nid])


def test_linear_reverse_trait_rescaling_formula():
    """Verify linear mapping [0, 10] -> [13, 5]: display = 13 - 0.8 * raw."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    # Raw trait range strictly [0.0, 10.0]
    trait_values = {
        "A": 0.0,
        "B": 2.5,
        "C": 7.5,
        "D": 10.0,
        id_ab: 1.25,
        id_cd: 8.75,
        id_root: 5.0,
    }

    plot_data = build_plot_data(
        tree,
        trait_values,
        trait_display_range=(13.0, 5.0),
    )

    assert plot_data.trait_display_range == (13.0, 5.0)
    assert plot_data.raw_trait_min == pytest.approx(0.0)
    assert plot_data.raw_trait_max == pytest.approx(10.0)
    assert plot_data.trait_min == pytest.approx(5.0)
    assert plot_data.trait_max == pytest.approx(13.0)

    # Specific point checks: display = 13 - 0.8 * raw
    assert plot_data.nodes["A"].raw_trait == pytest.approx(0.0)
    assert plot_data.nodes["A"].display_trait == pytest.approx(13.0)
    assert plot_data.nodes["A"].y == pytest.approx(13.0)

    assert plot_data.nodes["B"].raw_trait == pytest.approx(2.5)
    assert plot_data.nodes["B"].display_trait == pytest.approx(11.0)
    assert plot_data.nodes["B"].y == pytest.approx(11.0)

    assert plot_data.nodes[id_root].raw_trait == pytest.approx(5.0)
    assert plot_data.nodes[id_root].display_trait == pytest.approx(9.0)
    assert plot_data.nodes[id_root].y == pytest.approx(9.0)

    assert plot_data.nodes["C"].raw_trait == pytest.approx(7.5)
    assert plot_data.nodes["C"].display_trait == pytest.approx(7.0)
    assert plot_data.nodes["C"].y == pytest.approx(7.0)

    assert plot_data.nodes["D"].raw_trait == pytest.approx(10.0)
    assert plot_data.nodes["D"].display_trait == pytest.approx(5.0)
    assert plot_data.nodes["D"].y == pytest.approx(5.0)

    # Linearity check for all nodes
    for node in plot_data.nodes.values():
        expected_disp = 13.0 - 0.8 * node.raw_trait
        assert node.display_trait == pytest.approx(expected_disp)
        assert node.y == pytest.approx(expected_disp)


def test_source_range_from_nodes_only_not_baseline():
    """Verify source min/max is derived solely from node traits, and baseline_y=0 is preserved."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    # Raw range is [2.0, 8.0]
    trait_values = {
        "A": 2.0,
        "B": 8.0,
        id_root: 5.0,
    }

    plot_data = build_plot_data(
        tree,
        trait_values,
        baseline_y=0.0,
        trait_display_range=(13.0, 5.0),
    )

    # Source range must be [2.0, 8.0]
    assert plot_data.raw_trait_min == pytest.approx(2.0)
    assert plot_data.raw_trait_max == pytest.approx(8.0)

    # Mapping: raw 2.0 -> 13.0, raw 8.0 -> 5.0, raw 5.0 -> 9.0
    assert plot_data.nodes["A"].display_trait == pytest.approx(13.0)
    assert plot_data.nodes["B"].display_trait == pytest.approx(5.0)
    assert plot_data.nodes[id_root].display_trait == pytest.approx(9.0)

    # Baseline must remain exactly 0.0
    assert plot_data.baseline_y == pytest.approx(0.0)

    fig = build_figure(plot_data, baseline_y=0.0)

    # Check Mesh3d trace
    mesh_traces = [t for t in fig.data if isinstance(t, go.Mesh3d) or getattr(t, "type", None) == "mesh3d"]
    assert len(mesh_traces) == 1
    mesh = mesh_traces[0]

    # Vertex Y values should have min == 0.0 and max == 13.0
    assert min(mesh.y) == pytest.approx(0.0)
    assert max(mesh.y) == pytest.approx(13.0)

    # Color range should be [0.0, 13.0]
    assert mesh.cmin == pytest.approx(0.0)
    assert mesh.cmax == pytest.approx(13.0)

    # Invariant: mesh vertex intensity strictly equals vertex Y
    for y_val, intensity_val in zip(mesh.y, mesh.intensity):
        assert intensity_val == pytest.approx(y_val)


def test_hover_text_shows_both_raw_and_display_trait():
    """Verify hover text clearly displays both Raw Trait and Display Trait."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    trait_values = {
        "A": 0.0,
        "B": 10.0,
        id_root: 5.0,
    }

    plot_data = build_plot_data(
        tree,
        trait_values,
        trait_display_range=(13.0, 5.0),
    )

    fig = build_figure(plot_data, show_node_markers=True)

    # Tip trace
    tip_trace = [t for t in fig.data if t.name == "Terminal Taxa"][0]
    assert "Raw Trait:" in tip_trace.hovertemplate
    assert "Display Trait (Y):" in tip_trace.hovertemplate

    # Customdata has [label, node_id, raw_trait, display_trait]
    for row in tip_trace.customdata:
        label, nid, raw_t, disp_t = row
        if label == "A":
            assert raw_t == pytest.approx(0.0)
            assert disp_t == pytest.approx(13.0)
        elif label == "B":
            assert raw_t == pytest.approx(10.0)
            assert disp_t == pytest.approx(5.0)

    # Internal node trace
    internal_trace = [t for t in fig.data if t.name == "Internal Nodes"][0]
    assert "Raw Trait:" in internal_trace.hovertemplate
    assert "Display Trait (Y / Height):" in internal_trace.hovertemplate
    for row in internal_trace.customdata:
        nid, ntype, raw_t, disp_t, desc = row
        assert raw_t == pytest.approx(5.0)
        assert disp_t == pytest.approx(9.0)
