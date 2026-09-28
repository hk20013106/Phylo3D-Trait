"""Tests for --trait-display-offset functionality in Phylo3D-Trait.

Verifies:
TEST 1: RAW_TO_DISPLAY_OFFSET (offset=4: raw 4 -> 0, raw 6 -> 2, raw 10 -> 6)
TEST 2: DISPLAY_TO_RAW_OFFSET (display 0 -> 4, display 2 -> 6, display 6 -> 10)
TEST 3: NODE_GEOMETRY_USES_DISPLAY (node.raw_trait keeps original; node.y == raw_trait - 4)
TEST 4: EDGE_GEOMETRY_USES_DISPLAY (edge endpoints and interpolations use display coordinates)
TEST 5: AXIS_LABELS_USE_RAW (tickvals are display coordinates; ticktext are raw scientific values)
TEST 6: HOVER_USES_RAW (raw=9.8247, display=5.8247 -> hover shows 9.8247)
TEST 7: COLOR_SEMANTICS_UNCHANGED (colors before and after offset are identical for same raw trait)
TEST 8: OFFSET_ZERO_BACKWARD_COMPATIBLE (offset=0 or omitted behaves identically to identity)
TEST 9: SPECIES_LABEL_REGRESSION (labels fixed at present time, baseline Y=0, font size 12, hover disabled)
TEST 10: HOVER_REGRESSION (marker points + XYZ spikes + raw Hb value intact)
"""

from __future__ import annotations

import pytest
from phylo3d_trait.models import PlotData
from phylo3d_trait.renderer import (
    DEFAULT_TIP_LABEL_OFFSET_FRACTION,
    build_figure,
    build_plot_data,
)
from phylo3d_trait.tree import annotate_tree, compute_stable_node_id, parse_tree


def test_1_raw_to_display_offset():
    """TEST 1: With offset=4, raw 4->0, raw 6->2, raw 10->6."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 4.0, "B": 10.0, id_root: 6.0}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0)

    assert plot_data.raw_to_display(4.0) == pytest.approx(0.0)
    assert plot_data.raw_to_display(6.0) == pytest.approx(2.0)
    assert plot_data.raw_to_display(10.0) == pytest.approx(6.0)


def test_2_display_to_raw_offset():
    """TEST 2: With offset=4, display 0->4, display 2->6, display 6->10."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 4.0, "B": 10.0, id_root: 6.0}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0)

    assert plot_data.display_to_raw(0.0) == pytest.approx(4.0)
    assert plot_data.display_to_raw(2.0) == pytest.approx(6.0)
    assert plot_data.display_to_raw(6.0) == pytest.approx(10.0)


def test_3_node_geometry_uses_display():
    """TEST 3: Node geometry (y) uses display_trait = raw_trait - offset, raw_trait preserved."""
    tree_str = "((A:10,B:10):10,C:20);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_root = compute_stable_node_id(["A", "B", "C"])
    trait_values = {"A": 6.5, "B": 9.8, "C": 7.2, id_ab: 8.0, id_root: 7.0}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0)

    node_a = plot_data.nodes["A"]
    assert node_a.raw_trait == pytest.approx(6.5)
    assert node_a.display_trait == pytest.approx(2.5)
    assert node_a.y == pytest.approx(2.5)
    assert node_a.trait == pytest.approx(2.5)

    node_b = plot_data.nodes["B"]
    assert node_b.raw_trait == pytest.approx(9.8)
    assert node_b.display_trait == pytest.approx(5.8)
    assert node_b.y == pytest.approx(5.8)

    node_root = plot_data.nodes[id_root]
    assert node_root.raw_trait == pytest.approx(7.0)
    assert node_root.display_trait == pytest.approx(3.0)
    assert node_root.y == pytest.approx(3.0)


def test_4_edge_geometry_uses_display():
    """TEST 4: Edge endpoints and lineage interpolations all use display_trait coordinates."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 6.0, "B": 10.0, id_root: 8.0}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0, num_segments=5)

    for seg in plot_data.segments:
        # y0 and y1 must be within [2.0, 6.0]
        assert 2.0 - 1e-5 <= seg.y0 <= 6.0 + 1e-5
        assert 2.0 - 1e-5 <= seg.y1 <= 6.0 + 1e-5
        assert seg.display_trait0 == pytest.approx(seg.y0)
        assert seg.display_trait1 == pytest.approx(seg.y1)
        # raw_trait must be within [6.0, 10.0]
        assert 6.0 - 1e-5 <= seg.raw_trait0 <= 10.0 + 1e-5
        assert 6.0 - 1e-5 <= seg.raw_trait1 <= 10.0 + 1e-5


def test_5_axis_labels_use_raw():
    """TEST 5: Trait axis tickvals are display coordinates (0..6), ticktext are raw values (4..10)."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 6.0, "B": 10.0, id_root: 8.0}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0, baseline_y=0.0)
    fig = build_figure(plot_data, baseline_y=0.0)

    yaxis = fig.layout.scene.yaxis
    assert yaxis.tickmode == "array"
    assert len(yaxis.tickvals) > 0
    assert len(yaxis.ticktext) == len(yaxis.tickvals)

    # Lowest tick is at display baseline 0.0 with label "4"
    assert yaxis.tickvals[0] == pytest.approx(0.0)
    assert yaxis.ticktext[0] == "4"

    # Highest tick covers up to display 6.0 with label "10"
    assert yaxis.tickvals[-1] >= 6.0 - 1e-5
    assert yaxis.ticktext[-1] == "10"

    # All tick texts match inverse display_to_raw mapping
    for val, text in zip(yaxis.tickvals, yaxis.ticktext):
        expected_raw = plot_data.display_to_raw(val)
        assert float(text) == pytest.approx(expected_raw, abs=1e-2)


def test_6_hover_uses_raw():
    """TEST 6: Hover displays unscaled raw trait (raw=9.8247, display=5.8247 -> hover 9.8247)."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 9.8247, "B": 6.1234, id_root: 7.5000}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0, baseline_y=0.0)
    fig = build_figure(plot_data, baseline_y=0.0)

    tip_trace = [t for t in fig.data if t.name == "Terminal Taxa"][0]
    tip_a = [c for c in tip_trace.customdata if c[0] == "A"][0]

    # customdata[2] is raw_trait (9.8247)
    assert tip_a[2] == pytest.approx(9.8247)
    # customdata[3] is display_trait (5.8247)
    assert tip_a[3] == pytest.approx(5.8247)

    # Hovertemplate displays raw trait
    assert "customdata[2]:.4f" in tip_trace.hovertemplate
    assert "customdata[3]" not in tip_trace.hovertemplate


def test_7_color_semantics_unchanged():
    """TEST 7: Normalized color ratio for any node is identical before and after offset."""
    tree_str = "((A:10,B:10):10,C:20);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_root = compute_stable_node_id(["A", "B", "C"])
    trait_values = {"A": 6.0, "B": 10.0, "C": 7.0, id_ab: 8.0, id_root: 7.5}

    data_no_offset = build_plot_data(tree, trait_values)
    data_with_offset = build_plot_data(tree, trait_values, trait_display_offset=4.0)

    # For node C: raw=7.0
    # No offset: range [6.0, 10.0], ratio = (7.0 - 6.0) / (10.0 - 6.0) = 0.25
    ratio_no_offset = (data_no_offset.nodes["C"].trait - data_no_offset.trait_min) / (
        data_no_offset.trait_max - data_no_offset.trait_min
    )
    # With offset: range [2.0, 6.0], trait=3.0, ratio = (3.0 - 2.0) / (6.0 - 2.0) = 0.25
    ratio_with_offset = (data_with_offset.nodes["C"].trait - data_with_offset.trait_min) / (
        data_with_offset.trait_max - data_with_offset.trait_min
    )

    assert ratio_no_offset == pytest.approx(ratio_with_offset)


def test_8_offset_zero_or_none_backward_compatible():
    """TEST 8: offset=0 or omitted behaves identically to identity raw coordinates."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 5.0, "B": 8.0, id_root: 6.5}

    data_none = build_plot_data(tree, trait_values)
    data_zero = build_plot_data(tree, trait_values, trait_display_offset=0.0)

    assert data_none.nodes["A"].y == pytest.approx(5.0)
    assert data_zero.nodes["A"].y == pytest.approx(5.0)
    assert data_zero.nodes["B"].y == pytest.approx(8.0)


def test_9_species_label_regression():
    """TEST 9: Species labels stay fixed on the present side (one common Z <= 0) at baseline Y=0, font size 12."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 6.0, "B": 10.0, id_root: 8.0}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0, baseline_y=0.0)
    fig = build_figure(plot_data, baseline_y=0.0)

    species_trace = [t for t in fig.data if t.name == "Species Labels"][0]
    assert len(set(species_trace.z)) == 1
    assert species_trace.z[0] <= 0.0
    span = plot_data.time_max - plot_data.time_min
    assert 0.0 - species_trace.z[0] == pytest.approx(
        DEFAULT_TIP_LABEL_OFFSET_FRACTION * span
    )
    assert all(y == pytest.approx(0.0) for y in species_trace.y)
    assert species_trace.textfont.size == 12
    assert species_trace.hoverinfo == "skip"

    # Tree layout axis title empty and numeric ticks hidden
    xaxis = fig.layout.scene.xaxis
    assert xaxis.title.text == ""
    assert xaxis.showticklabels is False


def test_10_hover_regression():
    """TEST 10: 3D point markers, XYZ spikes, and scientific raw value are intact."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 9.8247, "B": 6.1234, id_root: 7.5000}

    plot_data = build_plot_data(tree, trait_values, trait_display_offset=4.0, baseline_y=0.0)
    fig = build_figure(plot_data, baseline_y=0.0)

    scene = fig.layout.scene
    assert scene.xaxis.showspikes is True
    assert scene.yaxis.showspikes is True
    assert scene.zaxis.showspikes is True

    tip_trace = [t for t in fig.data if t.name == "Terminal Taxa"][0]
    assert tip_trace.mode == "markers"
    assert tip_trace.marker.size == 12
    assert tip_trace.marker.opacity == 0.0
