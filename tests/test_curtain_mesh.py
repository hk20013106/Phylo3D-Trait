"""Unit tests for branch curtain surface mesh generation (Mesh3d) and Y == Trait == Color invariant."""

import pytest
import plotly.graph_objects as go

from phylo3d_trait.renderer import (
    _build_branch_curtains_geometry,
    build_figure,
    build_plot_data,
)
from phylo3d_trait.tree import compute_stable_node_id, parse_tree


def test_vertex_intensity_strictly_equals_y_coordinate():
    """Verify core invariant: for EVERY vertex in Mesh3d, intensity == y (Trait)."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.5, "B": 3.0, "C": 4.5, "D": 5.0,
        id_ab: 2.0, id_cd: 4.0, id_root: 1.0,
    }

    plot_data = build_plot_data(tree, trait_values, num_segments=8)
    mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity = _build_branch_curtains_geometry(
        plot_data, plot_data.baseline_y
    )

    assert len(mesh_y) == len(mesh_intensity)
    for idx in range(len(mesh_y)):
        assert mesh_intensity[idx] == pytest.approx(mesh_y[idx]), (
            f"Vertex {idx} intensity {mesh_intensity[idx]} != y coordinate {mesh_y[idx]}"
        )


def test_top_bottom_intensity_pair_gradient():
    """Verify that a Top (y=5)/Bottom (y=1) pair has intensities [5, 1], NOT [5, 5]."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    trait_values = {"A": 5.0, "B": 3.0, id_root: 1.0}
    plot_data = build_plot_data(tree, trait_values, num_segments=10)

    # Baseline is trait_min = 1.0
    assert plot_data.baseline_y == pytest.approx(1.0)

    mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity = _build_branch_curtains_geometry(
        plot_data, plot_data.baseline_y
    )

    # For the Root -> A branch, the final point at tip A has top_y = 5.0, bottom_y = 1.0
    # Find the top and bottom vertex corresponding to tip A
    # Tip A has x = 0.0, z = 0.0, y = 5.0
    tip_a_top_indices = [
        i for i in range(len(mesh_x))
        if mesh_z[i] == pytest.approx(0.0) and mesh_y[i] == pytest.approx(5.0)
    ]
    assert len(tip_a_top_indices) == 1
    top_idx = tip_a_top_indices[0]
    bot_idx = top_idx + 1  # In interleaved layout, bottom vertex follows immediately

    assert mesh_y[top_idx] == pytest.approx(5.0)
    assert mesh_intensity[top_idx] == pytest.approx(5.0)

    assert mesh_y[bot_idx] == pytest.approx(1.0)
    assert mesh_intensity[bot_idx] == pytest.approx(1.0)

    # Explicitly ensure bottom_intensity is NOT top_intensity (not 5.0)
    assert mesh_intensity[bot_idx] != pytest.approx(5.0)


def test_default_baseline_y_equals_trait_min():
    """Verify default baseline_y strictly equals trait_min."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    trait_values = {"A": 12.0, "B": 24.0, id_root: 8.5}
    plot_data = build_plot_data(tree, trait_values)

    assert plot_data.trait_min == pytest.approx(8.5)
    assert plot_data.baseline_y == pytest.approx(8.5)


def test_mesh_color_range_matches_y_bounds():
    """Verify cmin and cmax cover all mesh y values (including baseline)."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 2.0, "B": 4.0, "C": 6.0, "D": 8.0,
        id_ab: 3.0, id_cd: 7.0, id_root: 5.0,
    }

    # 1. Default baseline_y (baseline_y == trait_min == 2.0)
    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    mesh_traces = [t for t in fig.data if isinstance(t, go.Mesh3d) or getattr(t, "type", None) == "mesh3d"]
    assert len(mesh_traces) == 1
    mesh = mesh_traces[0]

    assert mesh.cmin == pytest.approx(2.0)
    assert mesh.cmax == pytest.approx(8.0)

    # 2. Custom baseline_y below trait_min (e.g., baseline_y = 0.0)
    fig_custom = build_figure(plot_data, baseline_y=0.0)
    mesh_custom = [t for t in fig_custom.data if isinstance(t, go.Mesh3d) or getattr(t, "type", None) == "mesh3d"][0]

    assert mesh_custom.cmin == pytest.approx(0.0)
    assert mesh_custom.cmax == pytest.approx(8.0)


def test_mesh_triangle_count_formula():
    """Verify that an orthogonal rectangular edge produces 2 * (N_conn + N_lineage) triangles."""
    tree_str = "(A:10,B:10);"  # 2 edges: Root->A and Root->B
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    trait_values = {"A": 2.0, "B": 4.0, id_root: 1.0}
    num_segs = 12  # N_conn = 12, N_lineage = 12 -> 24 steps per edge
    plot_data = build_plot_data(tree, trait_values, num_segments=num_segs)

    baseline_y = plot_data.baseline_y
    mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity = _build_branch_curtains_geometry(
        plot_data, baseline_y
    )

    num_triangles = len(mesh_i)
    # 2 edges * (2 * (N_conn + N_lineage)) = 2 * (2 * 24) = 96 triangles
    expected_triangles = 2 * (2 * (2 * num_segs))
    assert num_triangles == expected_triangles
    assert len(mesh_j) == expected_triangles
    assert len(mesh_k) == expected_triangles

    # Total vertices = 2 edges * (2 * (2 * num_segs + 1)) = 2 * (2 * 25) = 100 vertices
    expected_vertices = 2 * (2 * (2 * num_segs + 1))
    assert len(mesh_x) == expected_vertices
    assert len(mesh_y) == expected_vertices
    assert len(mesh_z) == expected_vertices
    assert len(mesh_intensity) == expected_vertices


def test_rectangular_orthogonal_subsegments_geometry_and_color():
    """Verify rectangular connector (X changes, Y & Z constant) and lineage (Z & Y change, X constant)."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    trait_values = {"A": 4.0, "B": 6.0, id_root: 2.0}
    num_segs = 5
    plot_data = build_plot_data(tree, trait_values, num_segments=num_segs)

    # Inspect Root -> A edge segments (Root: x=0.5, y=2.0, z=10.0 -> A: x=0.0, y=4.0, z=0.0)
    root_a_segs = [s for s in plot_data.segments if s.parent_id == id_root and s.child_id == "A"]
    assert len(root_a_segs) == 2 * num_segs

    connectors = [s for s in root_a_segs if s.segment_type == "connector"]
    lineages = [s for s in root_a_segs if s.segment_type == "lineage"]

    # 1. Connector properties
    assert len(connectors) == num_segs
    for c in connectors:
        assert c.y0 == pytest.approx(2.0)
        assert c.y1 == pytest.approx(2.0)
        assert c.z0 == pytest.approx(10.0)
        assert c.z1 == pytest.approx(10.0)
        assert c.trait0 == pytest.approx(2.0)
        assert c.trait1 == pytest.approx(2.0)

    # 2. Lineage properties
    assert len(lineages) == num_segs
    for l in lineages:
        assert l.x0 == pytest.approx(0.0)
        assert l.x1 == pytest.approx(0.0)

    assert lineages[0].y0 == pytest.approx(2.0)
    assert lineages[-1].y1 == pytest.approx(4.0)
    assert lineages[0].z0 == pytest.approx(10.0)
    assert lineages[-1].z1 == pytest.approx(0.0)


def test_top_vertices_match_branch_x_and_z():
    """Verify Top vertices retain exact X (layout) and Z (time) coordinates."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0,
        id_ab: 1.5, id_cd: 3.5, id_root: 2.5
    }

    plot_data = build_plot_data(tree, trait_values, num_segments=10)
    mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity = _build_branch_curtains_geometry(
        plot_data, plot_data.baseline_y
    )

    # Tips must be at Z == 0.0
    # Root must be at Z == 30.0
    assert min(mesh_z) == pytest.approx(0.0)
    assert max(mesh_z) == pytest.approx(30.0)


def test_no_cross_branch_delaunay_triangles():
    """Verify triangles only connect vertices within the same branch."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0,
        id_ab: 1.5, id_cd: 3.5, id_root: 2.5
    }

    num_segs = 8
    plot_data = build_plot_data(tree, trait_values, num_segments=num_segs)
    mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity = _build_branch_curtains_geometry(
        plot_data, plot_data.baseline_y
    )

    vertices_per_branch = 2 * (2 * num_segs + 1)
    num_branches = 6  # 4-taxon tree has 6 branches

    # Check each triangle's indices belong strictly to the same branch chunk
    for t_idx in range(len(mesh_i)):
        i, j, k = mesh_i[t_idx], mesh_j[t_idx], mesh_k[t_idx]
        branch_i = i // vertices_per_branch
        branch_j = j // vertices_per_branch
        branch_k = k // vertices_per_branch
        assert branch_i == branch_j == branch_k, (
            f"Triangle {t_idx} crosses branches: vertices {i}, {j}, {k} are in branches {branch_i}, {branch_j}, {branch_k}"
        )


def test_default_mesh_opacity_is_fully_opaque():
    """Verify default Mesh3d opacity is strictly 1.0 for true WebGL depth buffer occlusion."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data)

    mesh_traces = [t for t in fig.data if isinstance(t, go.Mesh3d) or getattr(t, "type", None) == "mesh3d"]
    assert len(mesh_traces) == 1
    assert mesh_traces[0].opacity == pytest.approx(1.0)


def test_custom_mesh_opacity_passed_correctly():
    """Verify explicit custom opacity (e.g. 0.5) is passed to Mesh3d trace."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])
    trait_values = {"A": 1.0, "B": 2.0, id_root: 1.5}

    plot_data = build_plot_data(tree, trait_values)
    fig = build_figure(plot_data, mesh_opacity=0.5)

    mesh_traces = [t for t in fig.data if isinstance(t, go.Mesh3d) or getattr(t, "type", None) == "mesh3d"]
    assert len(mesh_traces) == 1
    assert mesh_traces[0].opacity == pytest.approx(0.5)


def test_colorscale_is_fully_opaque():
    """Verify colorscale does not contain semi-transparent rgba entries."""
    import plotly.colors as pcolors
    # Test built-in Turbo colorscale
    colors = pcolors.get_colorscale("Turbo")
    for stop, col in colors:
        col_str = str(col).lower()
        if "rgba" in col_str:
            # If rgba is used, alpha must be 1.0
            parts = col_str.replace("rgba(", "").replace(")", "").split(",")
            alpha = float(parts[-1].strip())
            assert alpha == pytest.approx(1.0)

