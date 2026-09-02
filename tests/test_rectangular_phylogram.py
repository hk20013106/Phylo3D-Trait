"""Comprehensive unit tests for rectangular phylogram geometry and orthogonal curtain mesh."""

import pytest
import plotly.graph_objects as go

from phylo3d_trait.renderer import _build_branch_curtains_geometry, build_figure, build_plot_data
from phylo3d_trait.tree import annotate_tree, compute_stable_node_id, parse_tree


def test_biological_edge_split_into_orthogonal_subsegments():
    """Verify each parent->child biological edge decomposes into connector and lineage subsegments."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)
    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0,
        id_ab: 1.5, id_cd: 3.5, id_root: 2.5,
    }

    num_segs = 6
    plot_data = annotate_tree(tree, trait_values, num_segments=num_segs)

    # 4-taxon tree has 6 biological edges: (root->ab, root->cd, ab->A, ab->B, cd->C, cd->D)
    edge_keys = set((s.parent_id, s.child_id) for s in plot_data.segments)
    assert len(edge_keys) == 6

    for parent_id, child_id in edge_keys:
        edge_segs = [s for s in plot_data.segments if s.parent_id == parent_id and s.child_id == child_id]
        edge_segs.sort(key=lambda s: s.segment_index)

        p_node = plot_data.nodes[parent_id]
        c_node = plot_data.nodes[child_id]

        xp, yp, zp = p_node.x, p_node.y, p_node.z
        xc, yc, zc = c_node.x, c_node.y, c_node.z

        connectors = [s for s in edge_segs if s.segment_type == "connector"]
        lineages = [s for s in edge_segs if s.segment_type == "lineage"]

        if xp != xc:
            # 1. Connector starts at (Xp, Yp, Zp) and ends at elbow (Xc, Yp, Zp)
            assert len(connectors) == num_segs
            assert connectors[0].x0 == pytest.approx(xp)
            assert connectors[-1].x1 == pytest.approx(xc)
            for c in connectors:
                assert c.y0 == pytest.approx(yp)
                assert c.y1 == pytest.approx(yp)
                assert c.z0 == pytest.approx(zp)
                assert c.z1 == pytest.approx(zp)
                assert c.trait0 == pytest.approx(yp)
                assert c.trait1 == pytest.approx(yp)

            # 2. Lineage starts at elbow (Xc, Yp, Zp) and ends at child (Xc, Yc, Zc)
            assert len(lineages) == num_segs
            assert lineages[0].x0 == pytest.approx(xc)
            assert lineages[-1].x1 == pytest.approx(xc)
            assert lineages[0].y0 == pytest.approx(yp)
            assert lineages[-1].y1 == pytest.approx(yc)
            assert lineages[0].z0 == pytest.approx(zp)
            assert lineages[-1].z1 == pytest.approx(zc)
            for l in lineages:
                assert l.x0 == pytest.approx(xc)
                assert l.x1 == pytest.approx(xc)


def test_rectangular_curtain_mesh_vertex_intensities():
    """Verify curtain mesh intensities: connector top is constant Yp, lineage top interpolates to Yc."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    trait_values = {"A": 5.0, "B": 8.0, id_root: 2.0}
    plot_data = build_plot_data(tree, trait_values, num_segments=10)

    # Baseline is trait_min = 2.0
    baseline_y = plot_data.baseline_y
    assert baseline_y == pytest.approx(2.0)

    mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity = _build_branch_curtains_geometry(
        plot_data, baseline_y
    )

    # All mesh vertices must satisfy: intensity == y
    for idx in range(len(mesh_y)):
        assert mesh_intensity[idx] == pytest.approx(mesh_y[idx])

    # Along the connector top edge (where z == 10.0 and y == 2.0), top intensity is strictly 2.0
    connector_top_vertices = [
        mesh_intensity[i] for i in range(len(mesh_x))
        if mesh_z[i] == pytest.approx(10.0) and mesh_y[i] == pytest.approx(2.0)
    ]
    assert len(connector_top_vertices) > 0
    for val in connector_top_vertices:
        assert val == pytest.approx(2.0)

    # At tip A (x=0.0, z=0.0), top vertex is 5.0 and bottom is baseline_y (2.0)
    tip_a_indices = [
        i for i in range(len(mesh_x))
        if mesh_x[i] == pytest.approx(0.0) and mesh_z[i] == pytest.approx(0.0) and mesh_y[i] == pytest.approx(5.0)
    ]
    assert len(tip_a_indices) == 1
    top_idx = tip_a_indices[0]
    bot_idx = top_idx + 1

    assert mesh_y[top_idx] == pytest.approx(5.0)
    assert mesh_intensity[top_idx] == pytest.approx(5.0)

    assert mesh_y[bot_idx] == pytest.approx(2.0)
    assert mesh_intensity[bot_idx] == pytest.approx(2.0)


def test_rectangular_centerline_outline_path():
    """Verify branch top outline follows exact rectangular elbow path (parent -> elbow -> child)."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)
    id_root = compute_stable_node_id(["A", "B"])

    trait_values = {"A": 3.0, "B": 5.0, id_root: 1.0}
    plot_data = build_plot_data(tree, trait_values, num_segments=4)
    fig = build_figure(plot_data)

    line_traces = [t for t in fig.data if t.name == "Branch Centerlines"]
    assert len(line_traces) == 1
    line = line_traces[0]

    # Root is at (x=0.5, y=1.0, z=10.0)
    # Child A is at (x=0.0, y=3.0, z=0.0)
    # Elbow for A is at (x=0.0, y=1.0, z=10.0)
    xs, ys, zs = list(line.x), list(line.y), list(line.z)

    # Find the elbow point in the trace
    has_elbow = any(
        xs[i] == pytest.approx(0.0) and ys[i] == pytest.approx(1.0) and zs[i] == pytest.approx(10.0)
        for i in range(len(xs)) if xs[i] is not None
    )
    assert has_elbow, "Elbow vertex (Xc=0.0, Yp=1.0, Zp=10.0) must exist in rectangular centerline outline!"
