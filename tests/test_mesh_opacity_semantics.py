"""Feature C: mesh opacity semantics (standard opacity, not transparency).

Contract under test (existing ``--opacity`` / ``mesh_opacity`` engine, reused):
- ``--opacity`` is standard opacity: 1.0 = fully opaque, 0.9 = 90% opaque,
  0.7 = 70% opaque, 0.0 = fully transparent.
- Opacity applies ONLY to the curtain Mesh3d trace.
- Mesh geometry, triangles, intensity, color domain and colorscale never change.
- Species labels, hover traces and centerline are not affected.
"""

from __future__ import annotations

import pytest

from phylo3d_trait.renderer import build_figure, build_plot_data
from phylo3d_trait.tree import compute_stable_node_id, parse_tree

TREE_STR = "((A:10,B:10):20,(C:15,D:15):15);"


def _plot_data(**kwargs):
    tree = parse_tree(TREE_STR)
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
    return build_plot_data(tree, trait_values, **kwargs)


def _mesh(fig):
    traces = [t for t in fig.data if t.name == "Branch Curtains"]
    assert len(traces) == 1
    return traces[0]


def test_c1_opacity_1_is_fully_opaque():
    """TEST C1: opacity=1.0 -> Mesh3d opacity == 1.0 (0% transparent)."""
    fig = build_figure(_plot_data(), mesh_opacity=1.0)
    assert _mesh(fig).opacity == pytest.approx(1.0)


def test_c2_opacity_09_is_90_percent_opaque():
    """TEST C2: opacity=0.9 -> 90% opaque / 10% transparent."""
    fig = build_figure(_plot_data(), mesh_opacity=0.9)
    assert _mesh(fig).opacity == pytest.approx(0.9)


def test_c3_opacity_07_is_70_percent_opaque():
    """TEST C3: opacity=0.7 -> 70% opaque / 30% transparent."""
    fig = build_figure(_plot_data(), mesh_opacity=0.7)
    assert _mesh(fig).opacity == pytest.approx(0.7)


def test_c4_mesh_geometry_unchanged_by_opacity():
    """TEST C4: opacity changes nothing except mesh.opacity."""
    plot_data = _plot_data()
    figs = {
        opacity: build_figure(plot_data, mesh_opacity=opacity)
        for opacity in (1.0, 0.9, 0.7)
    }
    meshes = {opacity: _mesh(fig) for opacity, fig in figs.items()}

    baseline = meshes[1.0].to_plotly_json()
    baseline.pop("opacity", None)
    for opacity, mesh in meshes.items():
        if opacity == 1.0:
            continue
        payload = mesh.to_plotly_json()
        assert payload["opacity"] == pytest.approx(opacity)
        payload.pop("opacity", None)
        assert payload == baseline

    for opacity, mesh in meshes.items():
        assert list(mesh.x) == list(meshes[1.0].x)
        assert list(mesh.y) == list(meshes[1.0].y)
        assert list(mesh.z) == list(meshes[1.0].z)
        assert list(mesh.i) == list(meshes[1.0].i)
        assert list(mesh.j) == list(meshes[1.0].j)
        assert list(mesh.k) == list(meshes[1.0].k)
        assert list(mesh.intensity) == list(meshes[1.0].intensity)
        assert mesh.cmin == pytest.approx(meshes[1.0].cmin)
        assert mesh.cmax == pytest.approx(meshes[1.0].cmax)
        assert mesh.colorscale == meshes[1.0].colorscale


def test_c5_non_mesh_traces_unaffected_by_opacity():
    """TEST C5: labels, hover traces and centerline keep their own opacity."""
    plot_data = _plot_data()
    figs = {
        opacity: build_figure(plot_data, mesh_opacity=opacity)
        for opacity in (1.0, 0.9, 0.7)
    }

    def trace(fig, name):
        return [t for t in fig.data if t.name == name][0]

    baseline_labels = trace(figs[1.0], "Species Labels").to_plotly_json()
    baseline_tips = trace(figs[1.0], "Terminal Taxa").to_plotly_json()
    baseline_internal = trace(figs[1.0], "Internal Nodes").to_plotly_json()
    baseline_line = trace(figs[1.0], "Branch Centerlines").to_plotly_json()

    for opacity, fig in figs.items():
        assert trace(fig, "Species Labels").to_plotly_json() == baseline_labels
        assert trace(fig, "Terminal Taxa").to_plotly_json() == baseline_tips
        assert trace(fig, "Internal Nodes").to_plotly_json() == baseline_internal
        assert trace(fig, "Branch Centerlines").to_plotly_json() == baseline_line

    # Hover markers remain invisible-but-hoverable regardless of mesh opacity.
    for fig in figs.values():
        assert trace(fig, "Terminal Taxa").marker.opacity == pytest.approx(0.0)
        assert trace(fig, "Internal Nodes").marker.opacity == pytest.approx(0.0)
