"""Feature C / Bug #2: mesh opacity semantics and per-edge curtain traces.

Root cause fixed in this round:
- All biological parent->child curtain ribbons used to be merged into ONE
  Mesh3d trace. With opacity < 1.0, WebGL alpha blending cannot sort thousands
  of overlapping triangles inside a single draw call, so fronts blended over
  backs and the result looked far more transparent than requested.
- Curtains are now one Mesh3d trace per biological parent->child edge, while
  every scientific coordinate, triangle, intensity and color-domain value is
  preserved exactly (aggregate equivalence vs the legacy single-trace layout).

Contract under test:
- ``--opacity`` is standard opacity: 1.0 = fully opaque, 0.9 = 90% opaque,
  0.7 = 70% opaque, 0.0 = fully transparent.
- All curtain traces share ONE global color domain and ONE colorbar.
- Opacity applies ONLY to curtain Mesh3d traces.
- 0 / negative / >1 / NaN / inf are rejected loudly.
"""

from __future__ import annotations

import pytest
import plotly.graph_objects as go

from phylo3d_trait.renderer import (
    _build_branch_curtains_geometry,
    build_figure,
    build_plot_data,
)
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


def _curtain_meshes(fig):
    meshes = [t for t in fig.data if getattr(t, "type", None) == "mesh3d"]
    assert meshes, "figure must contain at least one curtain Mesh3d trace"
    return meshes


def _edge_keys(plot_data):
    return {(seg.parent_id, seg.child_id) for seg in plot_data.segments}


def _aggregate_mesh_geometry(meshes):
    """Concatenate per-edge mesh arrays, re-basing triangle indices globally."""
    xs, ys, zs, ii, jj, kk, inten = [], [], [], [], [], [], []
    offset = 0
    for mesh in meshes:
        xs.extend(mesh.x)
        ys.extend(mesh.y)
        zs.extend(mesh.z)
        ii.extend(idx + offset for idx in mesh.i)
        jj.extend(idx + offset for idx in mesh.j)
        kk.extend(idx + offset for idx in mesh.k)
        inten.extend(mesh.intensity)
        offset += len(mesh.x)
    return xs, ys, zs, ii, jj, kk, inten


# --------------------------------------------------------------------------- #
# TEST O1: one Mesh3d per biological parent->child edge                       #
# --------------------------------------------------------------------------- #
def test_o1_one_curtain_mesh_per_biological_edge():
    """TEST O1 / ONE_EDGE_ONE_MESH: curtain trace count == biological edge count."""
    plot_data = _plot_data()
    fig = build_figure(plot_data)
    meshes = _curtain_meshes(fig)

    expected_edges = _edge_keys(plot_data)
    assert len(meshes) == len(expected_edges) == 6

    names = [mesh.name for mesh in meshes]
    assert len(set(names)) == len(names), "each per-edge curtain must be uniquely named"
    assert set(names) == {
        f"Branch Curtain: {parent_id} -> {child_id}"
        for (parent_id, child_id) in expected_edges
    }


# --------------------------------------------------------------------------- #
# TEST O2: opacity propagates to every curtain trace                          #
# --------------------------------------------------------------------------- #
def test_o2_opacity_applied_to_every_curtain_mesh():
    """TEST O2: every per-edge curtain carries the requested standard opacity."""
    plot_data = _plot_data()
    fig = build_figure(plot_data, mesh_opacity=0.9)
    meshes = _curtain_meshes(fig)
    assert len(meshes) == 6
    for mesh in meshes:
        assert mesh.opacity == pytest.approx(0.9)


@pytest.mark.parametrize("opacity", [1.0, 0.9, 0.7, 0.0])
def test_o2_opacity_parameter_passthrough(opacity):
    """1.0 = fully opaque, 0.9 = 90% opaque, 0.7 = 70% opaque, 0.0 = transparent."""
    plot_data = _plot_data()
    fig = build_figure(plot_data, mesh_opacity=opacity)
    for mesh in _curtain_meshes(fig):
        assert mesh.opacity == pytest.approx(opacity)


# --------------------------------------------------------------------------- #
# TEST O3: one shared scientific color domain                                 #
# --------------------------------------------------------------------------- #
def test_o3_all_curtains_share_color_domain_and_scale():
    """TEST O3: identical cmin/cmax/colorscale/reversescale on every curtain."""
    plot_data = _plot_data()
    fig = build_figure(plot_data, mesh_opacity=0.9, reverse_colorscale=True)
    meshes = _curtain_meshes(fig)

    reference = meshes[0]
    for mesh in meshes[1:]:
        assert mesh.cmin == pytest.approx(reference.cmin)
        assert mesh.cmax == pytest.approx(reference.cmax)
        assert mesh.colorscale == reference.colorscale
        assert mesh.reversescale is True


# --------------------------------------------------------------------------- #
# TEST O4: exactly one colorbar                                               #
# --------------------------------------------------------------------------- #
def test_o4_exactly_one_colorbar():
    """TEST O4: only the first curtain trace owns the single shared colorbar."""
    plot_data = _plot_data()
    fig = build_figure(plot_data, mesh_opacity=0.9)
    meshes = _curtain_meshes(fig)

    scaled = [mesh for mesh in meshes if mesh.showscale]
    assert len(scaled) == 1
    assert scaled[0] is meshes[0]
    assert scaled[0].colorbar.title.text == "Trait Value"


# --------------------------------------------------------------------------- #
# TEST O5: aggregate geometry identical to the legacy single-trace layout      #
# --------------------------------------------------------------------------- #
def test_o5_aggregate_geometry_matches_legacy_aggregate():
    """TEST O5: per-edge traces recombine to the exact legacy aggregate geometry."""
    plot_data = _plot_data()
    baseline_y = plot_data.baseline_y
    fig = build_figure(plot_data, mesh_opacity=0.9)
    meshes = _curtain_meshes(fig)

    legacy = _build_branch_curtains_geometry(plot_data, baseline_y)
    combined = _aggregate_mesh_geometry(meshes)

    for got, want in zip(combined, legacy):
        assert list(got) == list(want)

    # Scientific color domain too
    assert meshes[0].cmin == pytest.approx(plot_data.trait_min)
    assert meshes[0].cmax == pytest.approx(plot_data.trait_max)


# --------------------------------------------------------------------------- #
# TEST O6: opacity changes nothing but the opacity property                   #
# --------------------------------------------------------------------------- #
def test_o6_opacity_changes_only_the_opacity_property():
    """TEST O6: 1.0 vs 0.9 keep identical geometry, color domain and metadata."""
    plot_data = _plot_data()
    fig_1 = build_figure(plot_data, mesh_opacity=1.0)
    fig_9 = build_figure(plot_data, mesh_opacity=0.9)

    meshes_1 = _curtain_meshes(fig_1)
    meshes_9 = _curtain_meshes(fig_9)
    assert len(meshes_1) == len(meshes_9)

    for mesh_1, mesh_9 in zip(meshes_1, meshes_9):
        payload_1 = mesh_1.to_plotly_json()
        payload_9 = mesh_9.to_plotly_json()
        assert payload_1["opacity"] == pytest.approx(1.0)
        assert payload_9["opacity"] == pytest.approx(0.9)
        payload_1.pop("opacity")
        payload_9.pop("opacity")
        assert payload_1 == payload_9


# --------------------------------------------------------------------------- #
# TEST O7: branch mode keeps top/bottom intensity pairs equal per edge        #
# --------------------------------------------------------------------------- #
def test_o7_branch_mode_intensity_pairs_per_edge():
    """TEST O7: 'branch' mode still projects each local trait color to the baseline."""
    plot_data = _plot_data(baseline_y=0.0)
    fig = build_figure(
        plot_data, baseline_y=0.0, curtain_color_mode="branch"
    )
    for mesh in _curtain_meshes(fig):
        assert len(mesh.intensity) == len(mesh.y)
        assert len(mesh.intensity) % 2 == 0
        for idx in range(0, len(mesh.intensity), 2):
            assert mesh.intensity[idx + 1] == pytest.approx(mesh.intensity[idx])
            assert mesh.y[idx + 1] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# TEST O8: height mode keeps intensity == y per edge                          #
# --------------------------------------------------------------------------- #
def test_o8_height_mode_intensity_equals_y_per_edge():
    """TEST O8: 'height' mode keeps the historical intensity == Y invariant per trace."""
    plot_data = _plot_data(baseline_y=0.0)
    fig = build_figure(plot_data, baseline_y=0.0, curtain_color_mode="height")
    for mesh in _curtain_meshes(fig):
        assert len(mesh.intensity) == len(mesh.y)
        for y_val, int_val in zip(mesh.y, mesh.intensity):
            assert int_val == pytest.approx(y_val)


# --------------------------------------------------------------------------- #
# Input validation                                                            #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bad", [-0.1, 1.1, float("nan"), float("inf"), float("-inf")])
def test_invalid_mesh_opacity_rejected(bad):
    """0..1 finite only; no silent clamping."""
    plot_data = _plot_data()
    with pytest.raises(ValueError):
        build_figure(plot_data, mesh_opacity=bad)


def test_mesh_opacity_bounds_inclusive():
    """Exact 0.0 and 1.0 are legal boundary values."""
    plot_data = _plot_data()
    for opacity in (0.0, 1.0):
        fig = build_figure(plot_data, mesh_opacity=opacity)
        for mesh in _curtain_meshes(fig):
            assert mesh.opacity == pytest.approx(opacity)


# --------------------------------------------------------------------------- #
# TEST C5 (kept): non-curtain traces unaffected by opacity                    #
# --------------------------------------------------------------------------- #
def test_c5_non_mesh_traces_unaffected_by_opacity():
    """TEST C5: labels, hover traces and the centerline keep their own styling.

    Updated contract (centerline transparency round): the Branch Centerlines
    geometry, colors and text are still opacity-independent, but at mesh
    opacity < 1.0 the trace deliberately carries the same standard opacity so
    foreground curtains can blend over background centerlines (see
    tests/test_centerline_transparency.py TEST CL1-CL10).
    """
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

    for opacity, fig in figs.items():
        assert trace(fig, "Species Labels").to_plotly_json() == baseline_labels
        assert trace(fig, "Terminal Taxa").to_plotly_json() == baseline_tips
        assert trace(fig, "Internal Nodes").to_plotly_json() == baseline_internal

    # Centerline: identical geometry and colors; standard opacity passthrough.
    baseline_line = trace(figs[1.0], "Branch Centerlines").to_plotly_json()
    assert baseline_line.get("opacity") in (None, 1.0)
    for opacity, fig in figs.items():
        line = trace(fig, "Branch Centerlines").to_plotly_json()
        for key in ("x", "y", "z", "line", "mode", "name"):
            assert line[key] == baseline_line[key]
        if opacity < 1.0:
            assert line["opacity"] == pytest.approx(opacity)
        else:
            assert line.get("opacity") in (None, 1.0)

    # Hover markers remain invisible-but-hoverable regardless of mesh opacity.
    for fig in figs.values():
        assert trace(fig, "Terminal Taxa").marker.opacity == pytest.approx(0.0)
        assert trace(fig, "Internal Nodes").marker.opacity == pytest.approx(0.0)
