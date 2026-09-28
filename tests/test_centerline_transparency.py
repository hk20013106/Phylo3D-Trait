"""Feature: centerline compositing with translucent curtains (TEST CL1-CL10).

Root cause fixed in this round:
- The Branch Centerlines Scatter3d trace was always added AFTER the curtain
  Mesh3d traces and was fully opaque (opacity 1.0). With ``--opacity < 1.0``
  the WebGL scene composites it independently of the translucent curtains
  (transparent mesh compositing is draw-order based), so background branch
  centerlines could remain visible as crisp 1-px lines running through
  foreground curtains.
- Fix: when the curtain meshes are rendered with opacity < 1.0, the centerline
  trace is inserted BEFORE all curtain Mesh3d traces and rendered with the same
  standard opacity. The foreground curtains then blend over the background
  centerline (correct alpha compositing) while uncovered top-edge centerlines
  stay visible. At opacity == 1.0 the trace order and payload are exactly the
  legacy ones, so opaque rendering is untouched.

Contract under test:
- CL1: centerline x/y/z (including None separators) never change.
- CL2: centerline colors never change (dark / trait / custom CSS modes).
- CL3: reordering only affects draw order; x/y/z and colors are identical.
- CL4: opacity 1.0 keeps the legacy order and an opaque centerline.
- CL5: curtain mesh payloads are unchanged at opacity 0.9.
- CL6: --no-centerline still renders without any centerline trace.
- CL7: trait-colored centerline still uses the global cmin/cmax/colorscale.
- CL8: hover traces are unchanged between opacity settings.
- CL9: species labels are unchanged; the camera hook may only move mesh3d traces.
- CL10: scientific geometry (curtain aggregate) is unchanged.
"""

from __future__ import annotations

import pytest
import plotly.graph_objects as go

from phylo3d_trait.renderer import (
    TIP_LABEL_CAMERA_ANCHOR_POST_SCRIPT,
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


def _centerline(fig):
    traces = [t for t in fig.data if t.name == "Branch Centerlines"]
    assert len(traces) == 1
    return traces[0]


def _mesh_indices(fig):
    return [i for i, t in enumerate(fig.data) if getattr(t, "type", None) == "mesh3d"]


def _reference_centerline(plot_data):
    """Independent reference construction: branch path from plot_data.segments."""
    edge_map = {}
    for seg in plot_data.segments:
        edge_map.setdefault((seg.parent_id, seg.child_id), []).append(seg)

    xs, ys, zs, colors = [], [], [], []
    for (p_id, c_id), segs in edge_map.items():
        sorted_segs = sorted(segs, key=lambda s: s.segment_index)
        if not sorted_segs:
            continue
        xs.append(sorted_segs[0].x0)
        ys.append(sorted_segs[0].y0)
        zs.append(sorted_segs[0].z0)
        colors.append(sorted_segs[0].trait0)
        for s in sorted_segs:
            xs.append(s.x1)
            ys.append(s.y1)
            zs.append(s.z1)
            colors.append(s.trait1)
        xs.append(None)
        ys.append(None)
        zs.append(None)
        colors.append(sorted_segs[-1].trait1)
    return xs, ys, zs, colors


def _mesh_payloads(fig):
    return [
        t.to_plotly_json()
        for t in fig.data
        if getattr(t, "type", None) == "mesh3d"
    ]


def _aggregate_mesh_geometry(fig):
    xs, ys, zs, ii, jj, kk, inten = [], [], [], [], [], [], []
    offset = 0
    for t in fig.data:
        if getattr(t, "type", None) != "mesh3d":
            continue
        xs.extend(t.x)
        ys.extend(t.y)
        zs.extend(t.z)
        ii.extend(idx + offset for idx in t.i)
        jj.extend(idx + offset for idx in t.j)
        kk.extend(idx + offset for idx in t.k)
        inten.extend(t.intensity)
        offset += len(t.x)
    return xs, ys, zs, ii, jj, kk, inten


# --------------------------------------------------------------------------- #
# TEST CL1: centerline geometry unchanged                                      #
# --------------------------------------------------------------------------- #
def test_cl1_centerline_geometry_unchanged():
    """CL1: x/y/z and separators equal the independent branch-path reference."""
    plot_data = _plot_data()
    for opacity in (1.0, 0.9, 0.7):
        fig = build_figure(plot_data, mesh_opacity=opacity)
        line = _centerline(fig)
        ref_x, ref_y, ref_z, _ = _reference_centerline(plot_data)
        assert list(line.x) == ref_x
        assert list(line.y) == ref_y
        assert list(line.z) == ref_z


# --------------------------------------------------------------------------- #
# TEST CL2: centerline colors unchanged                                        #
# --------------------------------------------------------------------------- #
def test_cl2_centerline_colors_unchanged():
    """CL2: dark mode keeps #2b2b2b; trait mode keeps the per-vertex colors."""
    plot_data = _plot_data()
    line_dark = _centerline(build_figure(plot_data, mesh_opacity=0.9))
    assert line_dark.line.color == "#2b2b2b"

    line_trait = _centerline(
        build_figure(plot_data, mesh_opacity=0.9, centerline_color="trait")
    )
    _, _, _, ref_colors = _reference_centerline(plot_data)
    assert list(line_trait.line.color) == ref_colors

    line_css = _centerline(
        build_figure(plot_data, mesh_opacity=0.9, centerline_color="#123456")
    )
    assert line_css.line.color == "#123456"


# --------------------------------------------------------------------------- #
# TEST CL3: reordering never changes coordinates or colors                     #
# --------------------------------------------------------------------------- #
def test_cl3_reorder_changes_only_draw_order():
    """CL3: opacity 0.9 puts the line before curtains; 1.0 keeps legacy order."""
    plot_data = _plot_data()
    fig_1 = build_figure(plot_data, mesh_opacity=1.0)
    fig_9 = build_figure(plot_data, mesh_opacity=0.9)

    line_1 = _centerline(fig_1)
    line_9 = _centerline(fig_9)
    assert list(line_1.x) == list(line_9.x)
    assert list(line_1.y) == list(line_9.y)
    assert list(line_1.z) == list(line_9.z)
    assert list(line_1.line.color) == list(line_9.line.color) if isinstance(line_1.line.color, (list, tuple)) else line_1.line.color == line_9.line.color

    # Draw-order placement
    mesh_idx_1 = _mesh_indices(fig_1)
    idx_1 = [i for i, t in enumerate(fig_1.data) if t.name == "Branch Centerlines"][0]
    assert idx_1 > max(mesh_idx_1), "opacity 1.0 keeps the legacy order (line after meshes)"

    mesh_idx_9 = _mesh_indices(fig_9)
    idx_9 = [i for i, t in enumerate(fig_9.data) if t.name == "Branch Centerlines"][0]
    assert idx_9 < min(mesh_idx_9), "opacity < 1.0 draws the line before the curtains"


# --------------------------------------------------------------------------- #
# TEST CL4: opacity 1.0 behavior unchanged                                     #
# --------------------------------------------------------------------------- #
def test_cl4_opacity_1_keeps_legacy_centerline_styling():
    """CL4: at 1.0 the centerline stays opaque and no trace is reordered."""
    plot_data = _plot_data()
    fig = build_figure(plot_data, mesh_opacity=1.0)
    line = _centerline(fig)
    assert line.opacity in (None, 1.0)
    assert line.visible in (None, True)


# --------------------------------------------------------------------------- #
# TEST CL5: opacity 0.9 mesh properties unchanged                              #
# --------------------------------------------------------------------------- #
def test_cl5_mesh_payloads_unchanged_by_centerline_reorder():
    """CL5: curtain payloads at 0.9 identical with and without centerline."""
    plot_data = _plot_data()
    fig_with = build_figure(plot_data, mesh_opacity=0.9, show_centerline=True)
    fig_without = build_figure(plot_data, mesh_opacity=0.9, show_centerline=False)
    assert _mesh_payloads(fig_with) == _mesh_payloads(fig_without)


# --------------------------------------------------------------------------- #
# TEST CL6: --no-centerline still works                                        #
# --------------------------------------------------------------------------- #
def test_cl6_no_centerline_renders_without_trace():
    """CL6: show_centerline=False removes the trace entirely (any opacity)."""
    plot_data = _plot_data()
    for opacity in (1.0, 0.9, 0.7):
        fig = build_figure(plot_data, mesh_opacity=opacity, show_centerline=False)
        assert [t for t in fig.data if t.name == "Branch Centerlines"] == []


# --------------------------------------------------------------------------- #
# TEST CL7: trait-colored centerline color domain                              #
# --------------------------------------------------------------------------- #
def test_cl7_trait_centerline_uses_global_color_domain():
    """CL7: trait color mode keeps the global cmin/cmax/colorscale/reversescale."""
    plot_data = _plot_data()
    fig = build_figure(
        plot_data,
        mesh_opacity=0.7,
        centerline_color="trait",
        reverse_colorscale=True,
    )
    line = _centerline(fig)
    mesh = [t for t in fig.data if getattr(t, "type", None) == "mesh3d"][0]
    assert line.line.cmin == pytest.approx(plot_data.trait_min)
    assert line.line.cmax == pytest.approx(plot_data.trait_max)
    # Plotly resolves the colorscale name to its full mapping; the centerline
    # must share the exact same global mapping as the curtain meshes.
    assert line.line.colorscale == mesh.colorscale
    assert line.line.reversescale is True


# --------------------------------------------------------------------------- #
# TEST CL8: hover traces unchanged                                             #
# --------------------------------------------------------------------------- #
def test_cl8_hover_traces_unchanged_across_opacities():
    """CL8: tip/internal hover payloads are identical for opacity 1.0 / 0.7."""
    plot_data = _plot_data()
    fig_1 = build_figure(plot_data, mesh_opacity=1.0)
    fig_7 = build_figure(plot_data, mesh_opacity=0.7)

    def by_name(fig, name):
        return [t for t in fig.data if t.name == name][0].to_plotly_json()

    for name in ("Internal Nodes", "Terminal Taxa"):
        assert by_name(fig_1, name) == by_name(fig_7, name)


# --------------------------------------------------------------------------- #
# TEST CL9: species labels unchanged; hook only moves meshes                   #
# --------------------------------------------------------------------------- #
def test_cl9_species_labels_unchanged_and_hook_mesh_only():
    """CL9: label payloads identical across opacities; hook moves mesh3d only."""
    plot_data = _plot_data()
    fig_1 = build_figure(plot_data, mesh_opacity=1.0)
    fig_9 = build_figure(plot_data, mesh_opacity=0.9)

    def labels(fig):
        return [t for t in fig.data if t.name == "Species Labels"][0].to_plotly_json()

    assert labels(fig_1) == labels(fig_9)

    script = TIP_LABEL_CAMERA_ANCHOR_POST_SCRIPT
    assert "if (gd.data[i].type === 'mesh3d') { idxs.push(i); }" in script


# --------------------------------------------------------------------------- #
# TEST CL10: scientific curtain geometry unchanged                             #
# --------------------------------------------------------------------------- #
def test_cl10_scientific_geometry_unchanged():
    """CL10: per-edge meshes still recombine into the legacy aggregate."""
    plot_data = _plot_data()
    baseline_y = plot_data.baseline_y
    fig = build_figure(plot_data, mesh_opacity=0.9)
    legacy = _build_branch_curtains_geometry(plot_data, baseline_y)
    combined = _aggregate_mesh_geometry(fig)
    for got, want in zip(combined, legacy):
        assert list(got) == list(want)
