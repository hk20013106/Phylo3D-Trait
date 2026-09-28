"""Plotly 3D renderer for Phylo3D-Trait visualization.

Builds interactive 3D WebGL plots mapping:
- X axis: Tree layout
- Y axis: Trait value (Height)
- Z axis: Evolutionary time before present

Visual representations:
- Continuous vertical curtain / ribbon surfaces (Mesh3d) descending from each
  branch's trait height down to a common trait baseline plane. Each biological
  parent->child edge gets its own Mesh3d trace so that, with opacity < 1.0,
  WebGL can depth-order the overlapping transparent ribbons safely (all traces
  share one global color domain and exactly one colorbar).
- Crisp branch top outlines (Scatter3d lines).
- Fixed text labels for terminal taxa anchored just outside the present time
  (Z = time_min) baseline plane, extending outward from each tip. World
  coordinates are fixed; the interactive HTML flips only the text anchor
  (middle right <-> middle left) when the camera orbits to the opposite side.
- Pure white / transparent background with clean axis gridlines.
- eLife-style camera preset with screen-vertical Y (Trait) and +Z foreground (MRCA).
- Global trait normalization with optional independent color reversal and
  branch-projected curtain coloring.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import math
import plotly.graph_objects as go

from phylo3d_trait.models import EdgeSegment, PlotData
from phylo3d_trait.tree import annotate_tree, parse_tree


CAMERA_PRESETS: Dict[str, Dict[str, Any]] = {
    "elife": dict(
        eye=dict(x=1.35, y=0.65, z=2.1),
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        projection=dict(type="orthographic"),
    ),
    "root_front": dict(
        eye=dict(x=0.0, y=0.5, z=2.4),
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        projection=dict(type="orthographic"),
    ),
    "tips_front": dict(
        eye=dict(x=1.35, y=0.65, z=-2.1),
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        projection=dict(type="orthographic"),
    ),
}

# Generic default outward offset for terminal species labels, expressed as a
# fraction of the Time-before-present span. Because Plotly normalizes scene
# aspects per axis, a fraction of the span is scale-free: it adapts to any
# phylogeny without hard-coded absolute distances.
DEFAULT_TIP_LABEL_OFFSET_FRACTION = 0.03

# Minimal camera-aware hook embedded into the generated interactive HTML.
#
# Species labels keep their exact world coordinates (x, y, z) at all times.
# Only the Scatter3d text anchor may switch between 'middle right' and
# 'middle left', because textposition is a screen-space (billboard) property:
# at the default eLife camera the tree-outward direction (beyond the present
# plane, world -Z) projects to screen-right, so text must hang right of its
# anchor; after orbiting to the opposite side it projects to screen-left, so
# the same text must hang left, otherwise it would extend into the tree.
#
# The sign of the outward direction's screen-x component is the sign of
# (eye.x - center.x) for a lookAt camera with world up (0, 1, 0) (see
# tests/test_tip_label_outward.py::test_l3_l4_* for the derivation).
# The hook listens to plotly.js gl3d camera events only and performs:
#   1. a single textposition restyle of the Species Labels trace;
#   2. a draw-order sort of the curtain Mesh3d traces (painter's algorithm,
#      camera-nearest last) so transparent curtains composite correctly.
# It never touches coordinates, text, fonts, scientific values or tip order.
#
# Camera-source notes (verified in a real browser):
# - On 'plotly_relayout' (commit: mouse release or programmatic relayout) the
#   layout camera is already updated, while the internal glplot camera may lag
#   one frame; prefer the layout camera there.
# - During a live drag, 'plotly_relayouting' fires before the layout camera is
#   saved, so prefer the live glplot camera there.
# - A one-frame requestAnimationFrame re-check converges regardless of event
#   ordering and lets the curtain sort use the settled camera.
#
# Draw-order note (verified with pixel measurements): splitting the curtains
# into one Mesh3d per edge is not enough on its own because plotly draws
# transparent mesh traces strictly in trace order. Painter's back-to-front
# order must be restored per camera: objects whose NEAREST point is farther
# from the camera are drawn first, so the camera-nearest curtains are painted
# last and dominate the alpha blend. This removed the "outer curtain almost
# fully transparent" artifact in both tested opposite camera orientations.
TIP_LABEL_CAMERA_ANCHOR_POST_SCRIPT = """/* Phylo3D tip-label camera anchor + transparent curtain draw order:
   species labels keep fixed world coordinates (only the text anchor flips),
   and curtain Mesh3d traces are repainted back-to-front for the current
   camera so transparency compositing stays correct. */
(function () {
    function comp(vec, key) {
        if (!vec) { return 0; }
        if (typeof vec[key] === 'number') { return vec[key]; }
        var idx = (key === 'x') ? 0 : (key === 'y') ? 1 : 2;
        return (typeof vec[idx] === 'number') ? vec[idx] : 0;
    }
    function tipLabelSide(eyeX) {
        return eyeX < 0 ? 'middle left' : 'middle right';
    }
    function readCamera(gd, preferLive) {
        var scene = gd._fullLayout && gd._fullLayout.scene;
        var state = scene && scene._scene;
        var live = state && state.glplot && state.glplot.camera;
        var layoutCam = scene && scene.camera;
        var primary = preferLive ? (live || layoutCam) : (layoutCam || live);
        var fallback = preferLive ? layoutCam : live;
        var eye = (primary && primary.eye) || (fallback && fallback.eye);
        var center = (primary && primary.center) || (fallback && fallback.center);
        return { eye: eye, center: center };
    }
    function applyTipLabelSide(gd, preferLive) {
        var idx = -1;
        for (var i = 0; i < gd.data.length; i++) {
            if (gd.data[i].name === 'Species Labels') { idx = i; break; }
        }
        if (idx < 0) { return; }
        var cam = readCamera(gd, preferLive);
        var eyeX = comp(cam.eye, 'x') - comp(cam.center, 'x');
        var side = tipLabelSide(eyeX);
        if (gd.data[idx].textposition !== side) {
            Plotly.restyle(gd, {textposition: side}, [idx]);
        }
    }
    function curtainTraceIndices(gd) {
        var idxs = [];
        for (var i = 0; i < gd.data.length; i++) {
            if (gd.data[i].type === 'mesh3d') { idxs.push(i); }
        }
        return idxs;
    }
    function sortCurtainTraces(gd) {
        if (typeof Plotly.moveTraces !== 'function') { return; }
        var scene = gd._fullLayout && gd._fullLayout.scene;
        var gl = scene && scene._scene && scene._scene.glplot;
        if (!gl || !gl.camera || !gl.camera.matrix || !gl.objects) { return; }
        var idxs = curtainTraceIndices(gd);
        if (idxs.length < 2 || gl.objects.length !== gd.data.length) { return; }
        var m = gl.camera.matrix;
        var entries = [];
        for (var k = 0; k < idxs.length; k++) {
            var obj = gl.objects[idxs[k]];
            var b = obj && obj.bounds;
            if (!b || b.length !== 2) { return; }
            // Depth (view-space z) of the object's NEAREST bounding-box corner.
            var nearest = -Infinity;
            for (var xi = 0; xi < 2; xi++) {
                for (var yi = 0; yi < 2; yi++) {
                    for (var zi = 0; zi < 2; zi++) {
                        var depth = m[2] * b[xi][0] + m[6] * b[yi][1] + m[10] * b[zi][2] + m[14];
                        if (depth > nearest) { nearest = depth; }
                    }
                }
            }
            entries.push({ trace: idxs[k], nearest: nearest });
        }
        entries.sort(function (a, b) { return a.nearest - b.nearest; });
        var order = [];
        for (var e = 0; e < entries.length; e++) { order.push(entries[e].trace); }
        var same = true;
        for (var c = 0; c < order.length; c++) {
            if (order[c] !== idxs[c]) { same = false; break; }
        }
        if (same) { return; }
        // Build destination array: for each trace in idxs (current order),
        // find its absolute target position in gd.data.
        // Curtains occupy contiguous slots starting at idxs[0].
        var dest = [];
        for (var d = 0; d < idxs.length; d++) {
            var rank = order.indexOf(idxs[d]);
            dest.push(idxs[rank]);
        }
        Plotly.moveTraces(gd, idxs, dest);
    }
    function applyCameraResponse(gd, preferLive, doSort) {
        applyTipLabelSide(gd, preferLive);
        if (doSort) { sortCurtainTraces(gd); }
    }
    function scheduleCameraResponse(gd, preferLive, doSort) {
        applyCameraResponse(gd, preferLive, doSort);
        if (typeof window.requestAnimationFrame === 'function') {
            window.requestAnimationFrame(function () {
                applyCameraResponse(gd, preferLive, doSort);
            });
        }
    }
    var gd = document.getElementById('{plot_id}');
    if (!gd) { return; }
    gd.on('plotly_relayout', function () { scheduleCameraResponse(gd, false, true); });
    gd.on('plotly_relayouting', function () { scheduleCameraResponse(gd, true, false); });
    scheduleCameraResponse(gd, false, true);
})();
"""


def build_plot_data(
    tree_input: Any,
    trait_values: Dict[str, float],
    num_segments: int = 10,
    colorscale: str = "Turbo",
    title: str = "3D Phylogenetic Tree with Continuous Trait Evolution",
    baseline_y: Optional[float] = None,
    trait_display_range: Optional[Tuple[float, float]] = None,
    baseline_raw_value: Optional[float] = None,
    trait_display_offset: Optional[float] = None,
) -> PlotData:
    """Convenience helper to parse tree, annotate traits, and construct PlotData.

    Args:
        tree_input: Newick/Nexus string, path, or Tree object.
        trait_values: Dictionary of node IDs to trait values.
        num_segments: Number of interpolation segments per branch.
        colorscale: Plotly colorscale name.
        title: Visualization title.
        baseline_y: Optional custom baseline Trait height for curtain meshes.
        trait_display_range: Optional custom (start, end) target display range for linear remapping.
        baseline_raw_value: Optional custom numeric trait value to label at baseline Y.

    Returns:
        PlotData object.
    """
    from Bio.Phylo.BaseTree import Clade, Tree

    if isinstance(tree_input, (Tree, Clade)) or hasattr(tree_input, "get_terminals"):
        tree = tree_input
    else:
        tree = parse_tree(tree_input)

    data = annotate_tree(
        tree_input=tree,
        trait_values=trait_values,
        num_segments=num_segments,
        colorscale=colorscale,
        title=title,
        trait_display_range=trait_display_range,
        trait_display_offset=trait_display_offset,
    )
    if trait_display_offset is not None:
        data.trait_display_offset = trait_display_offset
    if baseline_y is not None:
        data.baseline_y = baseline_y
    if baseline_raw_value is not None:
        data.baseline_raw_value = baseline_raw_value
    return data


def _build_single_edge_curtain_geometry(
    branch_pts: List[Tuple[float, float, float]],
    baseline_y: float,
    curtain_color_mode: str = "height",
) -> Tuple[List[float], List[float], List[float], List[int], List[int], List[int], List[float]]:
    """Triangulate ONE biological parent->child edge into a curtain ribbon.

    Vertices are interleaved as Top_k, Bottom_k along the sampled branch
    points. Triangle indices are LOCAL to the returned arrays (0-based), so
    the result can be used directly as a standalone Mesh3d trace.

    For each adjacent step (k, k+1) two triangles are generated:
        Triangle A: (Top_k, Bottom_k, Top_{k+1})
        Triangle B: (Bottom_k, Bottom_{k+1}, Top_{k+1})

    Args:
        branch_pts: Sampled (x, y, z) points along this edge, in order.
        baseline_y: The constant Y height of the baseline plane.
        curtain_color_mode: 'height' keeps the historical vertical gradient
            (vertex color intensity equals vertex Y); 'branch' projects the
            local branch trait color vertically to the baseline.

    Returns:
        Tuple of (x, y, z, i, j, k, intensity) for this single edge.
    """
    x: List[float] = []
    y: List[float] = []
    z: List[float] = []
    i: List[int] = []
    j: List[int] = []
    k: List[int] = []
    intensity: List[float] = []

    for xk, yk, zk in branch_pts:
        # Top vertex is always colored by the local branch trait.
        x.append(xk)
        y.append(yk)
        z.append(zk)
        intensity.append(yk)

        # Bottom vertex sits on the baseline plane; historical 'height' mode
        # colors it by geometric height, 'branch' mode extrudes the local
        # top-branch color straight down.
        x.append(xk)
        y.append(baseline_y)
        z.append(zk)
        intensity.append(baseline_y if curtain_color_mode == "height" else yk)

    for step in range(len(branch_pts) - 1):
        top_k = 2 * step
        bot_k = 2 * step + 1
        top_k1 = 2 * (step + 1)
        bot_k1 = 2 * (step + 1) + 1

        # Triangle A: (Top_k, Bottom_k, Top_{k+1})
        i.append(top_k)
        j.append(bot_k)
        k.append(top_k1)

        # Triangle B: (Bottom_k, Bottom_{k+1}, Top_{k+1})
        i.append(bot_k)
        j.append(bot_k1)
        k.append(top_k1)

    return x, y, z, i, j, k, intensity


def _build_branch_curtain_geometries(
    plot_data: PlotData,
    baseline_y: float,
    curtain_color_mode: str = "height",
) -> List[Dict[str, Any]]:
    """Construct one curtain geometry per biological parent->child edge.

    Each element of the returned list describes a single biological edge and
    can be rendered as its own Mesh3d trace (see ``build_figure``), which lets
    WebGL depth-order overlapping transparent curtains instead of blending
    thousands of triangles from every edge inside one draw call.

    Args:
        plot_data: PlotData containing edge segments and scaling bounds.
        baseline_y: The constant Y height of the baseline plane.
        curtain_color_mode: 'height' or 'branch' (see
            :func:`_build_single_edge_curtain_geometry`).

    Returns:
        List of dicts with keys: parent_id, child_id, x, y, z, i, j, k,
        intensity. Triangle indices are local to each edge's vertex arrays.
    """
    if curtain_color_mode not in {"height", "branch"}:
        raise ValueError(
            "curtain_color_mode must be 'height' or 'branch', "
            f"got {curtain_color_mode!r}"
        )

    # Group segments by parent-child edge, preserving first-seen edge order.
    edge_map: Dict[tuple, List[EdgeSegment]] = {}
    for seg in plot_data.segments:
        key = (seg.parent_id, seg.child_id)
        if key not in edge_map:
            edge_map[key] = []
        edge_map[key].append(seg)

    geometries: List[Dict[str, Any]] = []
    for (p_id, c_id), segs in edge_map.items():
        sorted_segs = sorted(segs, key=lambda s: s.segment_index)
        if not sorted_segs:
            continue

        branch_pts = [(sorted_segs[0].x0, sorted_segs[0].y0, sorted_segs[0].z0)]
        for s in sorted_segs:
            branch_pts.append((s.x1, s.y1, s.z1))

        x, y, z, i, j, k, intensity = _build_single_edge_curtain_geometry(
            branch_pts, baseline_y, curtain_color_mode
        )
        geometries.append(
            {
                "parent_id": p_id,
                "child_id": c_id,
                "x": x,
                "y": y,
                "z": z,
                "i": i,
                "j": j,
                "k": k,
                "intensity": intensity,
            }
        )

    return geometries


def _build_branch_curtains_geometry(
    plot_data: PlotData,
    baseline_y: float,
    curtain_color_mode: str = "height",
) -> Tuple[List[float], List[float], List[float], List[int], List[int], List[int], List[float]]:
    """Aggregate every per-edge curtain into one combined geometry.

    This is the legacy single-trace layout, kept as the reference
    implementation used by tests and validation to prove that per-edge
    Mesh3d traces recombine exactly into the historical aggregate vertices,
    triangles and intensities.

    Returns:
        Tuple of (mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity).
    """
    mesh_x: List[float] = []
    mesh_y: List[float] = []
    mesh_z: List[float] = []
    mesh_i: List[int] = []
    mesh_j: List[int] = []
    mesh_k: List[int] = []
    mesh_intensity: List[float] = []

    vertex_offset = 0
    for geom in _build_branch_curtain_geometries(
        plot_data, baseline_y, curtain_color_mode
    ):
        mesh_x.extend(geom["x"])
        mesh_y.extend(geom["y"])
        mesh_z.extend(geom["z"])
        mesh_intensity.extend(geom["intensity"])
        mesh_i.extend(idx + vertex_offset for idx in geom["i"])
        mesh_j.extend(idx + vertex_offset for idx in geom["j"])
        mesh_k.extend(idx + vertex_offset for idx in geom["k"])
        vertex_offset += len(geom["x"])

    return mesh_x, mesh_y, mesh_z, mesh_i, mesh_j, mesh_k, mesh_intensity


def _generate_rescaled_ticks(
    plot_data: PlotData,
    baseline_y: Optional[float] = None,
    baseline_raw_value: Optional[float] = None,
    num_ticks: int = 5,
    include_baseline: bool = True,
) -> Tuple[List[float], List[str]]:
    """Generate (tickvals, ticktext) in display space labeled with raw scientific trait values.

    Args:
        plot_data: PlotData container with raw trait bounds and trait_display_range.
        baseline_y: Custom baseline Y plane height.
        baseline_raw_value: Optional override for the numeric scientific value displayed at baseline Y.
        num_ticks: Number of trait tick steps across the display domain (default: 5).
        include_baseline: Whether to include an explicit baseline tick when the
            baseline lies below the displayed trait domain.

    Returns:
        Tuple of (tickvals, ticktext) for Plotly axis and colorbar.
    """
    if plot_data.trait_display_offset is not None:
        eff_baseline = baseline_y if baseline_y is not None else plot_data.baseline_y
        if eff_baseline is None:
            eff_baseline = 0.0

        d_min = plot_data.trait_min
        d_max = plot_data.trait_max

        if include_baseline and eff_baseline is not None and eff_baseline < d_min:
            raw_start = round(plot_data.display_to_raw(eff_baseline))
        else:
            raw_start = round(plot_data.display_to_raw(d_min))
        raw_end = math.ceil(plot_data.display_to_raw(d_max))

        span = raw_end - raw_start
        if span <= 8:
            step = 1.0
        elif span <= 16:
            step = 2.0
        else:
            step = max(1.0, round(span / float(num_ticks - 1)))

        raw_vals = []
        curr = float(raw_start)
        while curr <= raw_end + 1e-5:
            raw_vals.append(curr)
            curr += step

        tickvals = [round(plot_data.raw_to_display(r), 6) for r in raw_vals]
        ticktext = [str(int(r)) if abs(r - round(r)) < 1e-6 else f"{r:.4f}".rstrip("0").rstrip(".") for r in raw_vals]
        return tickvals, ticktext

    if plot_data.trait_display_range is None:
        return [], []

    d_start, d_end = plot_data.trait_display_range
    d_min = min(d_start, d_end)
    d_max = max(d_start, d_end)

    if d_max == d_min:
        disp_ticks = [d_min]
    else:
        step = (d_max - d_min) / float(num_ticks - 1)
        disp_ticks = [d_min + i * step for i in range(num_ticks)]

    tickvals: List[float] = []
    ticktext: List[str] = []

    eff_baseline = baseline_y if baseline_y is not None else plot_data.baseline_y
    eff_baseline_raw = (
        baseline_raw_value
        if baseline_raw_value is not None
        else plot_data.baseline_raw_value
    )

    # If baseline is explicitly below the display trait domain, add bottom tick with raw numeric baseline value
    if include_baseline and eff_baseline is not None and eff_baseline < d_min - 1e-4:
        if eff_baseline_raw is not None:
            raw_b = eff_baseline_raw
        elif d_start > d_end:  # reverse transform: lower display Y = higher raw trait
            raw_b = plot_data.raw_trait_max + 2.0
        else:  # forward transform: lower display Y = lower raw trait
            raw_b = plot_data.raw_trait_min - 2.0

        tickvals.append(float(eff_baseline))
        if abs(raw_b - round(raw_b)) < 1e-6:
            b_label = f"{int(round(raw_b))}"
        else:
            b_label = f"{raw_b:.4f}".rstrip("0").rstrip(".")
        ticktext.append(b_label)

    for d_val in disp_ticks:
        raw_val = plot_data.display_to_raw(d_val)
        tickvals.append(round(d_val, 6))
        if abs(raw_val - round(raw_val)) < 1e-6:
            label = f"{int(round(raw_val))}"
        else:
            label = f"{raw_val:.4f}".rstrip("0").rstrip(".")
        ticktext.append(label)

    return tickvals, ticktext


def build_figure(
    plot_data: PlotData,
    title: Optional[str] = None,
    branch_width: float = 1.0,
    show_tip_labels: bool = True,
    aspect_ratio: Optional[Dict[str, float]] = None,
    show_mesh: bool = True,
    mesh_opacity: float = 1.0,
    show_centerline: bool = True,
    centerline_color: str = "dark",
    baseline_y: Optional[float] = None,
    baseline_raw_value: Optional[float] = None,
    show_node_markers: bool = False,
    internal_marker_size: float = 4.0,
    background: str = "white",
    camera_preset: str = "elife",
    custom_camera: Optional[Dict[str, Any]] = None,
    reverse_colorscale: bool = False,
    curtain_color_mode: str = "height",
    trait_axis_scale: float = 1.0,
    tip_label_offset: Optional[float] = None,
) -> go.Figure:
    """Construct an interactive Plotly 3D Figure in clean publication style with eLife-style camera.

    Args:
        plot_data: PlotData containing annotated nodes, edge segments, and scaling limits.
        title: Optional title override.
        branch_width: Line width for 3D branch top outline (default: 1.0).
        show_tip_labels: Whether to display text labels for tip taxa on the present baseline plane, extending outward from each tip (default: True).
        aspect_ratio: Optional custom aspect ratio dictionary {'x': float, 'y': float, 'z': float}.
        show_mesh: Whether to render continuous vertical curtain meshes.
        mesh_opacity: Standard opacity for curtain meshes (0.0 to 1.0, default 1.0).
            Applied identically to every per-edge curtain trace. Must be finite
            and within 0.0..1.0 (rejected loudly otherwise).
        show_centerline: Whether to render top-edge outline along branches.
        centerline_color: Color mode for centerline ('dark', 'trait', or CSS color).
        baseline_y: Custom baseline Y plane height (defaults to plot_data.baseline_y).
        baseline_raw_value: Custom numeric scientific trait value to label at baseline Y.
        show_node_markers: Whether to render ancestral node markers (default: False).
        internal_marker_size: Marker size if show_node_markers is True.
        background: 'white' (default) or 'transparent'.
        camera_preset: Preset camera angle ('elife', 'root_front', 'tips_front', default: 'elife').
        custom_camera: Optional dict to override camera config completely.
        reverse_colorscale: Reverse only the color mapping while leaving trait
            heights and scientific values unchanged.
        curtain_color_mode: 'height' (default) colors curtains by geometric
            Y height; 'branch' projects each local branch trait color vertically
            to the baseline.
        trait_axis_scale: Visual aspect scale factor for the Trait (Y) dimension
            only (default: 1.0). Purely visual: it multiplies the scene Trait
            aspect ratio and changes nothing else (no trait values, ticks,
            hover, color domain, time or tree-layout geometry). Example: 0.5
            halves the Trait visual height. Must be finite and > 0.
        tip_label_offset: Outward offset of terminal species labels beyond the
            present plane (``time_min``), expressed as a fraction of the
            Time-before-present span. ``None`` (default) uses
            ``DEFAULT_TIP_LABEL_OFFSET_FRACTION`` (0.03). 0.0 keeps labels
            exactly on the present plane. Must be finite and >= 0.

    Returns:
        Plotly go.Figure configured for interactive 3D display.
    """
    plot_title = title if title is not None else plot_data.title

    try:
        trait_axis_scale = float(trait_axis_scale)
    except (TypeError, ValueError) as err:
        raise ValueError(
            f"trait_axis_scale must be a finite positive number, got {trait_axis_scale!r}"
        ) from err
    if not math.isfinite(trait_axis_scale) or trait_axis_scale <= 0:
        raise ValueError(
            f"trait_axis_scale must be a finite positive number, got {trait_axis_scale!r}"
        )

    if tip_label_offset is not None:
        try:
            tip_label_offset = float(tip_label_offset)
        except (TypeError, ValueError) as err:
            raise ValueError(
                "tip_label_offset must be a finite non-negative fraction of the "
                f"Time-before-present span, got {tip_label_offset!r}"
            ) from err
        if not math.isfinite(tip_label_offset) or tip_label_offset < 0:
            raise ValueError(
                "tip_label_offset must be a finite non-negative fraction of the "
                f"Time-before-present span, got {tip_label_offset!r}"
            )

    try:
        mesh_opacity = float(mesh_opacity)
    except (TypeError, ValueError) as err:
        raise ValueError(
            f"mesh_opacity must be a finite number between 0.0 and 1.0, got {mesh_opacity!r}"
        ) from err
    if not math.isfinite(mesh_opacity) or mesh_opacity < 0.0 or mesh_opacity > 1.0:
        raise ValueError(
            f"mesh_opacity must be a finite number between 0.0 and 1.0, got {mesh_opacity!r}"
        )

    eff_baseline_y = baseline_y if baseline_y is not None else plot_data.baseline_y
    if eff_baseline_y is None:
        eff_baseline_y = plot_data.trait_min

    if curtain_color_mode not in {"height", "branch"}:
        raise ValueError(
            "curtain_color_mode must be 'height' or 'branch', "
            f"got {curtain_color_mode!r}"
        )

    # Branch/marker colors are always normalized to the actual displayed trait
    # domain. In historical height mode the curtain additionally colors its
    # geometric baseline, so the mesh color domain must include baseline_y.
    trait_cmin = plot_data.trait_min
    trait_cmax = plot_data.trait_max
    if trait_cmin == trait_cmax:
        trait_cmin -= 0.5
        trait_cmax += 0.5

    if curtain_color_mode == "height":
        mesh_cmin = min(plot_data.trait_min, eff_baseline_y)
        mesh_cmax = max(plot_data.trait_max, eff_baseline_y)
        if mesh_cmin == mesh_cmax:
            mesh_cmin -= 0.5
            mesh_cmax += 0.5
    else:
        mesh_cmin = trait_cmin
        mesh_cmax = trait_cmax

    fig = go.Figure()

    # 0. Branch top centerline outlines (Scatter3d lines).
    #    The line data (x/y/z/colors) is always identical; only the DRAW ORDER
    #    and the trace opacity may change with the curtain opacity:
    #    - opacity == 1.0: the trace is appended after the curtain meshes,
    #      fully opaque (legacy behavior; verified pixel-identical).
    #    - opacity < 1.0: the trace is inserted BEFORE the curtain meshes and
    #      carries the same standard opacity (90% opaque for 0.9). WebGL
    #      composites transparent objects in draw order, so foreground curtains
    #      then blend over background centerlines instead of leaving them
    #      visible as crisp 1-px lines running through translucent curtains;
    #      uncovered top-edge centerlines stay visible.
    #    The camera hook never moves this trace (it only reorders mesh3d
    #    traces), so the line stays in its assigned slot after any camera
    #    change.
    branch_x: List[Optional[float]] = []
    branch_y: List[Optional[float]] = []
    branch_z: List[Optional[float]] = []
    branch_colors: List[Optional[float]] = []

    edge_map: Dict[tuple, List[EdgeSegment]] = {}
    for seg in plot_data.segments:
        key = (seg.parent_id, seg.child_id)
        if key not in edge_map:
            edge_map[key] = []
        edge_map[key].append(seg)

    for (p_id, c_id), segs in edge_map.items():
        sorted_segs = sorted(segs, key=lambda s: s.segment_index)
        if not sorted_segs:
            continue

        # Start vertex
        branch_x.append(sorted_segs[0].x0)
        branch_y.append(sorted_segs[0].y0)
        branch_z.append(sorted_segs[0].z0)
        branch_colors.append(sorted_segs[0].trait0)

        # End vertices
        for s in sorted_segs:
            branch_x.append(s.x1)
            branch_y.append(s.y1)
            branch_z.append(s.z1)
            branch_colors.append(s.trait1)

        # Disconnect line from next branch
        branch_x.append(None)
        branch_y.append(None)
        branch_z.append(None)
        branch_colors.append(sorted_segs[-1].trait1)

    centerline_before_meshes = (
        bool(show_mesh and plot_data.segments) and mesh_opacity < 1.0
    )

    def _add_centerline_trace() -> None:
        if centerline_color == "trait":
            line_cfg = dict(
                color=branch_colors,
                colorscale=plot_data.colorscale,
                cmin=trait_cmin,
                cmax=trait_cmax,
                reversescale=reverse_colorscale,
                width=branch_width,
            )
        elif centerline_color == "dark":
            line_cfg = dict(
                color="#2b2b2b",
                width=branch_width,
            )
        else:
            line_cfg = dict(
                color=centerline_color,
                width=branch_width,
            )

        trace_kwargs: Dict[str, Any] = dict(
            x=branch_x,
            y=branch_y,  # Y is TRAIT
            z=branch_z,  # Z is TIME
            mode="lines",
            line=line_cfg,
            hoverinfo="skip",
            name="Branch Centerlines",
            showlegend=False,
        )
        if centerline_before_meshes:
            trace_kwargs["opacity"] = mesh_opacity
        fig.add_trace(go.Scatter3d(**trace_kwargs))

    if show_centerline and branch_x and centerline_before_meshes:
        _add_centerline_trace()

    # 1. Build continuous curtain mesh surfaces (Mesh3d)
    is_transformed = (plot_data.trait_display_range is not None) or (plot_data.trait_display_offset is not None)
    colorbar_title = "Trait Value"
    y_axis_title = "Trait value"
    y_tickvals, y_ticktext = _generate_rescaled_ticks(
        plot_data=plot_data,
        baseline_y=eff_baseline_y,
        baseline_raw_value=baseline_raw_value,
        num_ticks=5,
        include_baseline=True,
    )
    colorbar_tickvals, colorbar_ticktext = _generate_rescaled_ticks(
        plot_data=plot_data,
        baseline_y=eff_baseline_y,
        baseline_raw_value=baseline_raw_value,
        num_ticks=5,
        include_baseline=(curtain_color_mode == "height"),
    )

    if show_mesh and plot_data.segments:
        # One Mesh3d trace per biological parent->child edge. This preserves
        # the exact aggregate geometry (see _build_branch_curtains_geometry)
        # while letting WebGL depth-order overlapping transparent ribbons;
        # a single giant trace blends thousands of interleaved triangles in
        # buffer order, which makes opacity=0.9 look far more transparent than
        # requested. All curtain traces still share ONE scientific color
        # domain, ONE opacity and exactly ONE colorbar.
        curtain_geometries = _build_branch_curtain_geometries(
            plot_data=plot_data,
            baseline_y=eff_baseline_y,
            curtain_color_mode=curtain_color_mode,
        )

        if curtain_geometries:
            cb_dict = dict(
                title=dict(text=colorbar_title, side="top", font=dict(size=12, color="#333333")),
                thickness=18,
                len=0.75,
                x=1.02,
            )
            if is_transformed and colorbar_tickvals:
                cb_dict["tickmode"] = "array"
                cb_dict["tickvals"] = colorbar_tickvals
                cb_dict["ticktext"] = colorbar_ticktext

            colorbar_attached = False
            for geom in curtain_geometries:
                if not geom["x"] or not geom["i"]:
                    continue
                fig.add_trace(
                    go.Mesh3d(
                        x=geom["x"],
                        y=geom["y"],  # Y is TRAIT (Top) and baseline_y (Bottom)
                        z=geom["z"],  # Z is TIME
                        i=geom["i"],
                        j=geom["j"],
                        k=geom["k"],
                        intensity=geom["intensity"],
                        colorscale=plot_data.colorscale,
                        cmin=mesh_cmin,
                        cmax=mesh_cmax,
                        reversescale=reverse_colorscale,
                        opacity=mesh_opacity,
                        flatshading=False,
                        lighting=dict(
                            ambient=0.85,
                            diffuse=0.5,
                            specular=0.08,
                            roughness=0.8,
                        ),
                        hoverinfo="skip",
                        name=f"Branch Curtain: {geom['parent_id']} -> {geom['child_id']}",
                        showlegend=False,
                        showscale=not colorbar_attached,
                        colorbar=cb_dict if not colorbar_attached else None,
                    )
                )
                colorbar_attached = True

    # 2. Branch top centerline outlines (Scatter3d lines) -- legacy append.
    #    When the legacy opaque ordering applies, the trace is appended here,
    #    after the curtain meshes (see block 0 for the transparent-curtain
    #    ordering where the trace is inserted before the meshes).
    if show_centerline and branch_x and not centerline_before_meshes:
        _add_centerline_trace()

    # 3. Ancestral Internal Nodes trace (Visible when show_node_markers=True, hoverable via invisible markers when False)
    internal_nodes = [n for n in plot_data.nodes.values() if not n.is_tip]
    if internal_nodes:
        if is_transformed:
            customdata_internal = [
                [
                    n.node_id,
                    "Ancestral Node",
                    n.raw_trait,
                    n.display_trait,
                    f"Descendants ({len(n.descendant_tips)} tips): {', '.join(n.descendant_tips[:3])}{'...' if len(n.descendant_tips) > 3 else ''}",
                ]
                for n in internal_nodes
            ]
            hovertemplate_internal = (
                "<b>Node: %{customdata[0]}</b><br>"
                "Type: %{customdata[1]}<br>"
                "Trait value (Y / Height): %{customdata[2]:.4f}<br>"
                "Time before present (Z): %{z:.4f}<br>"
                "Tree Layout (X): %{x:.2f}<br>"
                "%{customdata[4]}<extra></extra>"
            )
        else:
            customdata_internal = [
                [
                    n.node_id,
                    "Ancestral Node",
                    n.raw_trait,
                    f"Descendants ({len(n.descendant_tips)} tips): {', '.join(n.descendant_tips[:3])}{'...' if len(n.descendant_tips) > 3 else ''}",
                ]
                for n in internal_nodes
            ]
            hovertemplate_internal = (
                "<b>Node: %{customdata[0]}</b><br>"
                "Type: %{customdata[1]}<br>"
                "Trait value (Y / Height): %{customdata[2]:.4f}<br>"
                "Time before present (Z): %{z:.4f}<br>"
                "Tree Layout (X): %{x:.2f}<br>"
                "%{customdata[3]}<extra></extra>"
            )

        internal_marker_cfg = dict(
            size=internal_marker_size if show_node_markers else 12,
            color=[n.trait for n in internal_nodes],
            colorscale=plot_data.colorscale,
            cmin=trait_cmin,
            cmax=trait_cmax,
            reversescale=reverse_colorscale,
            symbol="diamond",
            opacity=0.95 if show_node_markers else 0.0,
        )

        fig.add_trace(
            go.Scatter3d(
                x=[n.x for n in internal_nodes],
                y=[n.y for n in internal_nodes],  # Y is TRAIT
                z=[n.z for n in internal_nodes],  # Z is TIME
                mode="markers",
                marker=internal_marker_cfg,
                customdata=customdata_internal,
                hovertemplate=hovertemplate_internal,
                name="Internal Nodes",
                showlegend=False,
            )
        )

    # 4. Terminal Nodes hover trace (Positions point indicator, XYZ spikes, and raw Hb buffer value at tree tips)
    tip_nodes = [n for n in plot_data.nodes.values() if n.is_tip]
    if tip_nodes:
        if is_transformed:
            customdata_tip = [
                [n.label, n.node_id, n.raw_trait, n.display_trait]
                for n in tip_nodes
            ]
            hovertemplate_tip = (
                "<b>Taxon: %{customdata[0]}</b><br>"
                "Node ID: %{customdata[1]}<br>"
                "Trait value (Y / Height): %{customdata[2]:.4f}<br>"
                "Time before present (Z): %{z:.4f}<br>"
                "Tree Layout (X): %{x:.2f}<extra></extra>"
            )
        else:
            customdata_tip = [
                [n.label, n.node_id, n.raw_trait]
                for n in tip_nodes
            ]
            hovertemplate_tip = (
                "<b>Taxon: %{customdata[0]}</b><br>"
                "Node ID: %{customdata[1]}<br>"
                "Trait value (Y / Height): %{customdata[2]:.4f}<br>"
                "Time before present (Z): %{z:.4f}<br>"
                "Tree Layout (X): %{x:.2f}<extra></extra>"
            )

        fig.add_trace(
            go.Scatter3d(
                x=[n.x for n in tip_nodes],
                y=[n.y for n in tip_nodes],  # Y is TRAIT
                z=[n.z for n in tip_nodes],  # Z is TIME
                mode="markers",
                marker=dict(
                    size=12,
                    color=[n.trait for n in tip_nodes],
                    colorscale=plot_data.colorscale,
                    cmin=trait_cmin,
                    cmax=trait_cmax,
                    reversescale=reverse_colorscale,
                    opacity=0.0,
                ),
                customdata=customdata_tip,
                hovertemplate=hovertemplate_tip,
                name="Terminal Taxa",
                showlegend=False,
            )
        )

    # 5. Fixed Terminal Species Labels trace (anchored strictly at present time
    #    Z = time_min / 0.0, constant baseline Y).
    #    Labels are pushed outward beyond the present plane (away from the tree
    #    body) by a generic fraction of the Time-before-present span, so the tip
    #    corresponds to the tree-facing end of the text rather than to its center.
    #    textposition='middle right' keeps the text body extending outward from
    #    that anchor. World coordinates are fixed: labels are never repositioned
    #    by camera changes or browser events.
    if tip_nodes and show_tip_labels:
        sorted_tips = sorted(tip_nodes, key=lambda n: n.x)
        time_span = abs(plot_data.time_max - plot_data.time_min)
        offset_fraction = (
            DEFAULT_TIP_LABEL_OFFSET_FRACTION
            if tip_label_offset is None
            else tip_label_offset
        )
        # The tree body extends from time_min (present) toward time_max (root),
        # so the outward direction beyond the present plane is away from the root.
        outward_sign = -1.0 if plot_data.time_max >= plot_data.time_min else 1.0
        label_time = plot_data.time_min + outward_sign * offset_fraction * time_span
        fig.add_trace(
            go.Scatter3d(
                x=[n.x for n in sorted_tips],
                y=[eff_baseline_y for _ in sorted_tips],
                z=[label_time for _ in sorted_tips],
                mode="text",
                text=[n.label for n in sorted_tips],
                textposition="middle right",
                textfont=dict(size=12, color="#222222"),
                hoverinfo="skip",
                name="Species Labels",
                showlegend=False,
            )
        )

    # 6. Configure Tree Layout axis (title removed, numeric ticks hidden, 3D spikes enabled)
    xaxis_cfg: Dict[str, Any] = dict(
        title=dict(text=""),
        showticklabels=False,
        ticks="",
        tickvals=[],
        ticktext=[],
        showbackground=False,
        gridcolor="#e5e5e5",
        zerolinecolor="#d0d0d0",
        showspikes=True,
        spikethickness=2,
        spikesides=True,
        spikecolor="#999999",
    )

    # Calculate default balanced manual aspect ratio.
    # trait_axis_scale multiplies ONLY the Trait (Y) dimension of the existing
    # aspect-ratio engine; Time (Z) and Tree Layout (X) visual aspects are
    # untouched, and no scientific coordinate is transformed.
    if aspect_ratio is None:
        ratio_x = 1.4
        ratio_y = 1.0
        ratio_z = 1.2
        ratio_dict = dict(x=ratio_x, y=ratio_y, z=ratio_z)
    else:
        ratio_dict = dict(aspect_ratio)
    ratio_dict["y"] = float(ratio_dict.get("y", 1.0)) * trait_axis_scale

    # Background colors
    is_transparent = background.lower() == "transparent"
    paper_bg = "rgba(0,0,0,0)" if is_transparent else "white"
    plot_bg = "rgba(0,0,0,0)" if is_transparent else "white"

    # Camera configuration
    if custom_camera is not None:
        camera_cfg = custom_camera
    else:
        camera_cfg = CAMERA_PRESETS.get(camera_preset.lower(), CAMERA_PRESETS["elife"])

    yaxis_cfg = dict(
        title=dict(text=y_axis_title, font=dict(size=13, color="#333333")),
        showbackground=False,
        gridcolor="#e5e5e5",
        zerolinecolor="#d0d0d0",
        showspikes=True,
        spikethickness=2,
        spikesides=True,
        spikecolor="#999999",
    )
    if is_transformed and y_tickvals:
        yaxis_cfg["tickmode"] = "array"
        yaxis_cfg["tickvals"] = y_tickvals
        yaxis_cfg["ticktext"] = y_ticktext

    # Scene and Camera Configuration
    fig.update_layout(
        title=dict(
            text=plot_title,
            x=0.5,
            xanchor="center",
            font=dict(size=18, color="#222222"),
        ),
        paper_bgcolor=paper_bg,
        plot_bgcolor=plot_bg,
        scene=dict(
            xaxis=xaxis_cfg,
            yaxis=yaxis_cfg,
            zaxis=dict(
                title=dict(text="Time before present", font=dict(size=13, color="#333333")),
                showbackground=False,
                gridcolor="#e5e5e5",
                zerolinecolor="#d0d0d0",
                showspikes=True,
                spikethickness=2,
                spikesides=True,
                spikecolor="#999999",
            ),
            aspectmode="manual",
            aspectratio=ratio_dict,
            camera=camera_cfg,
            annotations=[],
        ),
        margin=dict(l=20, r=20, t=50, b=20),
    )

    return fig
