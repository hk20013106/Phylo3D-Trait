import pytest

from phylo3d_trait.cli import build_parser
from phylo3d_trait.four_layer_renderer import (
    FOUR_LAYER_PEELS,
    _build_payload,
    _validate_opacity,
    write_four_layer_html,
)
from phylo3d_trait.renderer import build_plot_data
from phylo3d_trait.tree import compute_stable_node_id, parse_tree


def _data():
    tree = parse_tree("(A:10,B:10);")
    root = compute_stable_node_id(["A", "B"])
    return build_plot_data(
        tree,
        {"A": 5.0, "B": 9.0, root: 7.0},
        num_segments=4,
        baseline_y=0.0,
        trait_display_offset=4.0,
    )


def test_fixed_four_layer_contract_at_max_transparency():
    payload = _build_payload(
        _data(),
        opacity=0.5,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
    )

    assert FOUR_LAYER_PEELS == 4
    assert payload["peel_layers"] == 4
    assert payload["terminal_layer"] == 4
    assert payload["opacity"] == pytest.approx(0.5)
    assert payload["transparency"] == pytest.approx(0.5)
    assert payload["max_omitted_transmittance"] == pytest.approx(0.5 ** 4)
    assert payload["stats"]["triangles"] > 0
    assert len(payload["mesh"]["positions"]) // 3 == len(payload["mesh"]["colors"]) // 3


@pytest.mark.parametrize("value", [0.0, 0.49, -1.0, 1.01, float("inf")])
def test_four_layer_opacity_rejects_unsupported_values(value):
    with pytest.raises(ValueError):
        _validate_opacity(value)


@pytest.mark.parametrize("value", [0.5, 0.7, 0.9, 1.0])
def test_four_layer_opacity_accepts_supported_values(value):
    assert _validate_opacity(value) == pytest.approx(value)


def test_generated_html_uses_fragment_peeling_not_trace_sorting(tmp_path):
    out = tmp_path / "four_layer.html"
    payload = write_four_layer_html(
        _data(),
        out,
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        centerline_color="trait",
    )
    html = out.read_text(encoding="utf-8")

    assert payload["max_omitted_transmittance"] == pytest.approx(0.3 ** 4)
    assert "const MAX_PEEL_LAYERS=4;" in html
    assert "if(d.a>0.0)o=vec4(d.rgb,1.0)" in html
    assert "uPrevDepth" in html
    assert "gl_FragCoord.z" in html
    assert "Plotly.moveTraces" not in html
    assert "nearest-corner" not in html
    assert "gl.objects" not in html


def test_cli_exposes_four_layer_backend():
    parser = build_parser()
    args = parser.parse_args(
        [
            "plot",
            "--tree",
            "tree.nwk",
            "--values",
            "values.csv",
            "--output",
            "out.html",
            "--renderer",
            "four-layer",
            "--opacity",
            "0.7",
        ]
    )
    assert args.renderer == "four-layer"
    assert args.opacity == pytest.approx(0.7)


def test_generated_html_clears_samplers_before_peel_and_performs_symmetric_cleanup(tmp_path):
    out = tmp_path / "four_layer_cleanup.html"
    write_four_layer_html(
        _data(),
        out,
        opacity=0.7,
        baseline_y=0.0,
    )
    html = out.read_text(encoding="utf-8")

    # Guard against WebGL feedback loop (GL_INVALID_OPERATION 1282):
    # Must unbind TEXTURE0 before entering peel loop
    assert "gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);" in html
    # Must perform symmetric cleanup after line draw and at render exit
    assert "gl.bindFramebuffer(gl.FRAMEBUFFER,null);" in html
    assert "gl.bindVertexArray(null);" in html


def test_four_layer_axis_payload():
    payload = _build_payload(
        _data(),
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
    )

    assert "axes" in payload
    axes = payload["axes"]
    assert "y_axis" in axes
    assert "z_axis" in axes
    assert "x_axis" in axes
    assert axes["x_axis"]["show_numeric_labels"] is False

    y_axis = axes["y_axis"]
    assert y_axis["title"] == "Trait value"
    assert len(y_axis["ticks"]) >= 2
    # Verify Y raw labels are raw scientific traits (not corrupted by trait_display_offset 4.0)
    raw_texts = [t["text"] for t in y_axis["ticks"]]
    # In _data(), raw traits are 5.0, 7.0, 9.0; with baseline 0.0 (raw 4.0), raw ticks must include '4', '5', '9'
    assert "4" in raw_texts
    assert "5" in raw_texts
    assert "9" in raw_texts

    # For display Y = 1.0 (val == 1.0), text must be '5'
    tick_disp_1 = next((t for t in y_axis["ticks"] if abs(t["val"] - 1.0) < 1e-4), None)
    assert tick_disp_1 is not None
    assert tick_disp_1["text"] == "5"

    z_axis = axes["z_axis"]
    assert z_axis["title"] == "Time before present"
    assert len(z_axis["ticks"]) >= 2
    z_texts = [t["text"] for t in z_axis["ticks"]]
    # Root age is 10.0, present is 0.0
    assert "0" in z_texts
    assert "10" in z_texts


def test_four_layer_hover_payload():
    payload = _build_payload(
        _data(),
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
    )

    assert "nodes" in payload
    nodes = payload["nodes"]
    assert len(nodes) == 3

    tips = [n for n in nodes if n["is_tip"]]
    internals = [n for n in nodes if not n["is_tip"]]
    assert len(tips) == 2
    assert len(internals) == 1

    for tip in tips:
        assert "label" in tip and tip["label"] in ["A", "B"]
        assert "node_id" in tip and tip["node_id"] in ["A", "B"]
        assert "raw_trait" in tip and tip["raw_trait"] in [5.0, 9.0]
        assert "time" in tip and tip["time"] == pytest.approx(0.0)
        assert "x" in tip and isinstance(tip["x"], float)
        assert len(tip["position"]) == 3
        assert len(tip["baseline_pos"]) == 3

    internal = internals[0]
    assert "node_id" in internal and internal["node_id"].startswith("clade:")
    assert internal["raw_trait"] == pytest.approx(7.0)
    assert internal["time"] == pytest.approx(10.0)
    assert isinstance(internal["x"], float)
    assert "descendants" in internal
    assert internal["descendants"].startswith("Descendants (2 tips):")
    assert "A" in internal["descendants"] and "B" in internal["descendants"]


def test_four_layer_generated_html_static_features(tmp_path):
    out = tmp_path / "four_layer_features.html"
    write_four_layer_html(
        _data(),
        out,
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        centerline_color="trait",
    )
    html = out.read_text(encoding="utf-8")

    # Axes overlay
    assert '<svg id="axes">' in html
    assert "updateAxes" in html
    assert "svgYAxis" in html
    assert "svgZAxis" in html
    assert "frontAxisCorner" in html
    assert "outward2D" in html
    assert "svgYTicks" in html
    assert "svgZTicks" in html

    # Node picking and tooltip DOM
    assert '<div id="tooltip">' in html
    assert "updateHover" in html
    assert "projectedNodes" in html
    assert "hoverMarker" in html
    assert "hoverGuide" in html

    # Camera-aware label anchor
    assert "const eyeX=Math.cos(pitch)*Math.sin(yaw);" in html
    assert 'const tx=eyeX<0?"translate(calc(-100% - 3px),-50%)":"translate(3px,-50%)";' in html
    assert "e.style.transform=tx" in html

    # Guard against forbidden Plotly workarounds
    assert "Plotly.moveTraces" not in html
    assert "nearest-corner" not in html
    assert "gl.objects" not in html


def test_four_layer_axes_and_hover_visibility_config():
    data = _data()
    # 1. Defaults
    payload_default = _build_payload(
        data,
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
    )
    assert payload_default["axes"]["show_x"] is True
    assert payload_default["axes"]["show_y"] is True
    assert payload_default["axes"]["show_z"] is True
    assert payload_default["hover"]["show_tip"] is True
    assert payload_default["hover"]["show_internal"] is True

    # 2. Disabled
    payload_disabled = _build_payload(
        data,
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
        show_x_axis=False,
        show_y_axis=False,
        show_z_axis=False,
        show_tip_hover=False,
        show_internal_hover=False,
    )
    assert payload_disabled["axes"]["show_x"] is False
    assert payload_disabled["axes"]["show_y"] is False
    assert payload_disabled["axes"]["show_z"] is False
    assert payload_disabled["hover"]["show_tip"] is False
    assert payload_disabled["hover"]["show_internal"] is False


def test_four_layer_display_flags_do_not_alter_geometry_or_metadata():
    data = _data()
    p_on = _build_payload(
        data,
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
        show_x_axis=True,
        show_y_axis=True,
        show_z_axis=True,
        show_tip_hover=True,
        show_internal_hover=True,
    )
    p_off = _build_payload(
        data,
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
        show_x_axis=False,
        show_y_axis=False,
        show_z_axis=False,
        show_tip_hover=False,
        show_internal_hover=False,
    )

    # Geometry must be completely invariant
    assert p_on["mesh"] == p_off["mesh"]
    assert p_on["centerline"] == p_off["centerline"]
    assert p_on["stats"] == p_off["stats"]
    assert p_on["opacity"] == p_off["opacity"]
    assert p_on["camera_eye"] == p_off["camera_eye"]

    # Scientific node metadata and trait values must be invariant
    assert p_on["nodes"] == p_off["nodes"]





def test_four_layer_toolbar_exports_hybrid_svg_and_png(tmp_path):
    out = tmp_path / "four_layer_export.html"
    write_four_layer_html(
        _data(),
        out,
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        centerline_color="trait",
    )
    html = out.read_text(encoding="utf-8")

    assert 'id="toolbar"' in html
    assert 'id="reset-view"' in html
    assert 'id="download-png"' in html
    assert 'id="download-svg"' in html
    assert "function currentExportSvg()" in html
    assert 'canvas.toDataURL("image/png")' in html
    assert '<image href="' in html
    assert "function exportSvg()" in html
    assert "function exportPng()" in html


def test_four_layer_axes_do_not_draw_box_frame():
    payload = _build_payload(
        _data(),
        opacity=0.7,
        baseline_y=0.0,
        reverse_colorscale=True,
        curtain_color_mode="branch",
        trait_axis_scale=0.5,
        tip_label_offset=0.03,
        show_tip_labels=True,
        show_centerline=True,
        centerline_color="trait",
        background="white",
        camera_preset="elife",
    )
    axes = payload["axes"]
    assert "frame_lines" not in axes
    assert axes["y_axis"]["title"] == "Trait value"
    assert axes["z_axis"]["title"] == "Time before present"
    assert axes["x_axis"]["visible"] is False
