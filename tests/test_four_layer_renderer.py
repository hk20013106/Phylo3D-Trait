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

