"""Feature B: --trait-axis-scale changes only the Trait visual aspect ratio.

Contract under test:
- ``trait_axis_scale`` defaults to 1.0 and is fully backward compatible.
- It multiplies ONLY the Trait (Y) component of the Plotly scene aspect ratio.
- Time-before-present (Z) and Tree Layout (X) visual aspects are untouched.
- No scientific coordinate, tick, tooltip, color or mesh data is modified.
- ``--trait-display-offset`` (geometric zero shift) and ``--trait-axis-scale``
  (visual aspect scaling) are independent and compose.
- 0 / negative / NaN / inf are rejected loudly (never silently clamped).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from phylo3d_trait.cli import main
from phylo3d_trait.renderer import build_figure, build_plot_data
from phylo3d_trait.tree import compute_stable_node_id, parse_tree

TREE_STR = "((A:10,B:10):20,(C:15,D:15):15);"

DEFAULT_ASPECT = {"x": 1.4, "y": 1.0, "z": 1.2}


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


def _aspect(fig):
    ar = fig.layout.scene.aspectratio
    return {"x": float(ar.x), "y": float(ar.y), "z": float(ar.z)}


def _mesh(fig):
    return [t for t in fig.data if t.name == "Branch Curtains"][0]


def test_b1_default_scale():
    """TEST B1 / DEFAULT_SCALE: omitted or 1.0 keeps the historical aspect ratio."""
    plot_data = _plot_data()
    assert _aspect(build_figure(plot_data)) == DEFAULT_ASPECT
    assert _aspect(build_figure(plot_data, trait_axis_scale=1.0)) == DEFAULT_ASPECT

    fig_default = build_figure(plot_data)
    fig_explicit = build_figure(plot_data, trait_axis_scale=1.0)
    assert fig_default.to_dict() == fig_explicit.to_dict()


def test_b2_half_scale_only_trait_dimension():
    """TEST B2 / HALF_SCALE: 0.5 halves only the Trait (Y) aspect ratio."""
    plot_data = _plot_data()
    aspect = _aspect(build_figure(plot_data, trait_axis_scale=0.5))

    assert aspect["y"] == pytest.approx(DEFAULT_ASPECT["y"] * 0.5)
    assert aspect["x"] == pytest.approx(DEFAULT_ASPECT["x"])
    assert aspect["z"] == pytest.approx(DEFAULT_ASPECT["z"])


def test_b3_data_unchanged_by_scale():
    """TEST B3 / DATA_UNCHANGED: only scene aspect ratio may differ."""
    plot_data = _plot_data()
    fig_1 = build_figure(plot_data, trait_axis_scale=1.0)
    fig_05 = build_figure(plot_data, trait_axis_scale=0.5)

    aspect_1 = _aspect(fig_1)
    aspect_05 = _aspect(fig_05)
    assert aspect_1["y"] / aspect_05["y"] == pytest.approx(2.0)

    dict_1 = fig_1.to_dict()
    dict_05 = fig_05.to_dict()
    dict_1["layout"]["scene"].pop("aspectratio")
    dict_05["layout"]["scene"].pop("aspectratio")
    assert dict_1 == dict_05

    # Explicit spot checks: mesh geometry, triangles, intensities, color domain,
    # color scale, tick values and labels are all identical.
    mesh_1, mesh_05 = _mesh(fig_1), _mesh(fig_05)
    for attr in ("x", "y", "z", "i", "j", "k", "intensity"):
        assert list(getattr(mesh_1, attr)) == list(getattr(mesh_05, attr))
    assert mesh_1.cmin == pytest.approx(mesh_05.cmin)
    assert mesh_1.cmax == pytest.approx(mesh_05.cmax)
    assert mesh_1.colorscale == mesh_05.colorscale

    y1 = fig_1.layout.scene.yaxis
    y05 = fig_05.layout.scene.yaxis
    assert y1.to_plotly_json() == y05.to_plotly_json()

    for tip in (t for t in fig_1.data if t.name == "Terminal Taxa"):
        tip_05 = [t for t in fig_05.data if t.name == "Terminal Taxa"][0]
        assert list(tip.customdata) == list(tip_05.customdata)


def test_b4_time_axis_visual_unchanged():
    """TEST B4 / TIME_AXIS_UNCHANGED: Time-before-present aspect is not scaled."""
    plot_data = _plot_data()
    aspect = _aspect(build_figure(plot_data, trait_axis_scale=0.5))
    assert aspect["z"] == pytest.approx(DEFAULT_ASPECT["z"])

    zaxis = build_figure(plot_data, trait_axis_scale=0.5).layout.scene.zaxis
    assert zaxis.title.text == "Time before present"


def test_b5_tree_layout_visual_unchanged():
    """TEST B5 / TREE_LAYOUT_UNCHANGED: Tree Layout aspect and axis config unchanged."""
    plot_data = _plot_data()
    aspect = _aspect(build_figure(plot_data, trait_axis_scale=0.5))
    assert aspect["x"] == pytest.approx(DEFAULT_ASPECT["x"])

    xaxis_1 = build_figure(plot_data, trait_axis_scale=1.0).layout.scene.xaxis
    xaxis_05 = build_figure(plot_data, trait_axis_scale=0.5).layout.scene.xaxis
    assert xaxis_1.to_plotly_json() == xaxis_05.to_plotly_json()


def test_b6_offset_and_scale_are_independent():
    """TEST B6 / OFFSET_AND_SCALE_INDEPENDENT: offset maps raw 10 -> display 6 at scale 0.5."""
    plot_data = _plot_data(trait_display_offset=4.0)
    fig = build_figure(plot_data, trait_axis_scale=0.5)

    assert plot_data.raw_to_display(10.0) == pytest.approx(6.0)
    assert plot_data.raw_to_display(4.0) == pytest.approx(0.0)
    assert plot_data.nodes["D"].raw_trait == pytest.approx(8.0)
    assert plot_data.nodes["D"].display_trait == pytest.approx(4.0)

    # Display geometry still uses the offset mapping, while aspect is scaled.
    assert _aspect(fig)["y"] == pytest.approx(0.5)
    yaxis = fig.layout.scene.yaxis
    for val, text in zip(yaxis.tickvals, yaxis.ticktext):
        assert float(text) == pytest.approx(plot_data.display_to_raw(val), abs=1e-2)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf"), float("-inf")])
def test_b7_invalid_scale_rejected(bad):
    """TEST B7 / INVALID_SCALE: 0, negative, NaN and inf are rejected loudly."""
    plot_data = _plot_data()
    with pytest.raises(ValueError):
        build_figure(plot_data, trait_axis_scale=bad)


@pytest.mark.parametrize("bad", ["0", "-0.5", "nan", "inf"])
def test_b7_invalid_scale_rejected_by_cli(tmp_path, bad):
    """TEST B7 / INVALID_SCALE: the CLI rejects the same invalid values."""
    repo_root = Path(__file__).parent.parent
    tree_file = repo_root / "examples" / "example1" / "tree.nwk"
    values_file = repo_root / "examples" / "example1" / "node_values.csv"
    out_html = tmp_path / "invalid_scale.html"

    with pytest.raises(SystemExit):
        main(
            [
                "plot",
                "--tree",
                str(tree_file),
                "--values",
                str(values_file),
                "--output",
                str(out_html),
                "--trait-axis-scale",
                bad,
            ]
        )
    assert not out_html.exists()


def test_cli_valid_scale_end_to_end(tmp_path):
    """The CLI accepts a valid scale and applies it to the rendered scene."""
    repo_root = Path(__file__).parent.parent
    tree_file = repo_root / "examples" / "example1" / "tree.nwk"
    values_file = repo_root / "examples" / "example1" / "node_values.csv"
    out_html = tmp_path / "scaled.html"

    code = main(
        [
            "plot",
            "--tree",
            str(tree_file),
            "--values",
            str(values_file),
            "--output",
            str(out_html),
            "--trait-axis-scale",
            "0.5",
        ]
    )
    assert code == 0
    content = out_html.read_text(encoding="utf-8")
    assert '"aspectratio"' in content
    assert '"y":0.5' in content or '"y": 0.5' in content
