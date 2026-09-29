"""Feature A: terminal species labels aligned outward from the last tip.

Design contract under test:
- Labels stay on the fixed Trait plane (``baseline_y``) and at the present end
  of the Time-before-present axis.
- Labels are pushed *outward* (beyond ``time_min``, away from the tree body)
  along the Time-before-present axis by a generic fraction of the currently
  displayed time span (no magic absolute values, no per-dataset constants).
- ``textposition='middle right'`` makes the text body extend outward from the
  anchor, so the tree-facing (near) end of the label corresponds to the tip,
  with a small visible gap.
- Tree Layout (X) coordinates stay exactly equal to the terminal tip X.
- Tip ordering must not change.
- No camera callback / relayout JS may reposition labels.
"""

from __future__ import annotations

import inspect

import pytest

import phylo3d_trait.renderer as renderer_module
from phylo3d_trait.renderer import (
    DEFAULT_TIP_LABEL_OFFSET_FRACTION,
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


def _labels_trace(fig):
    traces = [t for t in fig.data if t.name == "Species Labels"]
    assert len(traces) == 1
    return traces[0]


def _sorted_tips(plot_data):
    return sorted(
        (n for n in plot_data.nodes.values() if n.is_tip), key=lambda n: n.x
    )


def test_a1_label_outside_tip():
    """TEST A1 / LABEL_OUTSIDE_TIP: label anchors no longer coincide with tips."""
    plot_data = _plot_data()
    fig = build_figure(plot_data)
    labels = _labels_trace(fig)

    span = plot_data.time_max - plot_data.time_min
    expected_z = plot_data.time_min - DEFAULT_TIP_LABEL_OFFSET_FRACTION * span
    assert expected_z < plot_data.time_min

    for idx, tip in enumerate(_sorted_tips(plot_data)):
        assert labels.z[idx] == pytest.approx(expected_z)
        # The label anchor is strictly outside the tree body (present side).
        assert labels.z[idx] < tip.z
        # The anchor no longer sits on the tip text-center position.
        assert (labels.x[idx], labels.z[idx]) != (tip.x, tip.z)


def test_a2_tip_matches_near_text_end():
    """TEST A2 / TIP_MATCHES_NEAR_TEXT_END: text body extends outward from the tip end."""
    plot_data = _plot_data()
    fig = build_figure(plot_data)
    labels = _labels_trace(fig)

    # The text hangs to the right of the anchor; combined with an outward
    # (present-side) anchor offset, the leftmost/tree-facing text end is the
    # end that corresponds to the tip, not the text center.
    assert labels.textposition == "middle right"

    span = plot_data.time_max - plot_data.time_min
    expected_offset = DEFAULT_TIP_LABEL_OFFSET_FRACTION * span
    assert expected_offset > 0
    for z in labels.z:
        assert plot_data.time_min - z == pytest.approx(expected_offset)


def test_a3_tip_layout_mapping():
    """TEST A3 / TIP_LAYOUT_MAPPING: label Tree Layout X still equals its tip X."""
    plot_data = _plot_data()
    fig = build_figure(plot_data)
    labels = _labels_trace(fig)
    tips = _sorted_tips(plot_data)

    assert len(labels.x) == len(tips)
    assert list(labels.text) == [tip.label for tip in tips]
    for idx, tip in enumerate(tips):
        assert labels.x[idx] == pytest.approx(tip.x)


def test_a4_tip_order_unchanged():
    """TEST A4 / TIP_ORDER_UNCHANGED: label ordering, uniqueness and completeness hold."""
    plot_data = _plot_data()
    fig = build_figure(plot_data)
    labels = _labels_trace(fig)
    tips = _sorted_tips(plot_data)

    assert list(labels.text) == ["A", "B", "C", "D"]
    assert len(labels.text) == len(set(labels.text)) == len(tips)
    xs = list(labels.x)
    assert xs == sorted(xs)


def test_a5_present_side_fixed():
    """TEST A5 / PRESENT_SIDE_FIXED: labels share one fixed present-side coordinate."""
    plot_data = _plot_data()
    fig = build_figure(plot_data)
    labels = _labels_trace(fig)

    assert len(set(labels.z)) == 1
    assert labels.z[0] <= plot_data.time_min
    # Fixed Trait plane, not per-tip trait height.
    assert len(set(labels.y)) == 1
    assert labels.y[0] == pytest.approx(plot_data.baseline_y)

    # Custom baseline plane keeps all labels on that same plane.
    fig_custom = build_figure(plot_data, baseline_y=4.0)
    labels_custom = _labels_trace(fig_custom)
    assert len(set(labels_custom.y)) == 1
    assert labels_custom.y[0] == pytest.approx(4.0)
    assert labels_custom.z[0] == pytest.approx(labels.z[0])


def test_a6_no_camera_dependent_reposition():
    """TEST A6 / NO_CAMERA_DEPENDENT_REPOSITION: coordinates are camera-independent."""
    plot_data = _plot_data()
    snapshots = []
    for preset in ("elife", "root_front", "tips_front"):
        labels = _labels_trace(build_figure(plot_data, camera_preset=preset))
        snapshots.append(
            (list(labels.x), list(labels.y), list(labels.z), labels.textposition)
        )
    assert all(snapshot == snapshots[0] for snapshot in snapshots)

    source = inspect.getsource(renderer_module).lower()
    for forbidden in ("addeventlistener", "plotly_relayout", "camerachanged", "on_camera"):
        assert forbidden not in source

    html = build_figure(plot_data).to_html(include_plotlyjs="cdn", full_html=False)
    assert "plotly_relayout" not in html.lower()
    assert "addeventlistener" not in html.lower()


def test_a7_label_trace_no_hover():
    """TEST A7 / LABEL_TRACE_NO_HOVER: label trace never steals hover."""
    plot_data = _plot_data()
    labels = _labels_trace(build_figure(plot_data))

    assert labels.hoverinfo == "skip"
    assert labels.mode == "text"
    assert labels.showlegend is False or labels.showlegend is None
    assert getattr(labels, "hovertemplate", None) in (None, "")


def test_tip_label_offset_parameter_is_generic_fraction_of_time_span():
    """The exposed offset parameter scales with the tree time span (generic)."""
    plot_data = _plot_data()
    span = plot_data.time_max - plot_data.time_min

    labels = _labels_trace(build_figure(plot_data, tip_label_offset=0.10))
    for z in labels.z:
        assert z == pytest.approx(plot_data.time_min - 0.10 * span)

    # A different tree scale yields a proportionally different distance,
    # never a hard-coded absolute value.
    plot_data_big = _plot_data()
    plot_data_big.time_max = plot_data_big.time_max * 10.0
    labels_big = _labels_trace(build_figure(plot_data_big, tip_label_offset=0.10))
    assert labels_big.z[0] == pytest.approx(plot_data.time_min - 0.10 * (plot_data_big.time_max - plot_data_big.time_min))
    assert abs(labels_big.z[0]) > abs(labels.z[0])


def test_tip_label_offset_zero_keeps_present_plane():
    """Offset 0 legitimately keeps labels exactly on the present plane."""
    plot_data = _plot_data()
    labels = _labels_trace(build_figure(plot_data, tip_label_offset=0.0))
    for z in labels.z:
        assert z == pytest.approx(plot_data.time_min)


@pytest.mark.parametrize("bad", [-0.1, float("nan"), float("inf"), float("-inf")])
def test_invalid_tip_label_offset_rejected(bad):
    """Negative / non-finite offsets fail loudly, never silently clamped."""
    plot_data = _plot_data()
    with pytest.raises(ValueError):
        build_figure(plot_data, tip_label_offset=bad)
