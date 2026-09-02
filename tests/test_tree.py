"""Unit tests for tree parsing, coordinate mapping, and branch interpolation."""

import pytest
from Bio.Phylo.BaseTree import Tree

from phylo3d_trait.models import AnnotatedNode, EdgeSegment, PlotData
from phylo3d_trait.tree import annotate_tree, compute_stable_node_id, parse_tree


def test_stable_node_id_invariance():
    """Verify stable clade ID is invariant to sister taxa permutation."""
    id1 = compute_stable_node_id(["A", "B"])
    id2 = compute_stable_node_id(["B", "A"])
    assert id1 == id2
    assert id1.startswith("clade:")

    # Multi-taxon clade
    id3 = compute_stable_node_id(["TaxonC", "TaxonA", "TaxonB"])
    id4 = compute_stable_node_id(["TaxonB", "TaxonC", "TaxonA"])
    assert id3 == id4


def test_tree_parsing_newick_and_nexus():
    """Verify parsing of Newick and Nexus strings."""
    nwk = "((A:5,B:5):10,(C:8,D:8):7);"
    tree_nwk = parse_tree(nwk)
    assert isinstance(tree_nwk, Tree)
    assert len(tree_nwk.get_terminals()) == 4

    nex = """#NEXUS
BEGIN TREES;
    TREE tree1 = ((A:5,B:5):10,(C:8,D:8):7);
END;"""
    tree_nex = parse_tree(nex)
    assert isinstance(tree_nex, Tree)
    assert len(tree_nex.get_terminals()) == 4


def test_trait_mapped_to_y_axis():
    """Verify trait values are strictly mapped to the Y axis (height)."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)

    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 10.5,
        "B": 20.0,
        "C": 30.2,
        "D": 40.8,
        id_ab: 15.0,
        id_cd: 35.0,
        id_root: 25.0,
    }

    plot_data = annotate_tree(tree, trait_values, num_segments=8)

    # Check each node: y MUST equal trait
    for node_id, node in plot_data.nodes.items():
        expected_trait = trait_values[node_id]
        assert node.y == expected_trait, f"Node {node_id} y={node.y} != trait={expected_trait}"
        assert node.trait == expected_trait
        assert node.y == node.trait


def test_time_mapped_to_z_axis():
    """Verify evolutionary time / divergence age is mapped to Z axis."""
    # Tree:
    # Root age: 30.0
    # Clade AB: age 10.0 (dist from root = 20.0, age = 30 - 20 = 10.0)
    # Clade CD: age 15.0 (dist from root = 15.0, age = 30 - 15 = 15.0)
    # Tips: age 0.0 (Z == 0)
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)

    id_ab = compute_stable_node_id(["A", "B"])
    id_cd = compute_stable_node_id(["C", "D"])
    id_root = compute_stable_node_id(["A", "B", "C", "D"])

    trait_values = {
        "A": 1.0,
        "B": 2.0,
        "C": 3.0,
        "D": 4.0,
        id_ab: 1.5,
        id_cd: 3.5,
        id_root: 2.5,
    }

    plot_data = annotate_tree(tree, trait_values)

    # Tips must have Z == 0.0
    for tip_name in ["A", "B", "C", "D"]:
        tip_node = plot_data.nodes[tip_name]
        assert tip_node.z == 0.0, f"Tip {tip_name} Z is {tip_node.z}, expected 0.0"
        assert tip_node.time == 0.0

    # Root must have Z == 30.0 (root age)
    root_node = plot_data.nodes[id_root]
    assert root_node.z == pytest.approx(30.0)
    assert root_node.time == pytest.approx(30.0)

    # Internal clades
    ab_node = plot_data.nodes[id_ab]
    assert ab_node.z == pytest.approx(10.0)
    cd_node = plot_data.nodes[id_cd]
    assert cd_node.z == pytest.approx(15.0)


def test_continuous_branch_linear_interpolation():
    """Verify branches are subdivided into rectangular orthogonal connector and lineage subsegments."""
    tree_str = "(A:10,B:10);"
    tree = parse_tree(tree_str)

    id_root = compute_stable_node_id(["A", "B"])

    # Parent (root) trait = 1.0, child (A) trait = 3.0
    trait_values = {
        "A": 3.0,
        "B": 5.0,
        id_root: 1.0,
    }

    num_segs = 10
    plot_data = annotate_tree(tree, trait_values, num_segments=num_segs)

    # Filter segments connecting root to A
    root_a_segs = [
        s for s in plot_data.segments if s.parent_id == id_root and s.child_id == "A"
    ]
    # Root (x=0.5) to A (x=0.0): xp != xc, so num_segs connector + num_segs lineage = 20
    assert len(root_a_segs) == 2 * num_segs

    # Sort by segment index
    root_a_segs.sort(key=lambda s: s.segment_index)

    connector_segs = [s for s in root_a_segs if s.segment_type == "connector"]
    lineage_segs = [s for s in root_a_segs if s.segment_type == "lineage"]

    assert len(connector_segs) == num_segs
    assert len(lineage_segs) == num_segs

    # 1. Verify connector properties (horizontal in X, constant Yp=1.0 and constant Zp=10.0)
    for c in connector_segs:
        assert c.y0 == pytest.approx(1.0)
        assert c.y1 == pytest.approx(1.0)
        assert c.z0 == pytest.approx(10.0)
        assert c.z1 == pytest.approx(10.0)
        assert c.trait0 == pytest.approx(1.0)
        assert c.trait1 == pytest.approx(1.0)

    assert connector_segs[0].x0 == pytest.approx(0.5)
    assert connector_segs[-1].x1 == pytest.approx(0.0)

    # 2. Verify lineage properties (through-time in Z from 10.0 to 0.0, constant Xc=0.0, Y interpolates 1.0 to 3.0)
    for l in lineage_segs:
        assert l.x0 == pytest.approx(0.0)
        assert l.x1 == pytest.approx(0.0)

    assert lineage_segs[0].z0 == pytest.approx(10.0)
    assert lineage_segs[-1].z1 == pytest.approx(0.0)
    assert lineage_segs[0].y0 == pytest.approx(1.0)
    assert lineage_segs[-1].y1 == pytest.approx(3.0)

    # Verify overall segment continuity
    for k in range(len(root_a_segs) - 1):
        assert root_a_segs[k].x1 == pytest.approx(root_a_segs[k + 1].x0)
        assert root_a_segs[k].y1 == pytest.approx(root_a_segs[k + 1].y0)
        assert root_a_segs[k].z1 == pytest.approx(root_a_segs[k + 1].z0)
        assert root_a_segs[k].trait1 == pytest.approx(root_a_segs[k + 1].trait0)


def test_missing_ancestral_values_fails_loudly():
    """Verify that missing ancestral or tip trait values raise a loud ValueError."""
    tree_str = "((A:10,B:10):20,(C:15,D:15):15);"
    tree = parse_tree(tree_str)

    # Missing ancestral nodes and tip D
    incomplete_traits = {
        "A": 1.0,
        "B": 2.0,
        "C": 3.0,
    }

    with pytest.raises(ValueError) as exc_info:
        annotate_tree(tree, incomplete_traits)

    err_msg = str(exc_info.value)
    assert "Trait values validation failed" in err_msg
    assert "Missing trait values" in err_msg
    assert "D" in err_msg
