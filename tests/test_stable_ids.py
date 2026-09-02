"""Unit tests for deterministic stable clade IDs under various tree permutations."""

from phylo3d_trait.tree import compute_stable_node_id, parse_tree, _collect_clade_descendant_tips


def test_sister_taxa_rotation():
    """Verify that rotating sister branches preserves exact internal node IDs."""
    tree1_str = "((A:1,B:1):2,(C:1,D:1):2);"
    tree2_str = "((B:1,A:1):2,(D:1,C:1):2);"
    tree3_str = "((C:1,D:1):2,(B:1,A:1):2);"

    t1 = parse_tree(tree1_str)
    t2 = parse_tree(tree2_str)
    t3 = parse_tree(tree3_str)

    def get_clade_ids(tree):
        ids = set()
        for c in tree.get_nonterminals():
            tips = _collect_clade_descendant_tips(c)
            ids.add(compute_stable_node_id(tips))
        return ids

    ids1 = get_clade_ids(t1)
    ids2 = get_clade_ids(t2)
    ids3 = get_clade_ids(t3)

    assert ids1 == ids2
    assert ids1 == ids3


def test_polytomy_clade_id():
    """Verify stable clade ID works on polytomies with 3+ children."""
    tips = ["Taxon_Alpha", "Taxon_Beta", "Taxon_Gamma"]
    id1 = compute_stable_node_id(tips)
    id2 = compute_stable_node_id(["Taxon_Gamma", "Taxon_Alpha", "Taxon_Beta"])
    assert id1 == id2
