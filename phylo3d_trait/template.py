"""Template generation for node trait values.

Generates a template CSV file containing all tip taxa and internal clades
with their deterministic stable IDs (e.g. 'clade:<hash>') and descendant lists,
allowing users to fill in ancestral and tip trait values.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, List, Optional, Union

from Bio.Phylo.BaseTree import Clade, Tree
from phylo3d_trait.tree import _collect_clade_descendant_tips, compute_stable_node_id, parse_tree


def generate_template_csv(
    tree_input: Any,
    output_path: Union[str, Path],
    default_value: str = "",
) -> Path:
    """Generate a template CSV file listing all nodes in the phylogenetic tree.

    Args:
        tree_input: File path, tree string, or Tree object.
        output_path: Path to save the template CSV file.
        default_value: Default placeholder for the trait column.

    Returns:
        Path to the generated CSV file.
    """
    if isinstance(tree_input, (Tree, Clade)) or hasattr(tree_input, "get_terminals"):
        tree = tree_input
    else:
        tree = parse_tree(tree_input)

    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []

    # First collect tips
    for term in tree.get_terminals():
        tip_name = term.name.strip() if term.name else "unnamed_tip"
        rows.append(
            {
                "node_id": tip_name,
                "label": tip_name,
                "node_type": "tip",
                "descendant_count": 1,
                "descendant_tips": tip_name,
                "trait": default_value,
            }
        )

    # Next collect internal clades
    for clade in tree.get_nonterminals():
        tips = _collect_clade_descendant_tips(clade)
        stable_id = compute_stable_node_id(tips)
        label = clade.name.strip() if (clade.name and clade.name.strip()) else stable_id
        is_root = clade == tree.root
        node_type = "root" if is_root else "internal"
        rows.append(
            {
                "node_id": stable_id,
                "label": label,
                "node_type": node_type,
                "descendant_count": len(tips),
                "descendant_tips": ";".join(tips),
                "trait": default_value,
            }
        )

    fieldnames = [
        "node_id",
        "label",
        "node_type",
        "descendant_count",
        "descendant_tips",
        "trait",
    ]

    with open(out_file, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return out_file
