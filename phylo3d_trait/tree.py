"""Tree parsing, topological layout, stable ID generation, and continuous trait interpolation.

Implements:
- Descendant-tip-set-based stable internal node IDs
- 2D tree skeleton layout (X = tip order / centroids, Z = divergence time before present)
- Fail-loud validation for tip and ancestral trait values
- Continuous branch interpolation (Y = linear trait interpolation)
"""

from __future__ import annotations

import hashlib
import io
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from phylo3d_trait.models import AnnotatedNode, EdgeSegment, PlotData


def compute_stable_node_id(descendant_tips: Iterable[str]) -> str:
    """Compute a deterministic, stable identifier for an internal clade.

    The ID is derived strictly from the sorted set of leaf/tip taxon names.
    This guarantees that topological re-orderings (e.g. (A,B) vs (B,A)) produce
    the identical stable ID.

    Args:
        descendant_tips: Collection of tip taxon names under this node.

    Returns:
        String of the format 'clade:<sha256_prefix>'.
    """
    sorted_tips = sorted(set(tip.strip() for tip in descendant_tips if tip.strip()))
    if not sorted_tips:
        raise ValueError("Cannot compute stable node ID for an empty tip set.")
    key = ",".join(sorted_tips).encode("utf-8")
    digest = hashlib.sha256(key).hexdigest()[:12]
    return f"clade:{digest}"


def parse_tree(
    tree_input: Union[str, Path, io.StringIO, io.TextIOBase],
    tree_format: Optional[str] = None,
) -> Tree:
    """Parse a phylogenetic tree from Newick or Nexus format.

    Args:
        tree_input: File path, tree string, or stream.
        tree_format: Explicit format ('newick' or 'nexus'). If None, auto-detected.

    Returns:
        Bio.Phylo Tree object.
    """
    if isinstance(tree_input, (str, Path)):
        p = Path(tree_input)
        if p.exists() and p.is_file():
            content = p.read_text(encoding="utf-8").strip()
        else:
            content = str(tree_input).strip()
    else:
        content = tree_input.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        content = content.strip()

    if tree_format is None:
        if content.upper().startswith("#NEXUS"):
            tree_format = "nexus"
        else:
            tree_format = "newick"

    stream = io.StringIO(content)
    try:
        tree = Phylo.read(stream, tree_format)
    except Exception as e:
        # Fallback format try
        alt_format = "newick" if tree_format == "nexus" else "nexus"
        stream.seek(0)
        try:
            tree = Phylo.read(stream, alt_format)
        except Exception:
            raise ValueError(f"Failed to parse phylogenetic tree as {tree_format} (or {alt_format}): {e}") from e

    return tree


def _collect_clade_descendant_tips(clade: Clade) -> List[str]:
    """Recursively collect all leaf tip names under a clade."""
    if clade.is_terminal():
        name = clade.name.strip() if clade.name else "unnamed_tip"
        return [name]
    tips: List[str] = []
    for child in clade.clades:
        tips.extend(_collect_clade_descendant_tips(child))
    return sorted(set(tips))


def _ensure_branch_lengths(tree: Tree) -> None:
    """Ensure all clades in tree have a valid numeric branch length."""
    for clade in tree.find_clades():
        if clade.branch_length is None or clade.branch_length < 0:
            clade.branch_length = 1.0


def compute_tree_layout(tree: Tree) -> Tuple[Dict[Clade, float], Dict[Clade, float], Dict[Clade, str], Dict[Clade, List[str]]]:
    """Calculate 2D skeleton coordinates (X = layout, Z = divergence time before present).

    X coordinates:
      - Terminals get X = 0, 1, ..., N-1 in leaf traversal order.
      - Internal nodes get X = mean(children X).

    Z coordinates:
      - For dated/ultrametric trees: tips have Z = 0.0, internal nodes have Z > 0.0,
        root has Z = max_age (time before present).

    Returns:
        Tuple of (x_coords, z_coords, node_ids, descendant_map).
    """
    _ensure_branch_lengths(tree)

    descendant_map: Dict[Clade, List[str]] = {}
    node_ids: Dict[Clade, str] = {}

    for clade in tree.find_clades():
        tips = _collect_clade_descendant_tips(clade)
        descendant_map[clade] = tips
        if clade.is_terminal():
            node_ids[clade] = clade.name.strip() if clade.name else tips[0]
        else:
            node_ids[clade] = compute_stable_node_id(tips)

    # 1. Compute X layout coordinates
    terminals = tree.get_terminals()
    x_coords: Dict[Clade, float] = {}
    for idx, term in enumerate(terminals):
        x_coords[term] = float(idx)

    # Post-order traversal to compute internal centroids
    def _assign_internal_x(c: Clade) -> float:
        if c in x_coords:
            return x_coords[c]
        child_xs = [_assign_internal_x(child) for child in c.clades]
        x_val = sum(child_xs) / len(child_xs) if child_xs else 0.0
        x_coords[c] = x_val
        return x_val

    _assign_internal_x(tree.root)

    # 2. Compute Z coordinates (Divergence time before present)
    # Calculate root-to-node distances
    root = tree.root
    root_dists: Dict[Clade, float] = {root: 0.0}

    def _calc_root_dists(c: Clade, current_dist: float) -> None:
        root_dists[c] = current_dist
        for child in c.clades:
            bl = child.branch_length if child.branch_length is not None else 1.0
            _calc_root_dists(child, current_dist + bl)

    _calc_root_dists(root, 0.0)

    max_dist_to_tips = max((root_dists[term] for term in terminals), default=1.0)
    z_coords: Dict[Clade, float] = {}

    for clade in tree.find_clades():
        dist = root_dists[clade]
        # Z is time before present (tips = 0, root = max_dist_to_tips)
        # Note: if tree has slight non-ultrametric tip distances, clamp tip to 0
        if clade.is_terminal():
            z_coords[clade] = 0.0
        else:
            time_val = max(0.0, max_dist_to_tips - dist)
            z_coords[clade] = time_val

    # Explicitly ensure root z equals max_dist_to_tips
    z_coords[root] = max_dist_to_tips

    return x_coords, z_coords, node_ids, descendant_map


def annotate_tree(
    tree: Tree,
    trait_values: Dict[str, float],
    num_segments: int = 10,
    colorscale: str = "Turbo",
    title: str = "3D Phylogenetic Tree with Continuous Trait Evolution",
) -> PlotData:
    """Map trait values onto tree nodes, validate completeness, and interpolate branch segments.

    Enforces:
    - X = Tree layout
    - Y = Trait value (Height)
    - Z = Evolutionary time before present
    - Linear continuous trait interpolation along branches
    - Global trait normalization across all nodes and segments
    - Loud failure on any missing tip or ancestral trait value

    Args:
        tree: Parsed Bio.Phylo Tree.
        trait_values: Dictionary mapping node IDs (tip names or 'clade:<hash>') to numeric traits.
        num_segments: Number of linear interpolation segments per branch (default 10).
        colorscale: Plotly continuous colorscale name.
        title: Plot title.

    Returns:
        PlotData container with annotated nodes and subdivided edge segments.

    Raises:
        ValueError: If any node lacks an explicit trait value or trait values are non-numeric.
    """
    x_coords, z_coords, node_ids, descendant_map = compute_tree_layout(tree)

    # Validate that all nodes have trait values provided
    missing_nodes: List[Tuple[str, str, List[str]]] = []
    resolved_traits: Dict[Clade, float] = {}

    for clade in tree.find_clades():
        nid = node_ids[clade]
        trait_val: Optional[float] = None

        # Check by stable node ID
        if nid in trait_values:
            trait_val = trait_values[nid]
        # Check by clade name (if internal node had a named label in tree)
        elif clade.name and clade.name.strip() in trait_values:
            trait_val = trait_values[clade.name.strip()]

        if trait_val is None:
            node_type = "Tip" if clade.is_terminal() else "Internal Node"
            missing_nodes.append((nid, node_type, descendant_map[clade]))
        else:
            try:
                resolved_traits[clade] = float(trait_val)
            except (ValueError, TypeError) as err:
                raise ValueError(f"Invalid non-numeric trait value '{trait_val}' for node '{nid}': {err}") from err

    if missing_nodes:
        sample_missing = "\n".join(
            f"  - [{ntype}] ID: {nid} (Descendant tips: {', '.join(tips[:4])}{'...' if len(tips) > 4 else ''})"
            for nid, ntype, tips in missing_nodes[:10]
        )
        extra = f"\n  ... and {len(missing_nodes) - 10} more" if len(missing_nodes) > 10 else ""
        raise ValueError(
            f"Trait values validation failed! Missing trait values for {len(missing_nodes)} node(s).\n"
            f"All tips and ancestral internal nodes MUST have explicit trait values.\n"
            f"{sample_missing}{extra}\n\n"
            f"Hint: Run 'phylo3d-trait template-values --tree <tree_file> --output <template.csv>' "
            f"to generate the complete list of required node IDs."
        )

    # Compute global trait and time boundaries
    all_traits = list(resolved_traits.values())
    trait_min = min(all_traits)
    trait_max = max(all_traits)

    all_times = list(z_coords.values())
    time_min = min(all_times)
    time_max = max(all_times)

    all_xs = list(x_coords.values())
    x_min = min(all_xs)
    x_max = max(all_xs)

    # Build parent-child mapping
    parent_map: Dict[Clade, Clade] = {}
    for parent in tree.find_clades():
        for child in parent.clades:
            parent_map[child] = parent

    # Build AnnotatedNode dictionary
    annotated_nodes: Dict[str, AnnotatedNode] = {}
    for clade in tree.find_clades():
        nid = node_ids[clade]
        trait_val = resolved_traits[clade]
        time_val = z_coords[clade]
        x_val = x_coords[clade]

        parent_clade = parent_map.get(clade)
        parent_id = node_ids[parent_clade] if parent_clade is not None else None
        children_ids = [node_ids[child] for child in clade.clades]

        label = clade.name.strip() if (clade.name and clade.name.strip()) else nid

        node = AnnotatedNode(
            node_id=nid,
            label=label,
            is_tip=clade.is_terminal(),
            x=x_val,
            y=trait_val,  # Y is TRAIT
            z=time_val,   # Z is TIME
            trait=trait_val,
            time=time_val,
            parent_id=parent_id,
            children_ids=children_ids,
            branch_length=clade.branch_length,
            descendant_tips=descendant_map[clade],
        )
        annotated_nodes[nid] = node

    # Subdivide branches into rectangular orthogonal subsegments (connector + lineage)
    segments: List[EdgeSegment] = []
    num_seg = max(1, int(num_segments))

    for child_clade, parent_clade in parent_map.items():
        parent_node = annotated_nodes[node_ids[parent_clade]]
        child_node = annotated_nodes[node_ids[child_clade]]

        xp, yp, zp = parent_node.x, parent_node.y, parent_node.z
        xc, yc, zc = child_node.x, child_node.y, child_node.z

        # 1. Connector subsegment: (xp, yp, zp) -> (xc, yp, zp)
        # Horizontal along X layout at constant parent time Zp and constant parent trait Yp
        n_conn = num_seg if xp != xc else 0
        n_lineage = num_seg
        total_edge_segments = n_conn + n_lineage

        current_idx = 0

        if n_conn > 0:
            for k in range(n_conn):
                t0 = k / float(n_conn)
                t1 = (k + 1) / float(n_conn)

                seg_x0 = xp + t0 * (xc - xp)
                seg_x1 = xp + t1 * (xc - xp)
                seg_y0 = yp
                seg_y1 = yp
                seg_z0 = zp
                seg_z1 = zp

                segment = EdgeSegment(
                    parent_id=parent_node.node_id,
                    child_id=child_node.node_id,
                    x0=seg_x0,
                    y0=seg_y0,
                    z0=seg_z0,
                    x1=seg_x1,
                    y1=seg_y1,
                    z1=seg_z1,
                    trait0=seg_y0,
                    trait1=seg_y1,
                    segment_index=current_idx,
                    total_segments=total_edge_segments,
                    segment_type="connector",
                )
                segments.append(segment)
                current_idx += 1

        # 2. Lineage subsegment: (xc, yp, zp) -> (xc, yc, zc)
        # Through-time along Z at constant child layout Xc; trait interpolates Yp -> Yc
        for m in range(n_lineage):
            t0 = m / float(n_lineage)
            t1 = (m + 1) / float(n_lineage)

            seg_x0 = xc
            seg_x1 = xc
            seg_y0 = yp + t0 * (yc - yp)
            seg_y1 = yp + t1 * (yc - yp)
            seg_z0 = zp + t0 * (zc - zp)
            seg_z1 = zp + t1 * (zc - zp)

            segment = EdgeSegment(
                parent_id=parent_node.node_id,
                child_id=child_node.node_id,
                x0=seg_x0,
                y0=seg_y0,
                z0=seg_z0,
                x1=seg_x1,
                y1=seg_y1,
                z1=seg_z1,
                trait0=seg_y0,
                trait1=seg_y1,
                segment_index=current_idx,
                total_segments=total_edge_segments,
                segment_type="lineage",
            )
            segments.append(segment)
            current_idx += 1

    return PlotData(
        nodes=annotated_nodes,
        segments=segments,
        trait_min=trait_min,
        trait_max=trait_max,
        time_min=time_min,
        time_max=time_max,
        x_min=x_min,
        x_max=x_max,
        colorscale=colorscale,
        title=title,
    )
