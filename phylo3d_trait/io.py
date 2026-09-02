"""Input/output utilities for reading trait values and tree files."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Optional, Union


def load_trait_values(values_path: Union[str, Path]) -> Dict[str, float]:
    """Load trait values from a CSV or TSV file.

    Supports headers such as 'node_id', 'id', 'label', 'taxon', 'name', 'clade' for node identifiers,
    and 'trait', 'trait_value', 'value', 'y', 'height' for trait values.
    If headers are unrecognized or missing, assumes the first column is node ID and second is trait value.

    Args:
        values_path: Path to the CSV or TSV file.

    Returns:
        Dictionary mapping node ID string to float trait value.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If file is empty or lines cannot be parsed as numeric floats.
    """
    path = Path(values_path)
    if not path.exists():
        raise FileNotFoundError(f"Trait values file not found: {path}")

    delimiter = "\t" if path.suffix.lower() in [".tsv", ".tab"] else ","

    trait_dict: Dict[str, float] = {}

    with open(path, mode="r", encoding="utf-8-sig") as f:
        # Detect delimiter if possible
        sample = f.read(2048)
        f.seek(0)
        if delimiter != "\t" and "\t" in sample and "," not in sample:
            delimiter = "\t"

        reader = csv.reader(f, delimiter=delimiter)
        rows = [row for row in reader if row and any(cell.strip() for cell in row)]

    if not rows:
        raise ValueError(f"Trait values file is empty: {path}")

    # Check header with priority ordering
    header = [col.strip().lower() for col in rows[0]]
    id_col_idx: Optional[int] = None
    trait_col_idx: Optional[int] = None
    label_col_idx: Optional[int] = None

    id_priority = ["node_id", "id", "node", "taxon", "taxa", "species", "clade", "name", "label"]
    for alias in id_priority:
        if alias in header:
            id_col_idx = header.index(alias)
            break

    if "label" in header:
        label_col_idx = header.index("label")

    trait_priority = ["trait", "trait_value", "traitvalue", "value", "val", "y", "height"]
    for alias in trait_priority:
        if alias in header:
            trait_col_idx = header.index(alias)
            break

    has_header = id_col_idx is not None or trait_col_idx is not None
    if id_col_idx is None:
        id_col_idx = 0
    if trait_col_idx is None:
        trait_col_idx = 1 if len(header) > 1 else 0

    start_idx = 1 if has_header else 0

    for line_num, row in enumerate(rows[start_idx:], start=start_idx + 1):
        if not row:
            continue
        # Skip comment rows
        if row[0].strip().startswith("#"):
            continue

        if len(row) <= max(id_col_idx, trait_col_idx):
            continue

        raw_id = row[id_col_idx].strip()
        raw_val = row[trait_col_idx].strip()

        if not raw_id or not raw_val:
            continue

        try:
            val = float(raw_val)
        except ValueError as e:
            raise ValueError(
                f"Non-numeric trait value '{raw_val}' at row {line_num} in '{path}': {e}"
            ) from e

        trait_dict[raw_id] = val
        if label_col_idx is not None and len(row) > label_col_idx:
            raw_label = row[label_col_idx].strip()
            if raw_label and raw_label not in trait_dict:
                trait_dict[raw_label] = val

    if not trait_dict:
        raise ValueError(f"No valid trait entries could be parsed from '{path}'.")

    return trait_dict
