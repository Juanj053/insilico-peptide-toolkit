"""JSON serialization for peptide analytical profiles and screening dossiers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Union


def export_to_json(
    data: Union[Dict[str, Any], List[Dict[str, Any]]],
    output_path: Union[str, Path],
    indent: int = 2,
) -> None:
    """
    Save dictionary or list of profiles to JSON file.

    Args:
        data: Serializable dictionary or list.
        output_path: Path to output JSON file.
        indent: Indentation spaces.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=indent, default=str)
