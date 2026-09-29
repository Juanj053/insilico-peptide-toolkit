"""
ChEMBL bioactivity and pharmacological target database connector with 3-tier fallback.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from insilico_peptide.connectors.base import BaseConnector


class ChEMBLConnector(BaseConnector):
    """Client for EMBL-EBI ChEMBL REST API with offline caching and curated fixtures."""

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/target"

    def fetch(self, target_chembl_id: str, offline: bool = False) -> Dict[str, Any]:
        """
        Fetch ChEMBL target data and bioactivities by target ChEMBL ID (e.g. 'CHEMBL1784', 'CHEMBL214').

        Args:
            target_chembl_id: Target identifier (e.g. 'CHEMBL1784').
            offline: Force offline mode.

        Returns:
            Normalized dictionary with target metadata and pharmacological bioactivities.
        """
        target_id = target_chembl_id.strip().upper()
        url = f"{self.BASE_URL}/{target_id}.json"
        cache_key = f"chembl_{target_id}"
        fixture_name = f"chembl_{target_id}.json"

        data = self.fetch_json(
            url=url,
            cache_key=cache_key,
            fixture_name=fixture_name,
            offline=offline,
        )

        return self._normalize_chembl_response(target_id, data)

    def _normalize_chembl_response(self, target_id: str, raw: Dict[str, Any]) -> Dict[str, Any]:
        """Normalize live API response or embedded fixture to uniform schema."""
        if "bioactivities" in raw and "pref_name" in raw:
            return raw

        target_chembl_id = raw.get("target_chembl_id", target_id)
        pref_name = raw.get("pref_name", "Unknown target")
        target_type = raw.get("target_type", "")
        organism = raw.get("organism", "")
        tax_id = raw.get("tax_id", 0)

        # Extract protein class
        classifications = raw.get("target_components", [])
        protein_class = ""
        if classifications and "target_component_synonyms" in classifications[0]:
            syns = classifications[0].get("target_component_synonyms", [])
            if syns:
                protein_class = syns[0].get("component_synonym", "")

        return {
            "target_chembl_id": target_chembl_id,
            "pref_name": pref_name,
            "target_type": target_type,
            "organism": organism,
            "tax_id": tax_id,
            "protein_class": protein_class,
            "bioactivities": [],
        }


def fetch_chembl_bioactivity(
    target_chembl_id: str,
    offline: bool = False,
    cache_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Convenience function to fetch ChEMBL target data and bioactivities.

    Args:
        target_chembl_id: ChEMBL target identifier (e.g. 'CHEMBL1784').
        offline: If True, operates offline from local cache/fixtures.
        cache_dir: Optional local cache directory path.

    Returns:
        Dictionary of target metadata and bioactivity values.
    """
    connector = ChEMBLConnector(cache_dir=cache_dir)
    return connector.fetch(target_chembl_id, offline=offline)
