"""
UniProt Knowledgebase connector with 3-tier fallback.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from insilico_peptide.connectors.base import (
    BaseConnector,
)


class UniProtConnector(BaseConnector):
    """Client for UniProt KB REST API with offline caching and curated fixtures."""

    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    def fetch(self, accession: str, offline: bool = False) -> Dict[str, Any]:
        """
        Fetch UniProt record by accession (e.g. 'P01275', 'P01308').

        Args:
            accession: UniProt accession ID or recognized benchmark identifier.
            offline: Force offline mode.

        Returns:
            Normalized dictionary containing sequence and annotation metadata.
        """
        acc = accession.strip().upper()
        url = f"{self.BASE_URL}/{acc}.json"
        cache_key = f"uniprot_{acc}"
        fixture_name = f"uniprot_{acc}.json"

        data = self.fetch_json(
            url=url,
            cache_key=cache_key,
            fixture_name=fixture_name,
            offline=offline,
        )

        return self._normalize_uniprot_response(acc, data)

    def _normalize_uniprot_response(self, accession: str, raw: Dict[str, Any]) -> Dict[str, Any]:
        """Normalize raw UniProt response or embedded fixture to uniform schema."""
        # Check if already normalized (fixture schema)
        if "sequence" in raw and isinstance(raw["sequence"], str) and "protein_name" in raw:
            return raw

        # Parse live UniProt KB JSON format
        primary_acc = raw.get("primaryAccession", accession)
        entry_name = raw.get("uniProtkbId", "")

        # Extract protein name
        desc = raw.get("proteinDescription", {})
        rec_name = desc.get("recommendedName", {}).get("fullName", {}).get("value")
        if not rec_name and "submissionNames" in desc and desc["submissionNames"]:
            rec_name = desc["submissionNames"][0].get("fullName", {}).get("value")
        protein_name = rec_name or "Unknown protein"

        # Extract gene name
        genes = raw.get("genes", [])
        gene_name = genes[0].get("geneName", {}).get("value", "") if genes else ""

        # Extract organism
        organism = raw.get("organism", {}).get("scientificName", "")
        tax_id = raw.get("organism", {}).get("taxonId", 0)

        # Extract sequence
        seq_obj = raw.get("sequence", {})
        sequence = seq_obj.get("value", "")
        length = seq_obj.get("length", len(sequence))

        # Extract features
        features: List[Dict[str, Any]] = []
        for feat in raw.get("features", []):
            feat_type = feat.get("type", "")
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            description = feat.get("description", "")
            features.append({
                "type": feat_type,
                "start": start,
                "end": end,
                "description": description,
            })

        return {
            "accession": primary_acc,
            "entry_name": entry_name,
            "protein_name": protein_name,
            "gene_name": gene_name,
            "organism": organism,
            "taxonomy_id": tax_id,
            "sequence": sequence,
            "length": length,
            "features": features,
        }


def fetch_uniprot_data(
    accession: str,
    offline: bool = False,
    cache_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Convenience function to fetch UniProt data.

    Args:
        accession: UniProt accession ID (e.g. 'P01308').
        offline: If True, operates offline from local cache/fixtures.
        cache_dir: Optional local cache directory path.

    Returns:
        Dictionary of sequence and metadata.
    """
    connector = UniProtConnector(cache_dir=cache_dir)
    return connector.fetch(accession, offline=offline)
