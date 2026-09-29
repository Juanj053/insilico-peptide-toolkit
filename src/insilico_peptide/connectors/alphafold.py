"""
AlphaFold Protein Structure Database connector with 3-tier fallback.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from insilico_peptide.connectors.base import BaseConnector


class AlphaFoldConnector(BaseConnector):
    """Client for EBI AlphaFold Structure Database API with offline fallback."""

    BASE_URL = "https://alphafold.ebi.ac.uk/api/prediction"

    def fetch(self, accession: str, offline: bool = False) -> Dict[str, Any]:
        """
        Fetch AlphaFold structure prediction metadata and pLDDT metrics.

        Args:
            accession: UniProt accession ID (e.g. 'P01275', 'P01308').
            offline: Force offline mode.

        Returns:
            Normalized dictionary containing pLDDT scores, URLs, and confidence statistics.
        """
        acc = accession.strip().upper()
        url = f"{self.BASE_URL}/{acc}"
        cache_key = f"alphafold_{acc}"
        fixture_name = f"alphafold_{acc}.json"

        data = self.fetch_json(
            url=url,
            cache_key=cache_key,
            fixture_name=fixture_name,
            offline=offline,
        )

        return self._normalize_alphafold_response(acc, data)

    def _normalize_alphafold_response(self, accession: str, raw: Any) -> Dict[str, Any]:
        """Normalize live API response array or fixture dict."""
        # Check if already normalized fixture dict
        if isinstance(raw, dict) and "mean_plddt" in raw and "pdb_url" in raw:
            return raw

        # Live API returns a list of prediction objects
        entry: Dict[str, Any] = {}
        if isinstance(raw, list) and len(raw) > 0:
            entry = raw[0]
        elif isinstance(raw, dict):
            entry = raw

        entry_id = entry.get("entryId", f"AF-{accession}-F1")
        uniprot_acc = entry.get("uniprotAccession", accession)
        pdb_url = entry.get("pdbUrl", "")
        cif_url = entry.get("cifUrl", "")
        bcif_url = entry.get("bcifUrl", "")
        mean_plddt = float(entry.get("globalMetricValue", 0.0))

        confidence_dist = entry.get("confidenceDistribution", {})
        if not confidence_dist and "mean_plddt" in entry:
            confidence_dist = entry.get("confidence_distribution", {})

        return {
            "entry_id": entry_id,
            "uniprot_accession": uniprot_acc,
            "pdb_url": pdb_url,
            "cif_url": cif_url,
            "bcif_url": bcif_url,
            "mean_plddt": mean_plddt,
            "confidence_distribution": confidence_dist,
            "sequence_length": entry.get("uniprotEnd", 0) - entry.get("uniprotStart", 1) + 1,
        }


def fetch_alphafold_plddt(
    accession: str,
    offline: bool = False,
    cache_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Convenience function to fetch AlphaFold pLDDT prediction data.

    Args:
        accession: UniProt accession ID (e.g. 'P01275').
        offline: If True, operates offline from local cache/fixtures.
        cache_dir: Optional local cache directory path.

    Returns:
        Dictionary of AlphaFold structure prediction metadata.
    """
    connector = AlphaFoldConnector(cache_dir=cache_dir)
    return connector.fetch(accession, offline=offline)
