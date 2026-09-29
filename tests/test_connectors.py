"""
Unit and integration tests for multi-omic connectors: UniProt, AlphaFold DB, ChEMBL.
"""

import pytest

from insilico_peptide.connectors.alphafold import (
    AlphaFoldConnector,
    fetch_alphafold_plddt,
)
from insilico_peptide.connectors.base import (
    BaseConnector,
    ConnectorNotFoundError,
)
from insilico_peptide.connectors.chembl import (
    ChEMBLConnector,
    fetch_chembl_bioactivity,
)
from insilico_peptide.connectors.uniprot import (
    UniProtConnector,
    fetch_uniprot_data,
)


def test_uniprot_connector_offline() -> None:
    """Verify UniProt connector loads embedded fixtures in offline mode."""
    # Test GLP-1 / Glucagon preproprotein
    data_glp = fetch_uniprot_data("P01275", offline=True)
    assert data_glp["accession"] == "P01275"
    assert data_glp["protein_name"] == "Glucagon"
    assert "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG" in data_glp["sequence"]

    # Test Insulin
    data_ins = fetch_uniprot_data("P01308", offline=True)
    assert data_ins["accession"] == "P01308"
    assert data_ins["protein_name"] == "Insulin"
    assert "chains" in data_ins or "sequence" in data_ins

    # Test BPC-157 fixture
    data_bpc = fetch_uniprot_data("BPC157", offline=True)
    assert data_bpc["sequence"] == "GEPPPGKPADDAGLV"


def test_alphafold_connector_offline() -> None:
    """Verify AlphaFold connector parses structure metadata and pLDDT scores."""
    af_data = fetch_alphafold_plddt("P01275", offline=True)
    assert af_data["entry_id"] == "AF-P01275-F1"
    assert af_data["mean_plddt"] > 70.0
    assert "pdb_url" in af_data
    assert "confidence_distribution" in af_data


def test_chembl_connector_offline() -> None:
    """Verify ChEMBL connector parses target bioactivities and IC50/Ki values."""
    ch_data = fetch_chembl_bioactivity("CHEMBL1784", offline=True)
    assert ch_data["target_chembl_id"] == "CHEMBL1784"
    assert "Glucagon-like peptide 1 receptor" in ch_data["pref_name"]
    assert len(ch_data["bioactivities"]) > 0

    # Verify bioactivity records
    semaglutide_act = [b for b in ch_data["bioactivities"] if "Semaglutide" in b.get("molecule_name", "")]
    assert len(semaglutide_act) > 0
    assert semaglutide_act[0]["standard_type"] == "EC50"


def test_connector_404_handling() -> None:
    """Verify unknown accession raises ConnectorNotFoundError."""
    connector = UniProtConnector()
    with pytest.raises(ConnectorNotFoundError):
        connector.fetch("NONEXISTENT_ACCESSION_99999", offline=True)


def test_connector_disk_cache(tmp_path) -> None:
    """Verify saving and retrieving data from local disk cache."""
    connector = BaseConnector(cache_dir=str(tmp_path))
    test_key = "test_item_123"
    payload = {"test": "data", "value": 42}

    assert connector.load_from_cache(test_key) is None
    connector.save_to_cache(test_key, payload)
    loaded = connector.load_from_cache(test_key)
    assert loaded == payload


def test_uniprot_normalization_live_format() -> None:
    """Verify UniProt normalization converts raw live REST format."""
    connector = UniProtConnector()
    raw_live = {
        "primaryAccession": "P99999",
        "uniProtkbId": "TEST_HUMAN",
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": "Test Protein"}}
        },
        "genes": [{"geneName": {"value": "TST"}}],
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "sequence": {"value": "ACDEF", "length": 5},
        "features": [
            {
                "type": "Region",
                "location": {"start": {"value": 1}, "end": {"value": 5}},
                "description": "Domain 1",
            }
        ],
    }
    norm = connector._normalize_uniprot_response("P99999", raw_live)
    assert norm["accession"] == "P99999"
    assert norm["protein_name"] == "Test Protein"
    assert norm["sequence"] == "ACDEF"
    assert len(norm["features"]) == 1


def test_alphafold_normalization_live_format() -> None:
    """Verify AlphaFold normalization converts raw list response."""
    connector = AlphaFoldConnector()
    raw_list = [
        {
            "entryId": "AF-P99999-F1",
            "uniprotAccession": "P99999",
            "pdbUrl": "https://example.com/test.pdb",
            "cifUrl": "https://example.com/test.cif",
            "globalMetricValue": 88.5,
            "confidenceDistribution": {"very_high": 80.0},
            "uniprotStart": 1,
            "uniprotEnd": 100,
        }
    ]
    norm = connector._normalize_alphafold_response("P99999", raw_list)
    assert norm["entry_id"] == "AF-P99999-F1"
    assert norm["mean_plddt"] == 88.5
    assert norm["sequence_length"] == 100


def test_chembl_normalization_live_format() -> None:
    """Verify ChEMBL normalization converts live target schema."""
    connector = ChEMBLConnector()
    raw_live = {
        "target_chembl_id": "CHEMBL999",
        "pref_name": "Test Kinase",
        "target_type": "SINGLE PROTEIN",
        "organism": "Homo sapiens",
        "tax_id": 9606,
        "target_components": [
            {
                "target_component_synonyms": [
                    {"component_synonym": "Protein Kinase C"}
                ]
            }
        ],
    }
    norm = connector._normalize_chembl_response("CHEMBL999", raw_live)
    assert norm["target_chembl_id"] == "CHEMBL999"
    assert norm["protein_class"] == "Protein Kinase C"

