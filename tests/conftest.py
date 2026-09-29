"""Pytest fixtures and shared reference standards for insilico-peptide-toolkit test suite."""

import pytest


@pytest.fixture
def ghk_sequence() -> str:
    """Human GHK tripeptide."""
    return "GHK"


@pytest.fixture
def bpc157_sequence() -> str:
    """Human gastric juice pentadecapeptide BPC-157."""
    return "GEPPPGKPADDAGLV"


@pytest.fixture
def glp1_canonical() -> str:
    """Human active GLP-1(7-37) sequence."""
    return "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG"


@pytest.fixture
def glp1_amidated_sequence() -> str:
    """Human circulating GLP-1(7-36) amide precursor sequence."""
    return "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR"


@pytest.fixture
def insulin_chain_a() -> str:
    """Human Insulin Chain A (21 AA)."""
    return "GIVEQCCTSICSLYQLENYCN"


@pytest.fixture
def insulin_chain_b() -> str:
    """Human Insulin Chain B (30 AA)."""
    return "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"


@pytest.fixture
def amphipathic_helix_model() -> str:
    """Model amphipathic alpha-helical peptide."""
    return "KLAKLAKKLAKLAK"
