"""
Unit tests for core biophysical calculations and reference peptide standards:
Human Insulin, GLP-1, GHK-Cu, BPC-157.
"""

import math

import pytest

from insilico_peptide.core.biophysics import (
    calculate_aliphatic_index,
    calculate_average_mw,
    calculate_gravy,
    calculate_hydrophobic_moment,
    calculate_instability_index,
    calculate_isoelectric_point,
    calculate_monoisotopic_mw,
    calculate_net_charge,
    calculate_properties,
    calculate_sliding_hydropathy,
    predict_proteolytic_cleavage,
    titrate_charge,
    validate_sequence,
)


def test_molecular_weight_reference_standards(
    ghk_sequence: str,
    bpc157_sequence: str,
    glp1_canonical: str,
    glp1_amidated_sequence: str,
    insulin_chain_a: str,
    insulin_chain_b: str,
) -> None:
    """Verify exact molecular weight calculation against benchmark standards."""
    # GHK (3 AA)
    mw_ghk_avg = calculate_average_mw(ghk_sequence)
    mw_ghk_mono = calculate_monoisotopic_mw(ghk_sequence)
    assert pytest.approx(mw_ghk_avg, abs=0.15) == 340.38
    assert pytest.approx(mw_ghk_mono, abs=0.05) == 340.1859

    # BPC-157 (15 AA)
    mw_bpc_avg = calculate_average_mw(bpc157_sequence)
    mw_bpc_mono = calculate_monoisotopic_mw(bpc157_sequence)
    assert pytest.approx(mw_bpc_avg, abs=0.15) == 1419.54
    assert pytest.approx(mw_bpc_mono, abs=0.05) == 1418.7042

    # GLP-1(7-37) (31 AA)
    mw_glp1_avg = calculate_average_mw(glp1_canonical)
    mw_glp1_mono = calculate_monoisotopic_mw(glp1_canonical)
    assert pytest.approx(mw_glp1_avg, abs=0.20) == 3355.67
    assert pytest.approx(mw_glp1_mono, abs=0.05) == 3353.6681

    # GLP-1(7-36) (30 AA)
    mw_glp36_avg = calculate_average_mw(glp1_amidated_sequence)
    mw_glp36_mono = calculate_monoisotopic_mw(glp1_amidated_sequence)
    assert pytest.approx(mw_glp36_avg, abs=0.20) == 3298.62
    assert pytest.approx(mw_glp36_mono, abs=0.05) == 3296.6466

    # Insulin Chains
    assert pytest.approx(calculate_average_mw(insulin_chain_a), abs=0.20) == 2383.70
    assert pytest.approx(calculate_average_mw(insulin_chain_b), abs=0.20) == 3429.93


def test_terminal_modifications_and_disulfides(glp1_amidated_sequence: str) -> None:
    """Verify N-terminal acetylation, C-terminal amidation, and disulfide bonds."""
    # Free acid vs C-terminal amide
    free_mw = calculate_average_mw(glp1_amidated_sequence, c_terminal="free")
    amide_mw = calculate_average_mw(glp1_amidated_sequence, c_terminal="amide")
    assert amide_mw < free_mw
    assert pytest.approx(free_mw - amide_mw, abs=0.05) == 0.98

    # N-terminal acetylation
    acetyl_mw = calculate_average_mw("GHK", n_terminal="acetyl")
    free_ghk_mw = calculate_average_mw("GHK", n_terminal="free")
    assert pytest.approx(acetyl_mw - free_ghk_mw, abs=0.05) == 42.04

    # Disulfide bond mass loss (loss of 2H per bond)
    oxidized_mw = calculate_average_mw("CGC", num_disulfide_bonds=1)
    reduced_mw = calculate_average_mw("CGC", num_disulfide_bonds=0)
    assert pytest.approx(reduced_mw - oxidized_mw, abs=0.05) == 2.016


def test_isoelectric_point_and_scales(
    ghk_sequence: str,
    bpc157_sequence: str,
    glp1_canonical: str,
) -> None:
    """Verify pI bisection solver across EMBOSS, Lehninger, and Bjellqvist scales."""
    # GHK: basic tripeptide
    pi_ghk_emboss = calculate_isoelectric_point(ghk_sequence, pka_scale="EMBOSS")
    pi_ghk_lehninger = calculate_isoelectric_point(ghk_sequence, pka_scale="LEHNINGER")
    assert pytest.approx(pi_ghk_emboss, abs=0.10) == 9.70
    assert pytest.approx(pi_ghk_lehninger, abs=0.10) == 10.11

    # BPC-157: acidic pentadecapeptide
    pi_bpc_emboss = calculate_isoelectric_point(bpc157_sequence, pka_scale="EMBOSS")
    pi_bpc_lehninger = calculate_isoelectric_point(bpc157_sequence, pka_scale="LEHNINGER")
    pi_bpc_bjellqvist = calculate_isoelectric_point(bpc157_sequence, pka_scale="BJELLQVIST")
    assert pytest.approx(pi_bpc_emboss, abs=0.05) == 3.88
    assert pytest.approx(pi_bpc_lehninger, abs=0.05) == 3.70
    assert pytest.approx(pi_bpc_bjellqvist, abs=0.05) == 4.05

    # GLP-1(7-37): pI ~ 5.4 - 5.6
    pi_glp_emboss = calculate_isoelectric_point(glp1_canonical, pka_scale="EMBOSS")
    pi_glp_lehninger = calculate_isoelectric_point(glp1_canonical, pka_scale="LEHNINGER")
    assert pytest.approx(pi_glp_emboss, abs=0.05) == 5.62
    assert pytest.approx(pi_glp_lehninger, abs=0.05) == 5.43


def test_net_charge_monotonicity_and_neutral_ph(
    bpc157_sequence: str,
    glp1_canonical: str,
) -> None:
    """Verify thermodynamic monotonicity of titration curve and physiological charge."""
    titration = titrate_charge(bpc157_sequence, pka_scale="EMBOSS", ph_min=0.0, ph_max=14.0, step=0.5)

    phs = sorted(titration.keys())
    for i in range(len(phs) - 1):
        ph1, ph2 = phs[i], phs[i + 1]
        assert titration[ph1] >= titration[ph2], (
            f"Titration monotonicity violation: Q({ph1})={titration[ph1]} < Q({ph2})={titration[ph2]}"
        )

    # Check physiological charges
    q74_bpc = calculate_net_charge(bpc157_sequence, 7.4, pka_scale="EMBOSS")
    q74_glp1 = calculate_net_charge(glp1_canonical, 7.4, pka_scale="EMBOSS")
    assert pytest.approx(q74_bpc, abs=0.10) == -2.06
    assert pytest.approx(q74_glp1, abs=0.10) == -0.95


def test_gravy_calculation_exactness(
    ghk_sequence: str,
    bpc157_sequence: str,
    glp1_canonical: str,
    insulin_chain_a: str,
) -> None:
    """Verify Kyte-Doolittle GRAVY hydropathicity exactness."""
    assert pytest.approx(calculate_gravy(ghk_sequence), abs=0.0005) == -2.5000
    assert pytest.approx(calculate_gravy(bpc157_sequence), abs=0.0005) == -0.6933
    assert pytest.approx(calculate_gravy(glp1_canonical), abs=0.0005) == -0.2355
    assert pytest.approx(calculate_gravy(insulin_chain_a), abs=0.0005) == 0.2143

    # Sliding window length check
    window_profile = calculate_sliding_hydropathy(bpc157_sequence, window_size=5)
    assert len(window_profile) == len(bpc157_sequence)


def test_instability_index_and_classification(
    bpc157_sequence: str,
    glp1_canonical: str,
    ghk_sequence: str,
    insulin_chain_b: str,
) -> None:
    """Verify Guruprasad et al. DIWV matrix instability index calculation and stability cutoff."""
    # BPC-157: Stable (39.32 <= 40.0)
    ii_bpc = calculate_instability_index(bpc157_sequence)
    assert pytest.approx(ii_bpc, abs=0.10) == 39.32
    assert ii_bpc <= 40.0

    # GLP-1(7-37): Highly stable (14.71 <= 40.0)
    ii_glp = calculate_instability_index(glp1_canonical)
    assert pytest.approx(ii_glp, abs=0.10) == 14.71
    assert ii_glp <= 40.0

    # GHK: Flexible tripeptide (85.60 > 40.0)
    ii_ghk = calculate_instability_index(ghk_sequence)
    assert pytest.approx(ii_ghk, abs=0.10) == 85.60
    assert ii_ghk > 40.0

    # Insulin B: Stable (9.85 <= 40.0)
    ii_ins_b = calculate_instability_index(insulin_chain_b)
    assert pytest.approx(ii_ins_b, abs=0.10) == 9.85
    assert ii_ins_b <= 40.0


def test_amphiphilic_hydrophobic_moment(
    amphipathic_helix_model: str,
    glp1_canonical: str,
) -> None:
    """Verify Eisenberg hydrophobic moment and amphipathic model cutoff."""
    mu_model = calculate_hydrophobic_moment(amphipathic_helix_model, angle_deg=100.0)
    assert mu_model > 0.35, f"Amphipathic model moment should exceed 0.35, got {mu_model}"

    # Homopolymer over complete helical turns (18 residues = 5 full 360-degree turns) cancels to zero
    mu_homo = calculate_hydrophobic_moment("A" * 18, angle_deg=100.0)
    assert mu_homo < 1e-4

    # Conformation string shortcuts
    mu_alpha = calculate_hydrophobic_moment(glp1_canonical, conformation="alpha")
    mu_310 = calculate_hydrophobic_moment(glp1_canonical, conformation="3_10")
    mu_beta = calculate_hydrophobic_moment(glp1_canonical, conformation="beta")
    assert mu_alpha > 0
    assert mu_310 > 0
    assert mu_beta > 0


def test_aliphatic_index_calculation(
    ghk_sequence: str,
    bpc157_sequence: str,
    insulin_chain_a: str,
    glp1_canonical: str,
) -> None:
    """Verify Ikai Aliphatic Index calculation against exact benchmark mole fractions."""
    assert calculate_aliphatic_index(ghk_sequence) == 0.0
    assert pytest.approx(calculate_aliphatic_index(bpc157_sequence), abs=0.05) == 58.67
    assert pytest.approx(calculate_aliphatic_index(insulin_chain_a), abs=0.05) == 88.10
    assert pytest.approx(calculate_aliphatic_index(glp1_canonical), abs=0.05) == 69.35


def test_proteolytic_cleavage_trypsin_proline_block(
    bpc157_sequence: str,
    glp1_canonical: str,
    insulin_chain_b: str,
) -> None:
    """Verify Trypsin cleavage and the physiological Proline-block rule."""
    # BPC-157: Lys7 is followed by Pro8 -> ZERO Trypsin cuts!
    cleavage_bpc = predict_proteolytic_cleavage(bpc157_sequence, proteases=["trypsin"])
    assert cleavage_bpc["cleavage_sites"]["trypsin"] == []
    assert cleavage_bpc["fragments"]["trypsin"] == [bpc157_sequence]

    # GLP-1(7-37): K20, K28, R30 -> Trypsin cuts at 20, 28, 30
    cleavage_glp1 = predict_proteolytic_cleavage(glp1_canonical, proteases=["trypsin"])
    assert cleavage_glp1["cleavage_sites"]["trypsin"] == [20, 28, 30]

    # Insulin B: R22, K29
    cleavage_ins = predict_proteolytic_cleavage(insulin_chain_b, proteases=["trypsin"])
    assert cleavage_ins["cleavage_sites"]["trypsin"] == [22, 29]


def test_proteolytic_cleavage_dpp4_and_chymotrypsin(glp1_canonical: str) -> None:
    """Verify DPP-IV N-terminal dipeptide cleavage and Chymotrypsin aromatic cleavage."""
    # GLP-1 has Ala at position 2 -> DPP-IV cleaves after residue 2
    cleavage_native = predict_proteolytic_cleavage(glp1_canonical, proteases=["dpp4", "chymotrypsin"])
    assert cleavage_native["cleavage_sites"]["dpp4"] == [2]
    assert cleavage_native["fragments"]["dpp4"][0] == "HA"

    # Mutate Ala2 to Gly2 -> Resistant to DPP-IV
    mutated_glp1 = "HGEGTFTSDVSSYLEGQAAKEFIAWLVKGRG"
    cleavage_mut = predict_proteolytic_cleavage(mutated_glp1, proteases=["dpp4"])
    assert cleavage_mut["cleavage_sites"]["dpp4"] == []

    # Chymotrypsin aromatic cleavage: F6, Y13, F22, W25
    assert cleavage_native["cleavage_sites"]["chymotrypsin"] == [6, 13, 22, 25]


def test_boundary_empty_and_single_residues() -> None:
    """Test boundary conditions for empty sequence and single amino acid."""
    with pytest.raises(ValueError, match="Sequence cannot be empty"):
        validate_sequence("")

    # Single residue Alanine
    mw_a = calculate_average_mw("A")
    assert pytest.approx(mw_a, abs=0.05) == 89.09
    assert calculate_gravy("A") == 1.8
    assert calculate_instability_index("A") == 0.0


def test_boundary_invalid_characters() -> None:
    """Verify sanitization and rejection of invalid amino acid characters."""
    with pytest.raises(ValueError, match="Invalid amino acid residue"):
        validate_sequence("GHK-Cu")

    with pytest.raises(ValueError, match="Invalid amino acid residue"):
        validate_sequence("ACDEFGH123")

    with pytest.raises(ValueError, match="Invalid amino acid residue"):
        validate_sequence("XYZB")


def test_extreme_ph_and_normalization() -> None:
    """Verify sequence whitespace/case normalization and extreme pH stability."""
    norm = validate_sequence("  ghk \n\t")
    assert norm == "GHK"

    # Extreme pH titration calculations should execute without overflow
    q_low = calculate_net_charge("GHK", 0.0)
    q_high = calculate_net_charge("GHK", 14.0)
    assert q_low > 0.0
    assert q_high < 0.0
    assert not math.isnan(q_low)
    assert not math.isnan(q_high)


def test_calculate_properties_facade(bpc157_sequence: str) -> None:
    """Verify full calculate_properties facade dataclass."""
    props = calculate_properties(bpc157_sequence)
    assert props.sequence == bpc157_sequence
    assert props.length == 15
    assert props.is_stable is True
    assert props.molecular_weight > 1000
    assert "trypsin" in props.proteolytic_cleavage_sites
    assert isinstance(props.to_dict(), dict)


def test_core_error_handling() -> None:
    """Verify exceptions for invalid scales, conformations, proteases, and window sizes."""
    # Invalid pKa scale
    with pytest.raises(ValueError, match="Unknown pKa scale"):
        calculate_net_charge("GHK", 7.4, pka_scale="UNKNOWN_SCALE")

    # Invalid conformation
    with pytest.raises(ValueError, match="Unknown conformation"):
        calculate_hydrophobic_moment("GHK", conformation="gamma_loop")

    # Invalid window size
    with pytest.raises(ValueError, match="Window size must be >= 1"):
        calculate_sliding_hydropathy("GHK", window_size=0)

    # Unknown protease
    with pytest.raises(ValueError, match="Unknown protease"):
        predict_proteolytic_cleavage("GHK", proteases=["unknown_enzyme"])

    # Non-string sequence
    with pytest.raises(TypeError, match="Sequence must be a string"):
        validate_sequence(12345)  # type: ignore

