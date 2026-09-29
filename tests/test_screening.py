"""
Unit and integration tests for in-silico screening, alanine scanning, and Pareto ranking.
"""

import pytest

from insilico_peptide.screening.alanine_scan import (
    MutantProfile,
    run_alanine_scan,
)
from insilico_peptide.screening.mutagenesis import (
    get_conservative_substitutions,
    run_homologous_scan,
    run_site_mutagenesis,
)
from insilico_peptide.screening.stability_delta import (
    predict_ddg,
    predict_solubility_delta,
    rank_pareto,
)


def test_alanine_scanning_deltas(bpc157_sequence: str) -> None:
    """Verify systematic alanine scanning produces delta vectors for all residues."""
    mutants = run_alanine_scan(bpc157_sequence)

    assert len(mutants) == len(bpc157_sequence)

    # Check first residue (G1A)
    g1 = mutants[0]
    assert g1.position == 1
    assert g1.wt_residue == "G"
    assert g1.mut_residue == "A"
    assert g1.delta_mw > 0  # Ala is heavier than Gly (+14.03 Da)
    assert pytest.approx(g1.delta_gravy, abs=0.01) == round((1.8 - (-0.4)) / 15, 4)

    # Residue already Ala (e.g. A9) mutates to Gly
    a9 = mutants[8]
    assert a9.position == 9
    assert a9.wt_residue == "A"
    assert a9.mut_residue == "G"
    assert a9.delta_mw < 0  # Gly is lighter than Ala

    # Asp10 mutation (D10A) changes net charge and hydropathy
    d10 = mutants[9]
    assert d10.position == 10
    assert d10.wt_residue == "D"
    assert d10.mut_residue == "A"
    assert d10.delta_pi > 0  # Loss of acidic residue raises pI
    assert d10.delta_charge > 0  # Loss of negative charge makes net charge more positive
    assert d10.delta_gravy > 0  # Ala is more hydrophobic than Asp (1.8 - (-3.5) = +5.3 / 15)


def test_alanine_scanning_hotspots(glp1_canonical: str) -> None:
    """Verify hotspot detection triggers for residues with major structural or charge disruption."""
    mutants = run_alanine_scan(glp1_canonical)
    hotspots = [m for m in mutants if m.is_hotspot]

    assert len(hotspots) > 0

    # Charged residues (K, R, D, E) should be flagged as hotspots due to |d_Q| >= 1.0
    charged_positions = {i for i, aa in enumerate(glp1_canonical, 1) if aa in ("K", "R", "D", "E")}
    hotspot_positions = {m.position for m in hotspots}
    assert charged_positions.issubset(hotspot_positions)


def test_homologous_substitutions_conservative() -> None:
    """Verify conservative substitutions map within physicochemical clusters."""
    subs_k = get_conservative_substitutions("K")
    assert "R" in subs_k or "H" in subs_k

    subs_d = get_conservative_substitutions("D")
    assert "E" in subs_d

    subs_l = get_conservative_substitutions("L")
    assert "I" in subs_l or "V" in subs_l or "M" in subs_l

    # Unknown residue raises ValueError
    with pytest.raises(ValueError):
        get_conservative_substitutions("X")


def test_site_mutagenesis_single_point(glp1_canonical: str) -> None:
    """Verify site-directed mutagenesis at specific position."""
    # Mutate Ala2 to Gly2 in GLP-1
    mut = run_site_mutagenesis(glp1_canonical, position=2, target_residue="G")
    assert mut.position == 2
    assert mut.wt_residue == "A"
    assert mut.mut_residue == "G"
    assert mut.mutant_sequence.startswith("HG")
    assert isinstance(mut.to_dict(), dict)

    # Out of bounds position raises IndexError
    with pytest.raises(IndexError):
        run_site_mutagenesis(glp1_canonical, position=99, target_residue="A")


def test_homologous_scan(ghk_sequence: str) -> None:
    """Verify homologous scan across peptide."""
    candidates = run_homologous_scan(ghk_sequence)
    assert len(candidates) > 0
    for cand in candidates:
        assert cand.pareto_rank >= 1


def test_stability_ddg_and_solubility() -> None:
    """Verify ddG thermodynamic and dSol solubility heuristics."""
    # Same residue should have 0.0 ddG
    assert predict_ddg("A", "A", 1, 10) == 0.0

    # Proline introduction inside helix incurs penalty
    ddg_pro = predict_ddg("A", "P", 5, 15)
    assert ddg_pro > 0.0

    # Adding net charge increases solubility score
    dsol = predict_solubility_delta(
        wt_sequence="AA",
        mut_sequence="AK",
        wt_charge=0.0,
        mut_charge=+1.0,
        wt_gravy=1.8,
        mut_gravy=-1.0,
    )
    assert dsol > 0.0


def test_pareto_ranking() -> None:
    """Verify Pareto frontier and ranking multi-objective optimization."""
    # Create sample profiles
    p1 = MutantProfile(1, "A", "G", "G", 0, 0, 0, 0, 0, 0, ddg_pred=-2.0, solubility_pred=1.5, is_hotspot=False)
    p2 = MutantProfile(2, "A", "V", "V", 0, 0, 0, 0, 0, 0, ddg_pred=-1.0, solubility_pred=0.5, is_hotspot=False)
    p3 = MutantProfile(3, "A", "I", "I", 0, 0, 0, 0, 0, 0, ddg_pred=1.0, solubility_pred=-1.0, is_hotspot=False)

    ranked = rank_pareto([p1, p2, p3], objectives=["ddg_pred", "solubility_pred"], maximize=[False, True])
    # p1 strictly dominates p2 and p3 in both ddg (lower) and solubility (higher)
    assert ranked[0].position == 1
    assert ranked[0].pareto_rank == 1
    assert ranked[-1].position == 3
