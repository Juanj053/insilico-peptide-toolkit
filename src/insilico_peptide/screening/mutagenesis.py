"""
Conservative and homologous mutational screening using BLOSUM62 and physicochemical clustering.
"""

from __future__ import annotations

from typing import List, Optional, Set

from insilico_peptide.core.biophysics import (
    calculate_average_mw,
    calculate_gravy,
    calculate_hydrophobic_moment,
    calculate_instability_index,
    calculate_isoelectric_point,
    calculate_net_charge,
    validate_sequence,
)
from insilico_peptide.core.scales import (
    BLOSUM62,
    CONSERVATIVE_CLUSTERS,
    STANDARD_AMINO_ACIDS,
)
from insilico_peptide.screening.alanine_scan import MutantProfile
from insilico_peptide.screening.stability_delta import (
    predict_ddg,
    predict_solubility_delta,
    rank_pareto,
)


def get_conservative_substitutions(residue: str, min_blosum62: int = 0) -> List[str]:
    """
    Find conservative amino acid substitutions for a given residue based on
    physicochemical clustering and positive or neutral BLOSUM62 substitution scores.

    Args:
        residue: Single-letter amino acid code.
        min_blosum62: Minimum BLOSUM62 substitution score threshold.

    Returns:
        Sorted list of candidate substitute residues (excluding wild-type).
    """
    res_upper = residue.upper()
    if res_upper not in STANDARD_AMINO_ACIDS:
        raise ValueError(f"Unknown amino acid residue '{residue}'")

    candidates: Set[str] = set()

    for cluster in CONSERVATIVE_CLUSTERS.values():
        if res_upper in cluster:
            for alt in cluster:
                if alt != res_upper:
                    score = BLOSUM62.get(res_upper, {}).get(alt, -4)
                    if score >= min_blosum62:
                        candidates.add(alt)

    # Fallback to any positive BLOSUM62 substitution if clustering didn't catch any
    if not candidates:
        for alt, score in BLOSUM62.get(res_upper, {}).items():
            if alt != res_upper and alt in STANDARD_AMINO_ACIDS and score >= min_blosum62:
                candidates.add(alt)

    return sorted(list(candidates))


def run_site_mutagenesis(
    sequence: str,
    position: int,
    target_residue: str,
    pka_scale: str = "EMBOSS",
) -> MutantProfile:
    """
    Evaluate a specific single-point mutation at position (1-based).

    Args:
        sequence: Wild-type sequence.
        position: 1-based residue position.
        target_residue: Target mutated amino acid.
        pka_scale: pKa scale name.

    Returns:
        MutantProfile with computed biophysical deltas.
    """
    seq = validate_sequence(sequence)
    mut_aa = target_residue.upper()
    if mut_aa not in STANDARD_AMINO_ACIDS:
        raise ValueError(f"Invalid target amino acid '{target_residue}'")

    if position < 1 or position > len(seq):
        raise IndexError(f"Position {position} is out of bounds for sequence of length {len(seq)}")

    wt_aa = seq[position - 1]
    mut_seq = seq[: position - 1] + mut_aa + seq[position:]

    # Compute baseline properties
    wt_mw = calculate_average_mw(seq)
    wt_pi = calculate_isoelectric_point(seq, pka_scale=pka_scale)
    wt_gravy = calculate_gravy(seq)
    wt_moment = calculate_hydrophobic_moment(seq, angle_deg=100.0)
    wt_instability = calculate_instability_index(seq)
    wt_charge = calculate_net_charge(seq, 7.4, pka_scale=pka_scale)

    # Compute mutant properties
    mut_mw = calculate_average_mw(mut_seq)
    mut_pi = calculate_isoelectric_point(mut_seq, pka_scale=pka_scale)
    mut_gravy = calculate_gravy(mut_seq)
    mut_moment = calculate_hydrophobic_moment(mut_seq, angle_deg=100.0)
    mut_instability = calculate_instability_index(mut_seq)
    mut_charge = calculate_net_charge(mut_seq, 7.4, pka_scale=pka_scale)

    d_mw = round(mut_mw - wt_mw, 2)
    d_pi = round(mut_pi - wt_pi, 2)
    d_gravy = round(mut_gravy - wt_gravy, 4)
    d_moment = round(mut_moment - wt_moment, 4)
    d_instability = round(mut_instability - wt_instability, 2)
    d_charge = round(mut_charge - wt_charge, 2)

    ddg = predict_ddg(wt_aa, mut_aa, position, len(seq))
    dsol = predict_solubility_delta(seq, mut_seq, wt_charge, mut_charge, wt_gravy, mut_gravy)

    reasons: List[str] = []
    if abs(d_moment) > 0.15:
        reasons.append(f"Disruption of amphiphilic moment (|d_mu_H|={abs(d_moment):.4f} > 0.15)")
    if abs(d_charge) >= 1.0:
        reasons.append(f"Alteration of formal net charge (|d_Q|={abs(d_charge):.2f} >= 1.0)")
    if d_instability > 10.0:
        reasons.append(f"Destabilization of backbone (d_II={d_instability:.2f} > 10.0)")
    if ddg > 1.5:
        reasons.append(f"Thermodynamic penalty (ddG={ddg:.2f} > 1.5 kcal/mol)")

    return MutantProfile(
        position=position,
        wt_residue=wt_aa,
        mut_residue=mut_aa,
        mutant_sequence=mut_seq,
        delta_mw=d_mw,
        delta_pi=d_pi,
        delta_gravy=d_gravy,
        delta_moment=d_moment,
        delta_instability=d_instability,
        delta_charge=d_charge,
        ddg_pred=ddg,
        solubility_pred=dsol,
        is_hotspot=len(reasons) > 0,
        hotspot_reasons=reasons,
    )


def run_homologous_scan(
    sequence: str,
    positions: Optional[List[int]] = None,
    pka_scale: str = "EMBOSS",
    min_blosum62: int = 0,
) -> List[MutantProfile]:
    """
    Perform homologous mutational scan across designated or all positions.

    Args:
        sequence: Wild-type sequence.
        positions: Optional list of 1-based positions to scan. If None, scans all positions.
        pka_scale: pKa scale name.
        min_blosum62: Minimum substitution score threshold.

    Returns:
        List of MutantProfile candidates sorted by Pareto rank.
    """
    seq = validate_sequence(sequence)
    target_positions = positions if positions is not None else list(range(1, len(seq) + 1))

    candidates: List[MutantProfile] = []

    for pos in target_positions:
        if pos < 1 or pos > len(seq):
            continue
        wt_aa = seq[pos - 1]
        subs = get_conservative_substitutions(wt_aa, min_blosum62=min_blosum62)

        for mut_aa in subs:
            profile = run_site_mutagenesis(seq, pos, mut_aa, pka_scale=pka_scale)
            candidates.append(profile)

    return rank_pareto(candidates, objectives=["ddg_pred", "solubility_pred"], maximize=[False, True])
