"""
Systematic in-silico Alanine Scanning and structural hotspot identification.

Mutates each position in a peptide sequence to Alanine (or Glycine if the position
is already Alanine) and computes differential biophysical vectors:
- delta MW
- delta pI
- delta GRAVY
- delta Hydrophobic Moment (mu_H)
- delta Instability Index (II)
- delta Net Charge (pH 7.4)
- Empirical ddG and dSol heuristics
- Hotspot classification
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List

from insilico_peptide.core.biophysics import (
    calculate_average_mw,
    calculate_gravy,
    calculate_hydrophobic_moment,
    calculate_instability_index,
    calculate_isoelectric_point,
    calculate_net_charge,
    validate_sequence,
)
from insilico_peptide.screening.stability_delta import (
    predict_ddg,
    predict_solubility_delta,
    rank_pareto,
)


@dataclass
class MutantProfile:
    """Differential biophysical profile of a single amino acid variant."""

    position: int
    wt_residue: str
    mut_residue: str
    mutant_sequence: str
    delta_mw: float
    delta_pi: float
    delta_gravy: float
    delta_moment: float
    delta_instability: float
    delta_charge: float
    ddg_pred: float
    solubility_pred: float
    is_hotspot: bool
    hotspot_reasons: List[str] = field(default_factory=list)
    pareto_rank: int = 1

    def to_dict(self) -> Dict[str, Any]:
        """Convert mutant profile to dictionary."""
        return {
            "position": self.position,
            "mutation": f"{self.wt_residue}{self.position}{self.mut_residue}",
            "wt_residue": self.wt_residue,
            "mut_residue": self.mut_residue,
            "mutant_sequence": self.mutant_sequence,
            "delta_mw": self.delta_mw,
            "delta_pi": self.delta_pi,
            "delta_gravy": self.delta_gravy,
            "delta_moment": self.delta_moment,
            "delta_instability": self.delta_instability,
            "delta_charge": self.delta_charge,
            "ddg_pred": self.ddg_pred,
            "solubility_pred": self.solubility_pred,
            "is_hotspot": self.is_hotspot,
            "hotspot_reasons": self.hotspot_reasons,
            "pareto_rank": self.pareto_rank,
        }


def run_alanine_scan(
    sequence: str,
    pka_scale: str = "EMBOSS",
) -> List[MutantProfile]:
    """
    Perform a complete, systematic in-silico Alanine Scanning on a peptide sequence.

    Args:
        sequence: Wild-type amino acid sequence.
        pka_scale: pKa scale for charge and pI calculations.

    Returns:
        List of MutantProfile objects for all residues in sequence order.
    """
    seq = validate_sequence(sequence)
    n = len(seq)

    # Compute baseline wild-type properties
    wt_mw = calculate_average_mw(seq)
    wt_pi = calculate_isoelectric_point(seq, pka_scale=pka_scale)
    wt_gravy = calculate_gravy(seq)
    wt_moment = calculate_hydrophobic_moment(seq, angle_deg=100.0)
    wt_instability = calculate_instability_index(seq)
    wt_charge = calculate_net_charge(seq, 7.4, pka_scale=pka_scale)

    mutants: List[MutantProfile] = []

    for i in range(1, n + 1):
        wt_aa = seq[i - 1]
        # Mutate to Alanine, or Glycine if already Alanine
        mut_aa = "G" if wt_aa == "A" else "A"
        mut_seq = seq[: i - 1] + mut_aa + seq[i:]

        # Compute mutant properties
        mut_mw = calculate_average_mw(mut_seq)
        mut_pi = calculate_isoelectric_point(mut_seq, pka_scale=pka_scale)
        mut_gravy = calculate_gravy(mut_seq)
        mut_moment = calculate_hydrophobic_moment(mut_seq, angle_deg=100.0)
        mut_instability = calculate_instability_index(mut_seq)
        mut_charge = calculate_net_charge(mut_seq, 7.4, pka_scale=pka_scale)

        # Differential delta vectors
        d_mw = round(mut_mw - wt_mw, 2)
        d_pi = round(mut_pi - wt_pi, 2)
        d_gravy = round(mut_gravy - wt_gravy, 4)
        d_moment = round(mut_moment - wt_moment, 4)
        d_instability = round(mut_instability - wt_instability, 2)
        d_charge = round(mut_charge - wt_charge, 2)

        ddg = predict_ddg(wt_aa, mut_aa, i, n)
        dsol = predict_solubility_delta(seq, mut_seq, wt_charge, mut_charge, wt_gravy, mut_gravy)

        # Hotspot classification criteria
        reasons: List[str] = []
        if abs(d_moment) > 0.15:
            reasons.append(f"Disruption of amphiphilic moment (|d_mu_H|={abs(d_moment):.4f} > 0.15)")
        if abs(d_charge) >= 1.0:
            reasons.append(f"Alteration of formal net charge (|d_Q|={abs(d_charge):.2f} >= 1.0)")
        if d_instability > 10.0:
            reasons.append(f"Severe destabilization of peptide backbone (d_II={d_instability:.2f} > 10.0)")
        if ddg > 1.5:
            reasons.append(f"Thermodynamic unfolding penalty (ddG={ddg:.2f} > 1.5 kcal/mol)")

        is_hotspot = len(reasons) > 0

        profile = MutantProfile(
            position=i,
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
            is_hotspot=is_hotspot,
            hotspot_reasons=reasons,
        )
        mutants.append(profile)

    # Compute Pareto ranking across the mutant candidates
    ranked = rank_pareto(mutants, objectives=["ddg_pred", "solubility_pred"], maximize=[False, True])

    # Return candidates preserved in original sequence position order
    mutants_by_pos = sorted(ranked, key=lambda m: m.position)
    return mutants_by_pos
