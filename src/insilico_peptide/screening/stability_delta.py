"""
Empirical thermodynamic stability (ddG) and aqueous solubility (dSol) heuristics,
along with Pareto multi-objective optimization for candidate ranking.
"""

from __future__ import annotations

import math
from typing import Any, List, Optional, Sequence, Tuple

from insilico_peptide.core.scales import (
    HELIX_PROPENSITY,
    KYTE_DOOLITTLE,
)

# Gatekeeper residues that actively disrupt amyloid beta-sheet aggregation
GATEKEEPER_RESIDUES = frozenset({"P", "K", "R", "E", "D"})


def predict_ddg(
    wt_residue: str,
    mut_residue: str,
    position: int,
    sequence_length: int,
    wt_is_helical: bool = True,
) -> float:
    """
    Predict empirical folding free energy change (delta delta G) in kcal/mol.

    Convention:
        ddG < 0: Stabilizing mutation (thermodynamically favorable)
        ddG > 0: Destabilizing mutation (unfavorable)

    Components:
        1. Helix propensity delta: based on Chou-Fasman / Pace relative propensities.
        2. Hydrophobic core penalty: based on Kyte-Doolittle hydropathy shift.
        3. Steric / charge boundary penalty.
    """
    if wt_residue == mut_residue:
        return 0.0

    p_wt = HELIX_PROPENSITY.get(wt_residue, 0.5)
    p_mut = HELIX_PROPENSITY.get(mut_residue, 0.5)

    # Free energy contribution from secondary structure propensity
    # RT at 298K = 0.592 kcal/mol
    rt = 0.592
    ddg_helix = -rt * math.log(max(1e-3, p_mut) / max(1e-3, p_wt))

    # Hydrophobic transfer delta: losing hydrophobic bulk destabilizes core
    h_wt = KYTE_DOOLITTLE.get(wt_residue, 0.0)
    h_mut = KYTE_DOOLITTLE.get(mut_residue, 0.0)

    # Mutating hydrophobic to hydrophilic in interior increases energy
    ddg_hydrophobic = 0.25 * (h_wt - h_mut)

    # Proline disruption penalty inside a helix
    ddg_proline = 0.0
    if mut_residue == "P" and 1 < position < sequence_length:
        ddg_proline = 1.5
    elif wt_residue == "P" and 1 < position < sequence_length:
        ddg_proline = -0.5

    ddg_total = ddg_helix + ddg_hydrophobic + ddg_proline
    return round(ddg_total, 3)


def predict_solubility_delta(
    wt_sequence: str,
    mut_sequence: str,
    wt_charge: float,
    mut_charge: float,
    wt_gravy: float,
    mut_gravy: float,
) -> float:
    """
    Calculate aqueous solubility change heuristic (delta Sol).

    Positive delta Sol indicates enhanced solubility; negative indicates aggregation propensity.

    Heuristic features:
        1. Net charge magnitude: higher absolute charge creates electrostatic repulsion.
        2. GRAVY: lower hydropathy decreases hydrophobic aggregation.
        3. Gatekeeper count: introduction of P, K, R, E, D disrupts intermolecular beta-aggregation.
    """
    delta_abs_charge = abs(mut_charge) - abs(wt_charge)
    delta_gravy = mut_gravy - wt_gravy

    wt_gatekeepers = sum(1 for aa in wt_sequence if aa in GATEKEEPER_RESIDUES)
    mut_gatekeepers = sum(1 for aa in mut_sequence if aa in GATEKEEPER_RESIDUES)
    delta_gatekeepers = mut_gatekeepers - wt_gatekeepers

    # Empirical weighting factors
    beta = 0.45   # electrostatic repulsion weight
    gamma = 0.35  # hydropathy penalty weight
    delta = 0.30  # structural gatekeeper weight

    sol_score = (beta * delta_abs_charge) - (gamma * delta_gravy) + (delta * delta_gatekeepers)
    return round(sol_score, 3)


def pareto_frontier(
    items: Sequence[Any],
    objectives: List[str],
    maximize: List[bool],
) -> List[int]:
    """
    Compute 0-based indices of items on the non-dominated Pareto frontier.

    Args:
        items: List of objects or dicts.
        objectives: Attribute names or dict keys to evaluate.
        maximize: List of booleans indicating whether each objective is maximized.

    Returns:
        List of indices corresponding to non-dominated items.
    """
    n = len(items)
    if n == 0:
        return []

    # Extract score matrix
    scores: List[List[float]] = []
    for item in items:
        row = []
        for obj in objectives:
            val = getattr(item, obj) if hasattr(item, obj) else item[obj]
            row.append(float(val))
        scores.append(row)

    dominated = [False] * n

    for i in range(n):
        if dominated[i]:
            continue
        for j in range(n):
            if i == j or dominated[i]:
                continue
            # Check if j dominates i
            # j is at least as good in all objectives, and strictly better in at least one
            better_or_equal = True
            strictly_better = False

            for k in range(len(objectives)):
                val_i = scores[i][k]
                val_j = scores[j][k]
                max_k = maximize[k]

                if max_k:
                    if val_j < val_i:
                        better_or_equal = False
                        break
                    if val_j > val_i:
                        strictly_better = True
                else:
                    if val_j > val_i:
                        better_or_equal = False
                        break
                    if val_j < val_i:
                        strictly_better = True

            if better_or_equal and strictly_better:
                dominated[i] = True
                break

    frontier = [i for i in range(n) if not dominated[i]]
    return frontier


def rank_pareto(
    items: List[Any],
    objectives: Optional[List[str]] = None,
    maximize: Optional[List[bool]] = None,
) -> List[Any]:
    """
    Assign Pareto ranks (1 = non-dominated, 2 = second front, etc.) to all items.

    Args:
        items: List of candidate objects (e.g. MutantProfile).
        objectives: Objective fields.
        maximize: Maximization directions.

    Returns:
        New list of items sorted by pareto_rank ascending, then primary objective.
    """
    if objectives is None:
        objectives = ["ddg_pred", "solubility_pred"]
    if maximize is None:
        maximize = [False, True]
    if not items:
        return []

    remaining_indices = list(range(len(items)))
    rank = 1

    while remaining_indices:
        subset = [items[idx] for idx in remaining_indices]
        frontier_sub = pareto_frontier(subset, objectives, maximize)
        if not frontier_sub:
            break

        frontier_indices = [remaining_indices[i] for i in frontier_sub]

        for idx in frontier_indices:
            item = items[idx]
            if hasattr(item, "pareto_rank"):
                item.pareto_rank = rank
            elif isinstance(item, dict):
                item["pareto_rank"] = rank

        remaining_indices = [idx for idx in remaining_indices if idx not in frontier_indices]
        rank += 1

    # Sort items by pareto_rank, then primary objective
    def sort_key(it: Any) -> Tuple[int, float]:
        r = getattr(it, "pareto_rank", 1) if hasattr(it, "pareto_rank") else it.get("pareto_rank", 1)
        val = getattr(it, objectives[0], 0.0) if hasattr(it, objectives[0]) else it.get(objectives[0], 0.0)
        mult = -1.0 if maximize[0] else 1.0
        return (r, mult * float(val))

    return sorted(items, key=sort_key)
