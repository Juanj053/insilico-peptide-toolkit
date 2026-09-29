"""
Pure-Python biophysical engine for peptide property modeling.

Implements mathematically rigorous algorithms for:
- Monoisotopic and average molecular weight with terminal modifications and disulfide bonds
- Net charge titration curve and bisection root-finding for isoelectric point (pI)
- Kyte-Doolittle hydropathicity (GRAVY) and sliding-window hydropathy profiles
- Eisenberg consensus hydrophobic moment (mu_H) across helical and extended conformations
- Guruprasad instability index (II) using the 20x20 DIWV matrix
- Ikai aliphatic index (AI)
- Proteolytic cleavage prediction (Trypsin with Proline-block, Chymotrypsin, Pepsin, DPP-IV, Proteinase K, Lys-C, Asp-N)
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from insilico_peptide.core.scales import (
    C_TERM_AMIDE_DELTA_AVG,
    C_TERM_AMIDE_DELTA_MONO,
    DIWV_MATRIX,
    EISENBERG_CONSENSUS,
    KYTE_DOOLITTLE,
    N_TERM_ACETYL_DELTA_AVG,
    N_TERM_ACETYL_DELTA_MONO,
    PKA_SCALES,
    RESIDUE_MASSES_AVERAGE,
    RESIDUE_MASSES_MONOISOTOPIC,
    STANDARD_AMINO_ACIDS,
    WATER_MW_AVERAGE,
    WATER_MW_MONOISOTOPIC,
)


def validate_sequence(sequence: str) -> str:
    """
    Sanitize and validate an amino acid sequence string.

    Args:
        sequence: Raw amino acid sequence string.

    Returns:
        Stripped, uppercase canonical amino acid string.

    Raises:
        ValueError: If sequence is empty or contains non-canonical amino acid codes.
    """
    if not isinstance(sequence, str):
        raise TypeError(f"Sequence must be a string, got {type(sequence).__name__}")

    sanitized = re.sub(r"\s+", "", sequence).upper()

    if not sanitized:
        raise ValueError("Sequence cannot be empty")

    invalid_chars = set(sanitized) - STANDARD_AMINO_ACIDS
    if invalid_chars:
        sorted_invalid = sorted(list(invalid_chars))
        raise ValueError(
            f"Invalid amino acid residue in sequence: {', '.join(sorted_invalid)}. "
            f"Only standard 20 amino acids are supported."
        )

    return sanitized


def calculate_molecular_weight(
    sequence: str,
    monoisotopic: bool = False,
    n_terminal: str = "free",
    c_terminal: str = "free",
    num_disulfide_bonds: int = 0,
) -> float:
    """
    Calculate molecular weight of a peptide in Daltons.

    Args:
        sequence: Amino acid sequence.
        monoisotopic: If True, compute monoisotopic mass; otherwise average molecular weight.
        n_terminal: 'free' or 'acetyl' / 'ac'.
        c_terminal: 'free' or 'amide' / 'nh2'.
        num_disulfide_bonds: Number of formed cystine bridges (-2.016 Da per bond).

    Returns:
        Molecular weight in Daltons.
    """
    seq = validate_sequence(sequence)

    mass_table = RESIDUE_MASSES_MONOISOTOPIC if monoisotopic else RESIDUE_MASSES_AVERAGE
    water_mass = WATER_MW_MONOISOTOPIC if monoisotopic else WATER_MW_AVERAGE

    mass = sum(mass_table[aa] for aa in seq) + water_mass

    # Terminal modifications
    n_term_lower = n_terminal.lower()
    if n_term_lower in ("acetyl", "ac"):
        mass += N_TERM_ACETYL_DELTA_MONO if monoisotopic else N_TERM_ACETYL_DELTA_AVG

    c_term_lower = c_terminal.lower()
    if c_term_lower in ("amide", "nh2"):
        mass += C_TERM_AMIDE_DELTA_MONO if monoisotopic else C_TERM_AMIDE_DELTA_AVG

    # Disulfide bond oxidation (loss of 2 Hydrogens per bridge)
    if num_disulfide_bonds > 0:
        h2_mass = 2.01565 if monoisotopic else 2.01588
        mass -= num_disulfide_bonds * h2_mass

    return round(mass, 4 if monoisotopic else 2)


def calculate_average_mw(sequence: str, **kwargs) -> float:
    """Compute average molecular weight in Daltons."""
    return calculate_molecular_weight(sequence, monoisotopic=False, **kwargs)


def calculate_monoisotopic_mw(sequence: str, **kwargs) -> float:
    """Compute monoisotopic molecular weight in Daltons."""
    return calculate_molecular_weight(sequence, monoisotopic=True, **kwargs)


def calculate_net_charge(
    sequence: str,
    ph: float,
    pka_scale: str = "EMBOSS",
    n_terminal: str = "free",
    c_terminal: str = "free",
) -> float:
    """
    Calculate net electrostatic charge of a peptide at a specific pH.

    Uses the multi-protic Henderson-Hasselbalch equation.

    Args:
        sequence: Amino acid sequence.
        ph: Solution pH (typically between 0.0 and 14.0).
        pka_scale: Name of pKa scale ('EMBOSS', 'LEHNINGER', 'BJELLQVIST').
        n_terminal: 'free' or 'acetyl'.
        c_terminal: 'free' or 'amide'.

    Returns:
        Net electrostatic charge in elementary charge units.
    """
    seq = validate_sequence(sequence)
    scale_name = pka_scale.upper()
    if scale_name not in PKA_SCALES:
        raise ValueError(
            f"Unknown pKa scale '{pka_scale}'. Supported scales: {list(PKA_SCALES.keys())}"
        )

    scale = PKA_SCALES[scale_name]

    # Positive basic groups: fraction protonated = 1 / (1 + 10^(pH - pKa))
    charge_pos = 0.0
    if n_terminal.lower() not in ("acetyl", "ac"):
        charge_pos += 1.0 / (1.0 + 10.0 ** (ph - scale["N_term"]))

    for aa in ("K", "R", "H"):
        count = seq.count(aa)
        if count > 0:
            charge_pos += count / (1.0 + 10.0 ** (ph - scale[aa]))

    # Negative acidic groups: fraction deprotonated = 1 / (1 + 10^(pKa - pH))
    charge_neg = 0.0
    if c_terminal.lower() not in ("amide", "nh2"):
        charge_neg += 1.0 / (1.0 + 10.0 ** (scale["C_term"] - ph))

    for aa in ("D", "E", "C", "Y"):
        count = seq.count(aa)
        if count > 0:
            charge_neg += count / (1.0 + 10.0 ** (scale[aa] - ph))

    return charge_pos - charge_neg


def titrate_charge(
    sequence: str,
    pka_scale: str = "EMBOSS",
    ph_min: float = 0.0,
    ph_max: float = 14.0,
    step: float = 0.1,
    n_terminal: str = "free",
    c_terminal: str = "free",
) -> Dict[float, float]:
    """
    Generate the complete net charge vs pH titration curve.

    Args:
        sequence: Amino acid sequence.
        pka_scale: pKa scale name.
        ph_min: Minimum pH (default 0.0).
        ph_max: Maximum pH (default 14.0).
        step: pH step interval (default 0.1).

    Returns:
        Dictionary mapping rounded pH value to net charge.
    """
    validate_sequence(sequence)
    titration: Dict[float, float] = {}

    current_ph = ph_min
    while current_ph <= ph_max + 1e-9:
        ph_round = round(current_ph, 2)
        charge = calculate_net_charge(
            sequence,
            ph_round,
            pka_scale=pka_scale,
            n_terminal=n_terminal,
            c_terminal=c_terminal,
        )
        titration[ph_round] = round(charge, 4)
        current_ph += step

    return titration


def calculate_isoelectric_point(
    sequence: str,
    pka_scale: str = "EMBOSS",
    n_terminal: str = "free",
    c_terminal: str = "free",
    tolerance: float = 1e-5,
    max_iterations: int = 100,
) -> float:
    """
    Compute isoelectric point (pI) where net charge Q(pH) == 0.

    Uses bisection root-finding over [0.0, 14.0] guaranteed by charge monotonicity.

    Args:
        sequence: Amino acid sequence.
        pka_scale: pKa scale name.
        n_terminal: Terminal cap state.
        c_terminal: Terminal cap state.
        tolerance: Absolute charge convergence threshold (|Q| < tolerance).
        max_iterations: Maximum bisection steps.

    Returns:
        Isoelectric point (pI) in pH units.
    """
    seq = validate_sequence(sequence)

    low = 0.0
    high = 14.0

    q_low = calculate_net_charge(seq, low, pka_scale, n_terminal, c_terminal)
    q_high = calculate_net_charge(seq, high, pka_scale, n_terminal, c_terminal)

    if q_low <= 0.0:
        return 0.0
    if q_high >= 0.0:
        return 14.0

    mid = 7.0
    for _ in range(max_iterations):
        mid = (low + high) / 2.0
        q_mid = calculate_net_charge(seq, mid, pka_scale, n_terminal, c_terminal)

        if abs(q_mid) < tolerance or (high - low) < 1e-5:
            return round(mid, 2)

        if q_mid > 0.0:
            low = mid
        else:
            high = mid

    return round(mid, 2)


def calculate_gravy(sequence: str) -> float:
    """
    Calculate Grand Average of Hydropathicity (GRAVY) via Kyte-Doolittle scale.

    Args:
        sequence: Amino acid sequence.

    Returns:
        GRAVY score (positive = hydrophobic, negative = hydrophilic).
    """
    seq = validate_sequence(sequence)
    total_hydropathy = sum(KYTE_DOOLITTLE[aa] for aa in seq)
    return round(total_hydropathy / len(seq), 4)


def calculate_sliding_hydropathy(sequence: str, window_size: int = 5) -> List[float]:
    """
    Calculate sliding-window hydropathy profile across the sequence.

    Args:
        sequence: Amino acid sequence.
        window_size: Window size (must be >= 1 and <= len(sequence)).

    Returns:
        List of window-averaged hydropathy values.
    """
    seq = validate_sequence(sequence)
    n = len(seq)
    if window_size < 1:
        raise ValueError(f"Window size must be >= 1, got {window_size}")
    if window_size > n:
        window_size = n

    profile: List[float] = []
    half = window_size // 2

    for i in range(n):
        start = max(0, i - half)
        end = min(n, start + window_size)
        if end - start < window_size and start > 0:
            start = max(0, end - window_size)
        window = seq[start:end]
        avg_h = sum(KYTE_DOOLITTLE[aa] for aa in window) / len(window)
        profile.append(round(avg_h, 3))

    return profile


def calculate_hydrophobic_moment(
    sequence: str,
    angle_deg: float = 100.0,
    conformation: Optional[str] = None,
) -> float:
    """
    Calculate the Eisenberg consensus amphiphilic hydrophobic moment (mu_H).

    Args:
        sequence: Amino acid sequence.
        angle_deg: Rotation angle per residue in degrees (default 100.0 for alpha-helix).
        conformation: Optional shortcut ('alpha': 100 deg, '3_10': 120 deg, 'beta': 160 deg).

    Returns:
        Hydrophobic moment magnitude normalized by sequence length.
    """
    seq = validate_sequence(sequence)

    if conformation:
        conf_lower = conformation.lower()
        if conf_lower in ("alpha", "alpha_helix", "a"):
            angle_deg = 100.0
        elif conf_lower in ("3_10", "310", "3_10_helix"):
            angle_deg = 120.0
        elif conf_lower in ("beta", "beta_sheet", "b"):
            angle_deg = 160.0
        else:
            raise ValueError(f"Unknown conformation '{conformation}'. Use 'alpha', '3_10', or 'beta'.")

    rad = math.radians(angle_deg)
    cos_sum = 0.0
    sin_sum = 0.0

    for i, aa in enumerate(seq, 1):
        h = EISENBERG_CONSENSUS[aa]
        cos_sum += h * math.cos(i * rad)
        sin_sum += h * math.sin(i * rad)

    moment = math.sqrt(cos_sum**2 + sin_sum**2) / len(seq)
    return round(moment, 4)


def calculate_instability_index(sequence: str) -> float:
    """
    Calculate the Guruprasad et al. (1990) Instability Index (II).

    A value <= 40.0 indicates a stable peptide; > 40.0 indicates unstable.

    Args:
        sequence: Amino acid sequence.

    Returns:
        Instability index score.
    """
    seq = validate_sequence(sequence)
    length = len(seq)

    if length < 2:
        return 0.0

    total_diwv = 0.0
    for i in range(length - 1):
        aa1 = seq[i]
        aa2 = seq[i + 1]
        weight = DIWV_MATRIX[aa1][aa2]
        total_diwv += weight

    instability_index = (10.0 / length) * total_diwv
    return round(instability_index, 2)


def calculate_aliphatic_index(sequence: str) -> float:
    """
    Calculate Ikai (1980) Aliphatic Index (AI).

    AI = X(Ala) + 2.9 * X(Val) + 3.9 * (X(Ile) + X(Leu))
    where X(aa) is the mole percentage of residue aa.

    Args:
        sequence: Amino acid sequence.

    Returns:
        Aliphatic index value.
    """
    seq = validate_sequence(sequence)
    length = len(seq)

    x_ala = (seq.count("A") / length) * 100.0
    x_val = (seq.count("V") / length) * 100.0
    x_ile = (seq.count("I") / length) * 100.0
    x_leu = (seq.count("L") / length) * 100.0

    ai = x_ala + 2.9 * x_val + 3.9 * (x_ile + x_leu)
    return round(ai, 2)


def predict_proteolytic_cleavage(
    sequence: str,
    proteases: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Predict enzymatic cleavage sites and resulting fragments for endopeptidases.

    Rules strictly enforced:
    - Trypsin: [KR](?!P) (Cleaves after K or R, except when followed by Proline!)
    - Chymotrypsin: [FYW](?!P) (Cleaves after aromatic residues, except before Proline)
    - Pepsin: [FLWY](?!P) (Primary hydrophobic cleavage)
    - Proteinase K: [AILVF](?!P) (Broad aliphatic/aromatic cleavage)
    - Endoproteinase Lys-C: K(?!P)
    - Endoproteinase Asp-N: Cleaves N-terminal to [DE]
    - DPP-IV: N-terminal dipeptide cleavage if second residue is P or A (^.[PA])

    Args:
        sequence: Amino acid sequence.
        proteases: List of protease names to evaluate (None evaluates all).

    Returns:
        Dict containing cleavage sites (1-based), resulting fragments, and vulnerability score.
    """
    seq = validate_sequence(sequence)
    length = len(seq)

    all_supported = [
        "trypsin",
        "chymotrypsin",
        "pepsin",
        "proteinase_k",
        "lys_c",
        "asp_n",
        "dpp4",
    ]

    selected = [p.lower() for p in (proteases or all_supported)]
    for p in selected:
        if p not in all_supported:
            raise ValueError(f"Unknown protease '{p}'. Supported proteases: {all_supported}")

    sites: Dict[str, List[int]] = {}
    fragments: Dict[str, List[str]] = {}

    for protease in selected:
        cut_positions: List[int] = []

        if protease == "trypsin":
            for i in range(length - 1):
                if seq[i] in ("K", "R") and seq[i + 1] != "P":
                    cut_positions.append(i + 1)

        elif protease == "chymotrypsin":
            for i in range(length - 1):
                if seq[i] in ("F", "Y", "W") and seq[i + 1] != "P":
                    cut_positions.append(i + 1)

        elif protease == "pepsin":
            for i in range(length - 1):
                if seq[i] in ("F", "L", "W", "Y") and seq[i + 1] != "P":
                    cut_positions.append(i + 1)

        elif protease == "proteinase_k":
            for i in range(length - 1):
                if seq[i] in ("A", "I", "L", "V", "F") and seq[i + 1] != "P":
                    cut_positions.append(i + 1)

        elif protease == "lys_c":
            for i in range(length - 1):
                if seq[i] == "K" and seq[i + 1] != "P":
                    cut_positions.append(i + 1)

        elif protease == "asp_n":
            for i in range(1, length):
                if seq[i] in ("D", "E"):
                    cut_positions.append(i)

        elif protease == "dpp4":
            if length >= 3 and seq[1] in ("P", "A"):
                cut_positions.append(2)

        # Generate fragments from cut positions
        frags: List[str] = []
        if not cut_positions:
            frags = [seq]
        else:
            prev = 0
            for pos in sorted(set(cut_positions)):
                frags.append(seq[prev:pos])
                prev = pos
            if prev < length:
                frags.append(seq[prev:length])

        sites[protease] = cut_positions
        fragments[protease] = frags

    total_sites = sum(len(s) for s in sites.values())
    vulnerability_score = round(total_sites / length, 4)

    return {
        "cleavage_sites": sites,
        "fragments": fragments,
        "total_cleavage_sites": total_sites,
        "vulnerability_score": vulnerability_score,
    }


@dataclass
class PeptideProperties:
    """Consolidated biophysical properties of a peptide."""

    sequence: str
    length: int
    molecular_weight: float
    monoisotopic_mass: float
    isoelectric_point: float
    gravy: float
    hydrophobic_moment: float
    instability_index: float
    is_stable: bool
    aliphatic_index: float
    charge_at_ph_7_4: float
    charge_at_ph_7_0: float
    proteolytic_cleavage_sites: Dict[str, List[int]]
    pka_scale: str
    n_terminal: str = "free"
    c_terminal: str = "free"
    fragments: Dict[str, List[str]] = field(default_factory=dict)
    vulnerability_score: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert properties to dictionary."""
        return {
            "sequence": self.sequence,
            "length": self.length,
            "molecular_weight": self.molecular_weight,
            "monoisotopic_mass": self.monoisotopic_mass,
            "isoelectric_point": self.isoelectric_point,
            "gravy": self.gravy,
            "hydrophobic_moment": self.hydrophobic_moment,
            "instability_index": self.instability_index,
            "is_stable": self.is_stable,
            "aliphatic_index": self.aliphatic_index,
            "charge_at_ph_7_4": self.charge_at_ph_7_4,
            "charge_at_ph_7_0": self.charge_at_ph_7_0,
            "pka_scale": self.pka_scale,
            "n_terminal": self.n_terminal,
            "c_terminal": self.c_terminal,
            "proteolytic_cleavage_sites": self.proteolytic_cleavage_sites,
            "fragments": self.fragments,
            "vulnerability_score": self.vulnerability_score,
        }


def calculate_properties(
    sequence: str,
    pka_scale: str = "EMBOSS",
    n_terminal: str = "free",
    c_terminal: str = "free",
    num_disulfide_bonds: int = 0,
) -> PeptideProperties:
    """
    Calculate full biophysical property profile for an amino acid sequence.

    Args:
        sequence: Amino acid sequence.
        pka_scale: pKa scale for pI and titration ('EMBOSS', 'LEHNINGER', 'BJELLQVIST').
        n_terminal: 'free' or 'acetyl'.
        c_terminal: 'free' or 'amide'.
        num_disulfide_bonds: Number of disulfide bridges.

    Returns:
        PeptideProperties dataclass instance.
    """
    seq = validate_sequence(sequence)

    avg_mw = calculate_molecular_weight(
        seq,
        monoisotopic=False,
        n_terminal=n_terminal,
        c_terminal=c_terminal,
        num_disulfide_bonds=num_disulfide_bonds,
    )
    mono_mw = calculate_molecular_weight(
        seq,
        monoisotopic=True,
        n_terminal=n_terminal,
        c_terminal=c_terminal,
        num_disulfide_bonds=num_disulfide_bonds,
    )

    pi = calculate_isoelectric_point(
        seq, pka_scale=pka_scale, n_terminal=n_terminal, c_terminal=c_terminal
    )
    gravy = calculate_gravy(seq)
    mu_h = calculate_hydrophobic_moment(seq, angle_deg=100.0)
    ii = calculate_instability_index(seq)
    ai = calculate_aliphatic_index(seq)

    q74 = round(calculate_net_charge(seq, 7.4, pka_scale=pka_scale, n_terminal=n_terminal, c_terminal=c_terminal), 2)
    q70 = round(calculate_net_charge(seq, 7.0, pka_scale=pka_scale, n_terminal=n_terminal, c_terminal=c_terminal), 2)

    cleavage_data = predict_proteolytic_cleavage(seq)

    return PeptideProperties(
        sequence=seq,
        length=len(seq),
        molecular_weight=avg_mw,
        monoisotopic_mass=mono_mw,
        isoelectric_point=pi,
        gravy=gravy,
        hydrophobic_moment=mu_h,
        instability_index=ii,
        is_stable=(ii <= 40.0),
        aliphatic_index=ai,
        charge_at_ph_7_4=q74,
        charge_at_ph_7_0=q70,
        proteolytic_cleavage_sites=cleavage_data["cleavage_sites"],
        fragments=cleavage_data["fragments"],
        vulnerability_score=cleavage_data["vulnerability_score"],
        pka_scale=pka_scale.upper(),
        n_terminal=n_terminal,
        c_terminal=c_terminal,
    )
