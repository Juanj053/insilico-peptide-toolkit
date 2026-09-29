"""
Headless Matplotlib scientific plotting engine for biophysical profiles and screening vectors.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import List

import matplotlib

# Guarantee headless backend for CLI / CI / background tasks
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from insilico_peptide.core.biophysics import (
    calculate_isoelectric_point,
    calculate_sliding_hydropathy,
    titrate_charge,
    validate_sequence,
)
from insilico_peptide.core.scales import EISENBERG_CONSENSUS
from insilico_peptide.screening.alanine_scan import MutantProfile


def plot_charge_titration_curve(
    sequence: str,
    output_path: str,
    pka_scale: str = "EMBOSS",
    dpi: int = 300,
) -> str:
    """
    Render net charge vs pH titration curve (pH 0 to 14).

    Args:
        sequence: Amino acid sequence.
        output_path: Output image path (.png, .pdf, .svg).
        pka_scale: pKa scale name.
        dpi: Output resolution.

    Returns:
        Absolute string path to saved plot.
    """
    seq = validate_sequence(sequence)
    titration = titrate_charge(seq, pka_scale=pka_scale, ph_min=0.0, ph_max=14.0, step=0.1)
    pi = calculate_isoelectric_point(seq, pka_scale=pka_scale)

    phs = sorted(titration.keys())
    charges = [titration[ph] for ph in phs]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(phs, charges, color="#0A3D36", linewidth=2.5, label=f"Net Charge ({pka_scale})")
    ax.axhline(0, color="#C5A059", linestyle="--", linewidth=1.2, alpha=0.8, label="Net Zero (pI)")
    ax.axvline(pi, color="#8F6E40", linestyle=":", linewidth=1.5, label=f"pI = {pi:.2f}")
    ax.axvline(7.4, color="#64748B", linestyle="--", linewidth=1.0, alpha=0.6, label="Physiological pH 7.4")

    ax.set_title(f"Net Charge Titration Curve: {seq}", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Solution pH", fontsize=11)
    ax.set_ylabel("Net Electrostatic Charge (e)", fontsize=11)
    ax.set_xlim(0, 14)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", frameon=True)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    plt.close(fig)
    return str(out.resolve())


def plot_hydropathy_profile(
    sequence: str,
    output_path: str,
    window_size: int = 5,
    dpi: int = 300,
) -> str:
    """
    Render Kyte-Doolittle sliding window hydropathy profile.

    Args:
        sequence: Amino acid sequence.
        output_path: Output image path.
        window_size: Sliding window length.
        dpi: Plot resolution.

    Returns:
        Absolute path to saved plot.
    """
    seq = validate_sequence(sequence)
    profile = calculate_sliding_hydropathy(seq, window_size=window_size)
    positions = list(range(1, len(seq) + 1))

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(positions, profile, color="#0A3D36", linewidth=2.0, marker="o", markersize=4)
    ax.axhline(0, color="#C5A059", linestyle="--", linewidth=1.2)

    # Fill hydrophobic vs hydrophilic zones
    ax.fill_between(
        positions,
        profile,
        0,
        where=[h >= 0 for h in profile],
        color="#FDE68A",
        alpha=0.5,
        interpolate=True,
        label="Hydrophobic core",
    )
    ax.fill_between(
        positions,
        profile,
        0,
        where=[h < 0 for h in profile],
        color="#93C5FD",
        alpha=0.4,
        interpolate=True,
        label="Hydrophilic / Soluble",
    )

    ax.set_title(f"Kyte-Doolittle Hydropathy Profile (W={window_size}): {seq}", fontsize=12, fontweight="bold", pad=12)
    ax.set_xlabel("Residue Position (1-based)", fontsize=11)
    ax.set_ylabel("Hydropathy Index", fontsize=11)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{seq[i-1]}\n{i}" for i in positions], fontsize=8)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="best", frameon=True)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    plt.close(fig)
    return str(out.resolve())


def plot_helical_wheel(
    sequence: str,
    output_path: str,
    dpi: int = 300,
) -> str:
    """
    Render 2D polar projection (helical wheel) with hydrophobic moment vector.

    Args:
        sequence: Amino acid sequence.
        output_path: Output image path.
        dpi: Plot resolution.

    Returns:
        Absolute path to saved plot.
    """
    seq = validate_sequence(sequence)
    n = len(seq)

    angles_rad: List[float] = []
    radii: List[float] = []
    colors: List[str] = []

    cos_sum = 0.0
    sin_sum = 0.0

    # Color palette
    basic_res = {"K", "R", "H"}
    acidic_res = {"D", "E"}
    hydrophobic_res = {"A", "V", "I", "L", "M", "F", "W", "Y"}

    for i, aa in enumerate(seq, 1):
        angle = math.radians((i * 100.0) % 360.0)
        radius = 1.0 + (i * 0.04)  # Slight outward spiral
        angles_rad.append(angle)
        radii.append(radius)

        h = EISENBERG_CONSENSUS.get(aa, 0.0)
        cos_sum += h * math.cos(angle)
        sin_sum += h * math.sin(angle)

        if aa in basic_res:
            colors.append("#2563EB")  # Blue
        elif aa in acidic_res:
            colors.append("#DC2626")  # Red
        elif aa in hydrophobic_res:
            colors.append("#D97706")  # Gold/Amber
        else:
            colors.append("#059669")  # Emerald polar

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": "polar"})
    ax.scatter(angles_rad, radii, c=colors, s=350, alpha=0.85, edgecolors="#0A1118", linewidth=1.5, zorder=3)

    for i, (ang, r) in enumerate(zip(angles_rad, radii), 1):
        aa = seq[i - 1]
        ax.text(ang, r, f"{aa}\n{i}", ha="center", va="center", fontsize=8, fontweight="bold", color="white", zorder=4)

    # Plot Hydrophobic Moment vector
    moment_mag = math.sqrt(cos_sum**2 + sin_sum**2) / n
    moment_angle = math.atan2(sin_sum, cos_sum)
    arrow_r = min(1.8, max(0.4, moment_mag * 3.0))

    ax.annotate(
        "",
        xy=(moment_angle, arrow_r),
        xytext=(0, 0),
        arrowprops=dict(facecolor="#7C2D12", edgecolor="#7C2D12", width=2.5, headwidth=8),
        zorder=5,
    )
    ax.text(
        moment_angle,
        arrow_r + 0.15,
        f"μ_H = {moment_mag:.3f}",
        fontsize=9,
        fontweight="bold",
        color="#7C2D12",
        ha="center",
    )

    ax.set_title(f"Helical Wheel Projection (α-helix, 100°)\n{seq[:15]}...", fontsize=11, fontweight="bold", pad=15)
    ax.set_yticklabels([])
    ax.grid(True, linestyle=":", alpha=0.5)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    plt.close(fig)
    return str(out.resolve())


def plot_alanine_deltas(
    mutants: List[MutantProfile],
    output_path: str,
    dpi: int = 300,
) -> str:
    """
    Render bar chart of stability and amphiphilic deltas across alanine scanning mutants.

    Args:
        mutants: List of MutantProfile objects.
        output_path: Path to save plot.
        dpi: Output resolution.

    Returns:
        Saved plot path.
    """
    positions = [m.position for m in mutants]
    labels = [f"{m.wt_residue}{m.position}" for m in mutants]
    ddgs = [m.ddg_pred for m in mutants]
    bar_colors = ["#DC2626" if m.is_hotspot else "#0A3D36" for m in mutants]

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(positions, ddgs, color=bar_colors, alpha=0.85, edgecolor="#0A1118", linewidth=1.0)
    ax.axhline(0, color="#C5A059", linestyle="--", linewidth=1.2)
    ax.axhline(1.5, color="#DC2626", linestyle=":", linewidth=1.0, alpha=0.7, label="Destabilization threshold (1.5 kcal/mol)")

    ax.set_title("In-Silico Alanine Scanning: Predicted Thermodynamic Stability Shift (ΔΔG)", fontsize=12, fontweight="bold", pad=12)
    ax.set_xlabel("Residue Position & Wild-Type Identity", fontsize=11)
    ax.set_ylabel("Predicted ΔΔG (kcal/mol)", fontsize=11)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper right", frameon=True)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    plt.close(fig)
    return str(out.resolve())
