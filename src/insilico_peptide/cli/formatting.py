"""
Rich ANSI terminal formatting and visualization for peptide profiles.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from insilico_peptide.core.biophysics import PeptideProperties
from insilico_peptide.screening.alanine_scan import MutantProfile

console = Console()


def print_biophysical_summary(properties: PeptideProperties, name: Optional[str] = None) -> None:
    """Print publication-grade summary card to the terminal using rich."""
    peptide_name = name or properties.sequence

    # Header Panel
    header_text = Text()
    header_text.append(f"Peptide: {peptide_name}\n", style="bold cyan")
    header_text.append(f"Sequence: {properties.sequence} ({properties.length} AA)\n", style="bold white")
    header_text.append(f"pKa Scale: {properties.pka_scale}", style="italic dim")

    console.print(Panel(header_text, title="[bold green]Insilico Peptide Analytical Profile[/bold green]", expand=False))

    # Core Biophysics Table
    table = Table(title="Core Physicochemical Metrics", header_style="bold magenta", expand=False)
    table.add_column("Property", style="dim", width=30)
    table.add_column("Value", justify="right", style="bold white", width=18)
    table.add_column("Interpretation", style="italic", width=36)

    table.add_row(
        "Average Molecular Weight",
        f"{properties.molecular_weight:.2f} Da",
        "Isotopic average mass",
    )
    table.add_row(
        "Monoisotopic Mass",
        f"{properties.monoisotopic_mass:.4f} Da",
        "High-resolution mass spec standard",
    )
    table.add_row(
        "Isoelectric Point (pI)",
        f"{properties.isoelectric_point:.2f}",
        "Neutral charge pH (|Q| < 10⁻⁵)",
    )
    table.add_row(
        "Net Charge at pH 7.4",
        f"{properties.charge_at_ph_7_4:+.2f} e",
        "Physiological pH charge state",
    )
    table.add_row(
        "Net Charge at pH 7.0",
        f"{properties.charge_at_ph_7_0:+.2f} e",
        "Neutral pH charge state",
    )

    gravy_style = "bold green" if properties.gravy <= 0 else "bold yellow"
    gravy_desc = "Net Hydrophilic (Water soluble)" if properties.gravy <= 0 else "Net Hydrophobic"
    table.add_row(
        "GRAVY Hydropathicity",
        f"[{gravy_style}]{properties.gravy:+.4f}[/{gravy_style}]",
        gravy_desc,
    )

    moment_desc = "Strong amphiphilic helix" if properties.hydrophobic_moment > 0.15 else "Low/moderate amphiphilicity"
    table.add_row(
        "Hydrophobic Moment (μ_H, 100°)",
        f"{properties.hydrophobic_moment:.4f}",
        moment_desc,
    )

    status_tag = "[bold green]Stable[/bold green]" if properties.is_stable else "[bold red]Unstable[/bold red]"
    table.add_row(
        "Guruprasad Instability Index (II)",
        f"{properties.instability_index:.2f}",
        f"{status_tag} (Cutoff: 40.0)",
    )
    table.add_row(
        "Aliphatic Index (AI)",
        f"{properties.aliphatic_index:.2f}",
        "Ikai thermostability metric",
    )
    table.add_row(
        "Vulnerability Score",
        f"{properties.vulnerability_score:.4f}",
        "Proteolytic cut sites / residue",
    )

    console.print(table)

    # Cleavage summary table
    cleave_table = Table(title="Endopeptidase Cleavage Vulnerability", header_style="bold yellow", expand=False)
    cleave_table.add_column("Protease", style="bold cyan", width=18)
    cleave_table.add_column("Cut Sites (1-based)", style="bold white", width=24)
    cleave_table.add_column("Status", width=20)

    for protease, sites in properties.proteolytic_cleavage_sites.items():
        if sites:
            sites_str = ", ".join(str(s) for s in sites)
            status = f"[bold red]Cleaves ({len(sites)} sites)[/bold red]"
        else:
            sites_str = "None"
            status = "[bold green]Resistant[/bold green]"
        cleave_table.add_row(protease.upper(), sites_str, status)

    console.print(cleave_table)


def print_alanine_scan_table(mutants: List[MutantProfile], top_n: Optional[int] = None) -> None:
    """Print alanine scanning mutational table with hotspot highlights."""
    table = Table(title="In-Silico Mutational Screening & Hotspot Scan", header_style="bold blue", expand=False)
    table.add_column("Pos", justify="center", style="dim", width=5)
    table.add_column("Mutation", justify="center", style="bold white", width=10)
    table.add_column("ΔMW", justify="right", width=9)
    table.add_column("ΔpI", justify="right", width=8)
    table.add_column("ΔGRAVY", justify="right", width=10)
    table.add_column("Δμ_H", justify="right", width=10)
    table.add_column("ΔII", justify="right", width=8)
    table.add_column("ΔΔG (kcal/mol)", justify="right", width=14)
    table.add_column("ΔSol", justify="right", width=8)
    table.add_column("Rank", justify="center", width=8)
    table.add_column("Classification", width=16)

    display_mutants = mutants[:top_n] if top_n else mutants

    for m in display_mutants:
        status = "[bold red]HOTSPOT[/bold red]" if m.is_hotspot else "[green]Neutral[/green]"
        table.add_row(
            str(m.position),
            f"`{m.wt_residue}{m.position}{m.mut_residue}`",
            f"{m.delta_mw:+.2f}",
            f"{m.delta_pi:+.2f}",
            f"{m.delta_gravy:+.4f}",
            f"{m.delta_moment:+.4f}",
            f"{m.delta_instability:+.2f}",
            f"{m.ddg_pred:+.2f}",
            f"{m.solubility_pred:+.2f}",
            f"Rank {m.pareto_rank}",
            status,
        )

    console.print(table)


def print_multi_omics_card(data: Dict[str, Any]) -> None:
    """Print multi-omic summary card."""
    card_text = Text()
    if "uniprot" in data:
        u = data["uniprot"]
        card_text.append("UniProt Knowledgebase:\n", style="bold cyan")
        card_text.append(f"  Accession: {u.get('accession')} ({u.get('entry_name')})\n")
        card_text.append(f"  Protein: {u.get('protein_name')}\n")
        card_text.append(f"  Organism: {u.get('organism')}\n\n")

    if "alphafold" in data:
        af = data["alphafold"]
        card_text.append("AlphaFold Structure DB:\n", style="bold green")
        card_text.append(f"  Entry: {af.get('entry_id')}\n")
        card_text.append(f"  Mean pLDDT: {af.get('mean_plddt'):.2f}\n")
        card_text.append(f"  PDB URL: {af.get('pdb_url')}\n\n")

    if "chembl" in data:
        ch = data["chembl"]
        card_text.append("ChEMBL Target Bioactivity:\n", style="bold yellow")
        card_text.append(f"  Target: {ch.get('target_chembl_id')} — {ch.get('pref_name')}\n")
        card_text.append(f"  Type: {ch.get('target_type')} | Organism: {ch.get('organism')}\n")

    console.print(Panel(card_text, title="[bold magenta]Multi-Omic Integration[/bold magenta]", expand=False))
