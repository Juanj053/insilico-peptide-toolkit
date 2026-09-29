"""
Command-line interface (CLI) for insilico-peptide-toolkit.

Subcommands:
- analyze: Biophysical characterization of peptides from string or FASTA.
- scan: In-silico alanine scanning and homologous mutational screening.
- fetch: Resilient retrieval of multi-omic data from UniProt, AlphaFold DB, and ChEMBL.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from rich.console import Console

from insilico_peptide.cli.export_json import export_to_json
from insilico_peptide.cli.export_markdown import export_to_markdown
from insilico_peptide.cli.formatting import (
    print_alanine_scan_table,
    print_biophysical_summary,
    print_multi_omics_card,
)
from insilico_peptide.cli.plotting import (
    plot_alanine_deltas,
    plot_charge_titration_curve,
    plot_helical_wheel,
    plot_hydropathy_profile,
)
from insilico_peptide.connectors.alphafold import fetch_alphafold_plddt
from insilico_peptide.connectors.chembl import fetch_chembl_bioactivity
from insilico_peptide.connectors.uniprot import fetch_uniprot_data
from insilico_peptide.core.biophysics import calculate_properties, validate_sequence
from insilico_peptide.screening.alanine_scan import run_alanine_scan
from insilico_peptide.screening.mutagenesis import run_homologous_scan

console = Console()


def parse_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    """
    Parse a FASTA file and return list of (header, sequence) tuples.

    Args:
        fasta_path: Path to FASTA formatted file.

    Returns:
        List of (header, sequence) tuples.

    Raises:
        FileNotFoundError: If fasta_path does not exist.
        ValueError: If file is empty or malformed.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    records: List[Tuple[str, str]] = []
    current_header: Optional[str] = None
    current_seq_parts: List[str] = []

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line_no, raw_line in enumerate(f, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    full_seq = "".join(current_seq_parts)
                    if not full_seq:
                        raise ValueError(f"Empty sequence for entry '{current_header}' in {fasta_path}")
                    records.append((current_header, full_seq))
                current_header = line[1:].strip()
                current_seq_parts = []
            else:
                if current_header is None:
                    raise ValueError(f"Invalid FASTA format: sequence found before header at line {line_no}")
                current_seq_parts.append(line)

        if current_header is not None:
            full_seq = "".join(current_seq_parts)
            if not full_seq:
                raise ValueError(f"Empty sequence for entry '{current_header}' in {fasta_path}")
            records.append((current_header, full_seq))

    if not records:
        raise ValueError(f"FASTA file is empty or contains no records: {fasta_path}")

    return records


def handle_analyze(args: argparse.Namespace) -> int:
    """Execute 'analyze' subcommand."""
    sequences: List[Tuple[str, str]] = []

    if args.fasta:
        try:
            records = parse_fasta(Path(args.fasta))
            sequences.extend(records)
        except Exception as e:
            console.print(f"[bold red]Error parsing FASTA file:[/bold red] {e}", style="red")
            return 1
    elif args.sequence:
        name = args.name or args.sequence
        sequences.append((name, args.sequence))
    else:
        console.print("[bold red]Error:[/bold red] Either --sequence or --fasta must be provided.", style="red")
        return 1

    all_results = []

    for name, raw_seq in sequences:
        try:
            props = calculate_properties(raw_seq, pka_scale=args.scale)
        except Exception as e:
            console.print(f"[bold red]Error analyzing peptide '{name}':[/bold red] {e}")
            return 1

        all_results.append((name, props))

        if args.output_format == "table":
            print_biophysical_summary(props, name=name)
        elif args.output_format == "json" and not args.export_json:
            import json
            print(json.dumps(props.to_dict(), indent=2))
        elif args.output_format == "markdown" and not args.export_markdown:
            from insilico_peptide.cli.export_markdown import generate_markdown_dossier
            print(generate_markdown_dossier(props, name=name))

        # Plots
        if args.plot_charge:
            p_out = plot_charge_titration_curve(props.sequence, args.plot_charge, pka_scale=args.scale)
            console.print(f"[green]Saved charge titration plot to:[/green] {p_out}")

        if args.plot_hydropathy:
            p_out = plot_hydropathy_profile(props.sequence, args.plot_hydropathy)
            console.print(f"[green]Saved hydropathy profile to:[/green] {p_out}")

        if args.plot_helical:
            p_out = plot_helical_wheel(props.sequence, args.plot_helical)
            console.print(f"[green]Saved helical wheel to:[/green] {p_out}")

    # Exports
    if args.export_json:
        if len(all_results) == 1:
            export_to_json(all_results[0][1].to_dict(), args.export_json)
        else:
            export_to_json([p.to_dict() for _, p in all_results], args.export_json)
        console.print(f"[bold green]Exported JSON dossier to:[/bold green] {args.export_json}")

    if args.export_markdown:
        if len(all_results) == 1:
            export_to_markdown(all_results[0][1], args.export_markdown, name=all_results[0][0])
        else:
            # Multi-entry markdown
            from insilico_peptide.cli.export_markdown import generate_markdown_dossier
            md_full = "\n\n---\n\n".join(generate_markdown_dossier(p, name=n) for n, p in all_results)
            with open(args.export_markdown, "w", encoding="utf-8") as f:
                f.write(md_full)
        console.print(f"[bold green]Exported Markdown dossier to:[/bold green] {args.export_markdown}")

    return 0


def handle_scan(args: argparse.Namespace) -> int:
    """Execute 'scan' subcommand."""
    sequence = ""
    name = "Peptide"

    if args.fasta:
        try:
            records = parse_fasta(Path(args.fasta))
            name, sequence = records[0]
        except Exception as e:
            console.print(f"[bold red]Error parsing FASTA file:[/bold red] {e}")
            return 1
    elif args.sequence:
        sequence = args.sequence
        name = args.name or sequence
    else:
        console.print("[bold red]Error:[/bold red] Either --sequence or --fasta must be provided.", style="red")
        return 1

    try:
        seq = validate_sequence(sequence)
    except Exception as e:
        console.print(f"[bold red]Sequence validation error:[/bold red] {e}")
        return 1

    mode = args.mode.lower()
    if mode == "alanine":
        mutants = run_alanine_scan(seq, pka_scale=args.scale)
    elif mode == "homologous":
        mutants = run_homologous_scan(seq, pka_scale=args.scale)
    else:
        console.print(f"[bold red]Unknown scan mode '{args.mode}'. Use 'alanine' or 'homologous'.[/bold red]")
        return 1

    print_alanine_scan_table(mutants, top_n=args.top)

    if args.export_json:
        data = [m.to_dict() for m in mutants]
        export_to_json(data, args.export_json)
        console.print(f"[bold green]Exported scan JSON to:[/bold green] {args.export_json}")

    if args.export_markdown:
        props = calculate_properties(seq, pka_scale=args.scale)
        export_to_markdown(props, args.export_markdown, name=name, mutants=mutants)
        console.print(f"[bold green]Exported scan Markdown report to:[/bold green] {args.export_markdown}")

    if args.plot_deltas:
        p_out = plot_alanine_deltas(mutants, args.plot_deltas)
        console.print(f"[bold green]Saved mutational delta plot to:[/bold green] {p_out}")

    return 0


def handle_fetch(args: argparse.Namespace) -> int:
    """Execute 'fetch' subcommand."""
    if not (args.uniprot or args.chembl):
        console.print("[bold red]Error:[/bold red] Provide at least --uniprot or --chembl to fetch data.", style="red")
        return 1

    results = {}

    if args.uniprot:
        try:
            u_data = fetch_uniprot_data(args.uniprot, offline=args.offline)
            results["uniprot"] = u_data
            if args.alphafold:
                af_data = fetch_alphafold_plddt(args.uniprot, offline=args.offline)
                results["alphafold"] = af_data
        except Exception as e:
            console.print(f"[bold red]Error fetching UniProt data for '{args.uniprot}':[/bold red] {e}")
            return 1

    if args.chembl:
        try:
            ch_data = fetch_chembl_bioactivity(args.chembl, offline=args.offline)
            results["chembl"] = ch_data
        except Exception as e:
            console.print(f"[bold red]Error fetching ChEMBL data for '{args.chembl}':[/bold red] {e}")
            return 1

    print_multi_omics_card(results)

    if args.export_json:
        export_to_json(results, args.export_json)
        console.print(f"[bold green]Exported multi-omics JSON to:[/bold green] {args.export_json}")

    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct CLI argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog="insilico-peptide",
        description="In-Silico Peptide Biophysical Modeling & Mutational Screening Toolkit (@Juanj053)",
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # Subcommand: analyze
    p_analyze = subparsers.add_parser("analyze", help="Analyze biophysical properties of peptide sequence")
    p_analyze.add_argument("--sequence", "-s", type=str, help="Peptide sequence string")
    p_analyze.add_argument("--fasta", "-f", type=str, help="Path to input FASTA file")
    p_analyze.add_argument("--name", "-n", type=str, help="Peptide identifier or common name")
    p_analyze.add_argument("--ph", type=float, default=7.4, help="Solution pH for net charge (default: 7.4)")
    p_analyze.add_argument("--scale", type=str, default="EMBOSS", choices=["EMBOSS", "LEHNINGER", "BJELLQVIST"], help="pKa scale (default: EMBOSS)")
    p_analyze.add_argument("--output-format", choices=["table", "json", "markdown"], default="table", help="Terminal stdout format")
    p_analyze.add_argument("--export-json", type=str, help="Save full numerical JSON report to path")
    p_analyze.add_argument("--export-markdown", type=str, help="Save publication-grade Markdown dossier to path")
    p_analyze.add_argument("--plot-charge", type=str, help="Save charge vs pH titration plot (.png)")
    p_analyze.add_argument("--plot-hydropathy", type=str, help="Save sliding window hydropathy plot (.png)")
    p_analyze.add_argument("--plot-helical", type=str, help="Save helical wheel projection plot (.png)")

    # Subcommand: scan
    p_scan = subparsers.add_parser("scan", help="In-silico mutational screening and alanine scanning")
    p_scan.add_argument("--sequence", "-s", type=str, help="Wild-type peptide sequence")
    p_scan.add_argument("--fasta", "-f", type=str, help="Path to input FASTA file")
    p_scan.add_argument("--name", "-n", type=str, help="Peptide identifier")
    p_scan.add_argument("--mode", choices=["alanine", "homologous"], default="alanine", help="Scan mode (default: alanine)")
    p_scan.add_argument("--scale", type=str, default="EMBOSS", help="pKa scale (default: EMBOSS)")
    p_scan.add_argument("--top", type=int, default=10, help="Display top N candidates in terminal (default: 10)")
    p_scan.add_argument("--export-json", type=str, help="Save mutant profiles to JSON file")
    p_scan.add_argument("--export-markdown", type=str, help="Save screening report to Markdown file")
    p_scan.add_argument("--plot-deltas", type=str, help="Save mutational delta bar chart (.png)")

    # Subcommand: fetch
    p_fetch = subparsers.add_parser("fetch", help="Fetch multi-omic annotations from UniProt, AlphaFold, ChEMBL")
    p_fetch.add_argument("--uniprot", "-u", type=str, help="UniProt accession ID (e.g. P01275, P01308)")
    p_fetch.add_argument("--alphafold", "-a", action="store_true", help="Fetch AlphaFold structure pLDDT data")
    p_fetch.add_argument("--chembl", "-c", type=str, help="ChEMBL target identifier (e.g. CHEMBL1784)")
    p_fetch.add_argument("--offline", action="store_true", help="Operate in offline mode using local cache and fixtures")
    p_fetch.add_argument("--export-json", type=str, help="Export fetched data to JSON file")

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.subcommand:
        parser.print_help(sys.stderr)
        return 1

    if args.subcommand == "analyze":
        return handle_analyze(args)
    elif args.subcommand == "scan":
        return handle_scan(args)
    elif args.subcommand == "fetch":
        return handle_fetch(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
