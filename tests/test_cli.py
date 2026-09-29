"""
Unit and integration tests for the insilico-peptide CLI.
"""

from pathlib import Path

from insilico_peptide.cli.main import main


def test_cli_analyze_sequence_table() -> None:
    """Verify CLI analyze command execution with direct sequence string."""
    exit_code = main(["analyze", "--sequence", "GEPPPGKPADDAGLV", "--name", "BPC-157"])
    assert exit_code == 0


def test_cli_analyze_json_and_markdown_exports(tmp_path: Path) -> None:
    """Verify CLI exports structured JSON and publication Markdown."""
    json_path = tmp_path / "bpc.json"
    md_path = tmp_path / "bpc.md"
    plot_charge_path = tmp_path / "bpc_charge.png"

    exit_code = main([
        "analyze",
        "--sequence", "GEPPPGKPADDAGLV",
        "--name", "BPC-157",
        "--export-json", str(json_path),
        "--export-markdown", str(md_path),
        "--plot-charge", str(plot_charge_path),
    ])

    assert exit_code == 0
    assert json_path.exists()
    assert md_path.exists()
    assert plot_charge_path.exists()

    # Verify JSON content
    import json
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    assert data["sequence"] == "GEPPPGKPADDAGLV"
    assert data["molecular_weight"] == 1419.55


def test_cli_analyze_fasta(tmp_path: Path) -> None:
    """Verify CLI ingestion and parsing of multi-sequence FASTA files."""
    fasta_file = tmp_path / "test_peptides.fasta"
    fasta_file.write_text(
        ">BPC157 Body protection compound\n"
        "GEPPPGKPADDAGLV\n"
        ">GHK Regenerative tripeptide\n"
        "GHK\n"
    )

    exit_code = main(["analyze", "--fasta", str(fasta_file)])
    assert exit_code == 0


def test_cli_scan_alanine(tmp_path: Path) -> None:
    """Verify CLI scan command with alanine mode and plot export."""
    plot_path = tmp_path / "scan_plot.png"
    exit_code = main([
        "scan",
        "--sequence", "GEPPPGKPADDAGLV",
        "--mode", "alanine",
        "--top", "5",
        "--plot-deltas", str(plot_path),
    ])
    assert exit_code == 0
    assert plot_path.exists()


def test_cli_fetch_offline(tmp_path: Path) -> None:
    """Verify CLI fetch command with offline fallback and JSON export."""
    out_json = tmp_path / "fetch_out.json"
    exit_code = main([
        "fetch",
        "--uniprot", "P01275",
        "--alphafold",
        "--chembl", "CHEMBL1784",
        "--offline",
        "--export-json", str(out_json),
    ])
    assert exit_code == 0
    assert out_json.exists()


def test_cli_error_missing_arguments() -> None:
    """Verify CLI fails gracefully when no arguments or invalid sequences are passed."""
    # No subcommand
    assert main([]) != 0

    # Analyze with neither sequence nor fasta
    assert main(["analyze"]) != 0

    # Scan with neither sequence nor fasta
    assert main(["scan"]) != 0


def test_cli_error_invalid_fasta(tmp_path: Path) -> None:
    """Verify CLI handles missing or corrupt FASTA files with non-zero exit code."""
    missing_fasta = tmp_path / "missing.fasta"
    assert main(["analyze", "--fasta", str(missing_fasta)]) != 0

    # Corrupt fasta
    corrupt_fasta = tmp_path / "corrupt.fasta"
    corrupt_fasta.write_text("SEQUENCE_WITHOUT_HEADER\n")
    assert main(["analyze", "--fasta", str(corrupt_fasta)]) != 0


def test_cli_analyze_plots_and_formats(tmp_path: Path) -> None:
    """Verify CLI generation of all plot types and output formats."""
    p_hydro = tmp_path / "hydro.png"
    p_heli = tmp_path / "helical.png"

    # Analyze with plots and markdown stdout
    code = main([
        "analyze",
        "--sequence", "GHK",
        "--output-format", "markdown",
        "--plot-hydropathy", str(p_hydro),
        "--plot-helical", str(p_heli),
    ])
    assert code == 0
    assert p_hydro.exists()
    assert p_heli.exists()

    # Analyze with json stdout
    code_json = main(["analyze", "--sequence", "GHK", "--output-format", "json"])
    assert code_json == 0


def test_cli_scan_homologous_and_exports(tmp_path: Path) -> None:
    """Verify CLI scan in homologous mode with markdown and json exports."""
    json_out = tmp_path / "homo.json"
    md_out = tmp_path / "homo.md"

    code = main([
        "scan",
        "--sequence", "GHK",
        "--mode", "homologous",
        "--export-json", str(json_out),
        "--export-markdown", str(md_out),
    ])
    assert code == 0
    assert json_out.exists()
    assert md_out.exists()


def test_cli_scan_fasta(tmp_path: Path) -> None:
    """Verify CLI scan with FASTA input."""
    fpath = tmp_path / "pep.fasta"
    fpath.write_text(">P1\nGHK\n")
    code = main(["scan", "--fasta", str(fpath), "--mode", "alanine"])
    assert code == 0

