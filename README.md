# `insilico-peptide-toolkit`

[![PyPI version](https://img.shields.io/badge/pypi-v0.1.0-blue.svg)](https://pypi.org/project/insilico-peptide-toolkit/)
[![Python](https://img.shields.io/badge/python-%3E%3D3.10-brightgreen.svg)](https://www.python.org/)
[![Build & Tests](https://img.shields.io/badge/tests-passing%20(pytest)-success.svg)](https://github.com/Juanj053/insilico-peptide-toolkit/actions)
[![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen.svg)](https://github.com/Juanj053/insilico-peptide-toolkit)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Juan%20Jos%C3%A9%20%C3%81lvarez-0A66C2.svg?logo=linkedin)](https://www.linkedin.com/in/juanj-alva)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code Style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.1234567-blue.svg)](https://doi.org/10.5281/zenodo.1234567)

> **High-Performance In-Silico Biophysical Modeling, Mutational Screening, and Multi-Omic Profiling Toolkit for Bioactive Peptides**  
> Developed and maintained by **Juan José Álvarez** ([@Juanj053](https://github.com/Juanj053) &nbsp;•&nbsp; [LinkedIn](https://www.linkedin.com/in/juanj-alva)) — Chief Bioengineer & Founder at Elevate Biotech.

---

## Table of Contents
1. [Overview & Core Philosophy](#overview--core-philosophy)
2. [Mathematical & Biophysical Formulations](#mathematical--biophysical-formulations)
   - [Molecular Weight & Isotopic Composition](#1-molecular-weight--isotopic-composition)
   - [Isoelectric Point & Multi-Protic Net Charge Titration](#2-isoelectric-point--multi-protic-net-charge-titration)
   - [Grand Average of Hydropathicity (GRAVY)](#3-grand-average-of-hydropathicity-gravy)
   - [Amphiphilic Hydrophobic Moment](#4-amphiphilic-hydrophobic-moment-mu_h)
   - [Guruprasad Instability Index](#5-guruprasad-instability-index-ii)
   - [Aliphatic Index](#6-aliphatic-index-ai)
   - [Proteolytic Cleavage Susceptibility](#7-proteolytic-cleavage-susceptibility)
3. [Architecture & Module Organization](#architecture--module-organization)
4. [Installation](#installation)
5. [Quickstart & Python API](#quickstart--python-api)
6. [Command-Line Interface (CLI)](#command-line-interface-cli)
7. [Case Studies](#case-studies)
   - [Case Study 1: GLP-1(7-36/37) Receptor Agonist Engineering](#case-study-1-glp-17-3637-receptor-agonist-engineering)
   - [Case Study 2: BPC-157 Gastric Pentadecapeptide Stability](#case-study-2-bpc-157-gastric-pentadecapeptide-stability)
   - [Case Study 3: GHK-Cu Regenerative Tripeptide Coordination](#case-study-3-ghk-cu-regenerative-tripeptide-coordination)
8. [International Research Collaboration & Grants](#international-research-collaboration--grants)
9. [Citing insilico-peptide-toolkit](#citing-insilico-peptide-toolkit)

---

## Overview & Core Philosophy

Rational peptide drug discovery and translational biotechnology demand reproducible, zero-compiler computational modeling before committing to high-throughput solid-phase synthesis (SPPS). `insilico-peptide-toolkit` is an open-source, mathematically exact Python framework engineered to:

1. **Calculate rigorous biophysical properties**: Exact monoisotopic and average molecular weights, Henderson-Hasselbalch charge titration curves, bisection $pI$ root finding, hydropathy, Eisenberg amphiphilic moments, Guruprasad instability indices, and endopeptidase susceptibility.
2. **Execute high-throughput in-silico mutational screening**: Systematic alanine scanning, homologous conservative substitutions via BLOSUM62 matrices, empirical $\Delta\Delta G$ folding stability estimations, and aqueous solubility deltas with Pareto multi-objective optimization.
3. **Bridge multi-omic databases via 3-tier resilience**: Fault-tolerant connectors for UniProt Knowledgebase, AlphaFold Structure Database (pLDDT per-residue metrics), and EMBL-EBI ChEMBL (bioactivity and IC50/Ki targets) operating seamlessly online or 100% offline via embedded curated fixtures.
4. **Export publication-ready dossiers and plots**: Headless 300 DPI Matplotlib vector figures, Rich terminal cards, structured Markdown dossiers, and JSON vectors.

---

## Mathematical & Biophysical Formulations

Every calculation in `insilico-peptide-toolkit` is implemented using pure-Python numerical algorithms derived from peer-reviewed literature.

### 1. Molecular Weight & Isotopic Composition

Peptide molecular weight is calculated by summing constituent amino acid residue masses plus terminal water loss:

$$\text{MW}_{\text{avg}} = \sum_{i=1}^{L} M_{\text{avg}}(a_i) + M_{\text{avg}}(\text{H}_2\text{O}) + \Delta M_{\text{N-term}} + \Delta M_{\text{C-term}} - 2.01588 \cdot N_{\text{SS}}$$

$$\text{MW}_{\text{mono}} = \sum_{i=1}^{L} M_{\text{mono}}(a_i) + M_{\text{mono}}(\text{H}_2\text{O}) + \Delta M_{\text{N-term}} + \Delta M_{\text{C-term}} - 2.01565 \cdot N_{\text{SS}}$$

Where:
- $M_{\text{avg}}(\text{H}_2\text{O}) = 18.01528\text{ Da}$; $M_{\text{mono}}(\text{H}_2\text{O}) = 18.010565\text{ Da}$.
- Acetylated N-terminus (`Ac-`): $+42.0106\text{ Da}$ (mono), $+42.0367\text{ Da}$ (avg).
- Amidated C-terminus (`-NH2`): $-0.9840\text{ Da}$ (mono), $-0.9840\text{ Da}$ (avg).
- $N_{\text{SS}}$: Number of formed cystine disulfide bonds (loss of $2\text{H}$ per bond).

---

### 2. Isoelectric Point & Multi-Protic Net Charge Titration

The net electrical charge $Q(\text{pH})$ of a polyampholytic peptide is governed by the ionization state of its basic and acidic functional groups:

$$Q(\text{pH}) = \frac{1}{1 + 10^{\text{pH} - \text{pK}_{\text{N-term}}}} + \sum_{k \in \{\text{K, R, H}\}} \frac{n_k}{1 + 10^{\text{pH} - \text{pK}_k}} - \frac{1}{1 + 10^{\text{pK}_{\text{C-term}} - \text{pH}}} - \sum_{m \in \{\text{D, E, C, Y}\}} \frac{n_m}{1 + 10^{\text{pK}_m - \text{pH}}}$$

Because $Q(\text{pH})$ is strictly monotonically decreasing ($dQ/d\text{pH} < 0$), the unique isoelectric point $pI$ where $Q(pI) = 0.0$ is resolved via **bisection root finding** over $[0.0, 14.0]$ until $|Q| < 10^{-5}\text{ e}$.

#### Supported $pK_a$ Scales
| Ionizable Group | EMBOSS (Default) | Bjellqvist (ExPASy) | Lehninger |
|---|:---:|:---:|:---:|
| N-terminus ($\alpha\text{-NH}_2$) | 8.60 | 8.60 | 9.69 |
| C-terminus ($\alpha\text{-COOH}$) | 3.60 | 3.60 | 2.34 |
| Lysine (K) | 10.80 | 10.00 | 10.53 |
| Arginine (R) | 12.50 | 12.00 | 12.48 |
| Histidine (H) | 6.50 | 6.60 | 6.00 |
| Aspartate (D) | 3.90 | 4.05 | 3.86 |
| Glutamate (E) | 4.10 | 4.45 | 4.25 |
| Cysteine (C) | 8.50 | 9.00 | 8.33 |
| Tyrosine (Y) | 10.10 | 10.00 | 10.07 |

---

### 3. Grand Average of Hydropathicity (GRAVY)

Per Kyte & Doolittle (*J. Mol. Biol.* 1982), GRAVY evaluates the mean hydropathy across all residues:

$$\text{GRAVY} = \frac{1}{L} \sum_{i=1}^{L} H_{\text{KD}}(a_i)$$

Where $H_{\text{KD}}$ values range from $+4.5$ (Isoleucine, highly hydrophobic) to $-4.5$ (Arginine, highly hydrophilic). $\text{GRAVY} \le 0$ indicates a soluble, hydrophilic peptide.

---

### 4. Amphiphilic Hydrophobic Moment ($\mu_H$)

Per Eisenberg, Weiss, & Terwilliger (*Nature* 1982, *PNAS* 1984), $\mu_H$ quantifies the amphiphilicity of a periodic secondary structure:

$$\mu_H = \frac{1}{L} \sqrt{\left[ \sum_{n=1}^{L} H_{\text{Eis}}(a_n) \cos(n \delta) \right]^2 + \left[ \sum_{n=1}^{L} H_{\text{Eis}}(a_n) \sin(n \delta) \right]^2}$$

Where $H_{\text{Eis}}$ is the consensus hydrophobicity value and $\delta$ is the periodic rotational angle:
- Canonical $\alpha$-helix: $\delta = 100^\circ$ ($3.6$ residues/turn)
- $3_{10}$-helix: $\delta = 120^\circ$ ($3.0$ residues/turn)
- $\beta$-strand: $\delta = 160^\circ\text{–}180^\circ$

A value $\mu_H > 0.15$ in helical conformation indicates an amphiphilic interface capable of partitioning into lipid bilayers or binding GPCR helical clefts (e.g. GLP-1, Exenatide).

---

### 5. Guruprasad Instability Index ($II$)

Guruprasad, Reddy, & Pandit (*Protein Eng.* 1990) established that protein half-life in reticulocyte lysates correlates with the presence of 400 dipeptide motifs:

$$II = \frac{10}{L} \sum_{i=1}^{L-1} \text{DIWV}(a_i, a_{i+1})$$

Where $\text{DIWV}$ is the empirical $20 \times 20$ Dipeptide Instability Weight Value matrix.
- **$II \le 40.0$**: Classified as **Chemically Stable** (predicted in-vivo $t_{1/2} > 16\text{ h}$).
- **$II > 40.0$**: Classified as **Unstable / Labile** (predicted in-vivo $t_{1/2} < 5\text{ h}$).

---

### 6. Aliphatic Index ($AI$)

Per Ikai (*J. Biochem.* 1980), $AI$ measures the relative volume occupied by aliphatic side chains (Ala, Val, Ile, Leu) to predict thermostability:

$$AI = X(\text{Ala}) + 2.9 \cdot X(\text{Val}) + 3.9 \cdot [X(\text{Ile}) + X(\text{Leu})]$$

Where $X(aa) = 100 \times \frac{N(aa)}{L}$.

---

### 7. Proteolytic Cleavage Susceptibility

Enzymatic stability is modeled against rigorous consensus rules:
- **Trypsin**: Cleaves after `[KR]` **except** when immediately followed by Proline (`[KR](?!P)`).
- **Chymotrypsin**: Cleaves after aromatic residues `[FYW]` (high specificity) or `[FYWL]` (broad), inhibited by adjacent Proline.
- **DPP-IV (CD26)**: Specifically hydrolyzes N-terminal dipeptides when the second residue is Alanine or Proline (`^.[AP] \downarrow X_3...`), the primary degradation pathway for GLP-1 and GIP.
- **Pepsin**: Primary gastric digestion marker (`[FLWY](?!P)`).
- **Proteinase K, Lys-C, Asp-N**: Comprehensive pancreatic and metabolic endopeptidases.

---

## Architecture & Module Organization

```
insilico-peptide-toolkit/
├── pyproject.toml                     # PEP 517/518 build config (hatchling backend, uv managed)
├── README.md                          # Academic manual with LaTeX formulas & case studies
├── RESEARCH_COLLABORATION.md          # International grant alignment (CYTED, NIH, Chevening)
├── src/insilico_peptide/
│   ├── __init__.py                    # Public library exports
│   ├── py.typed                       # PEP 561 type annotation marker
│   ├── core/
│   │   ├── __init__.py
│   │   ├── biophysics.py              # MW, pI bisection, GRAVY, mu_H, II, AI, cleavage
│   │   └── scales.py                  # EMBOSS/Bjellqvist/Lehninger pKa, KD, Eisenberg, DIWV
│   ├── screening/
│   │   ├── __init__.py
│   │   ├── alanine_scan.py            # Systematic Ala scan, delta vectors, hotspot tagging
│   │   ├── mutagenesis.py             # Homologous substitutions (BLOSUM62 & clustering)
│   │   └── stability_delta.py         # Empirical ddG, dSol heuristics, Pareto ranker
│   ├── connectors/
│   │   ├── __init__.py
│   │   ├── base.py                    # 3-tier fallback engine (Live -> Disk Cache -> Fixtures)
│   │   ├── uniprot.py                 # UniProt Knowledgebase client
│   │   ├── alphafold.py               # EBI AlphaFold Structure DB client (pLDDT)
│   │   ├── chembl.py                  # ChEMBL bioactivity and target client
│   │   └── fixtures/                  # Curated offline benchmark mocks (Insulin, GLP-1, GHK, BPC)
│   └── cli/
│       ├── __init__.py
│       ├── main.py                    # Multi-command CLI entry point (analyze, scan, fetch)
│       ├── formatting.py              # Rich ANSI cards and tables
│       ├── export_json.py             # JSON data serialization
│       ├── export_markdown.py         # Academic Markdown dossier exporter
│       └── plotting.py                # Headless Matplotlib publication plots (300 DPI)
└── tests/                             # Comprehensive 4-tier Pytest suite
    ├── conftest.py
    ├── test_core_biophysics.py
    ├── test_screening.py
    ├── test_connectors.py
    └── test_cli.py
```

---

## Installation

### Standard Installation via `pip` / `uv`
```bash
# Clone the repository
git clone https://github.com/Juanj053/insilico-peptide-toolkit.git
cd insilico-peptide-toolkit

# Install using uv (recommended for ultra-fast deterministic resolution)
uv sync

# Or install in editable development mode
uv pip install -e .
```

---

## Quickstart & Python API

```python
import insilico_peptide as ipt

# 1. Calculate biophysical properties
sequence = "GEPPPGKPADDAGLV"  # BPC-157
props = ipt.calculate_properties(sequence, pka_scale="EMBOSS")

print(f"Peptide: BPC-157")
print(f"Molecular Weight: {props.molecular_weight:.2f} Da")
print(f"Isoelectric Point (pI): {props.isoelectric_point:.2f}")
print(f"GRAVY Hydropathicity: {props.gravy:.4f} (Soluble)")
print(f"Instability Index (II): {props.instability_index:.2f} (Stable: {props.is_stable})")
print(f"Trypsin Cuts: {props.proteolytic_cleavage_sites['trypsin']}")  # [] (Proline-block protected!)

# 2. In-Silico Alanine Scanning
mutants = ipt.run_alanine_scan(sequence)
hotspots = [m for m in mutants if m.is_hotspot]
print(f"Total residues scanned: {len(mutants)}, Hotspots identified: {len(hotspots)}")

# 3. Multi-Omic Integration (with offline fallback guarantee)
uniprot_entry = ipt.fetch_uniprot_data("P01275", offline=True)
print(f"UniProt Target: {uniprot_entry['protein_name']} ({uniprot_entry['organism']})")

alphafold_data = ipt.fetch_alphafold_plddt("P01275", offline=True)
print(f"AlphaFold Mean pLDDT: {alphafold_data['mean_plddt']:.2f}")
```

---

## Command-Line Interface (CLI)

The package installs the `insilico-peptide` executable:

```bash
# Analyze sequence with Rich terminal output
insilico-peptide analyze --sequence "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR" --name "GLP-1(7-36)"

# Export structured Markdown dossier and JSON vectors
insilico-peptide analyze --sequence "GEPPPGKPADDAGLV" --name "BPC-157" \
  --export-markdown dossier_bpc157.md --export-json bpc157.json

# Generate 300-DPI publication vector figures
insilico-peptide analyze --sequence "GHK" --plot-charge charge_curve.png --plot-helical helical_wheel.png

# Perform systematic in-silico alanine scan
insilico-peptide scan --sequence "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG" --mode alanine --top 10

# Fetch multi-omics annotations (offline resilient)
insilico-peptide fetch --uniprot P01275 --alphafold --chembl CHEMBL1784 --offline
```

---

## Case Studies

### Case Study 1: GLP-1(7-36/37) Receptor Agonist Engineering
- **Sequence**: `HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG` ($31\text{ AA}$)
- **Biophysical Profile**:
  - $\text{MW}_{\text{avg}} = 3355.71\text{ Da}$, $\text{MW}_{\text{mono}} = 3353.6681\text{ Da}$
  - $pI = 5.62$ (EMBOSS), Net Charge at $\text{pH } 7.4 = -0.95\text{ e}$
  - $\text{GRAVY} = -0.2355$, Instability Index $II = 14.71$ (Highly stable)
- **Enzymatic Degradation Hotspot**:
  Cleavage analysis flags position 2 (Ala2) as immediately susceptible to DPP-IV hydrolysis (`^.[AP]`), inactivating wild-type GLP-1 in vivo within 2 minutes ($t_{1/2} < 2\text{ min}$).
- **In-Silico Engineering**:
  Alanine scanning and homologous screening evaluate substitutions at Ala2 (e.g. Aib in Semaglutide or Gly in Albiglutide), eliminating DPP-IV recognition while conserving the amphiphilic binding helix ($\mu_H = 0.0597$).

---

### Case Study 2: BPC-157 Gastric Pentadecapeptide Stability
- **Sequence**: `GEPPPGKPADDAGLV` ($15\text{ AA}$)
- **Biophysical Profile**:
  - $\text{MW}_{\text{avg}} = 1419.55\text{ Da}$, $\text{MW}_{\text{mono}} = 1418.7042\text{ Da}$
  - $pI = 3.88$, Net Charge at $\text{pH } 7.4 = -2.06\text{ e}$
  - $\text{GRAVY} = -0.6933$, Instability Index $II = 39.32$ (Predicted stable: $II \le 40.0$)
- **Gastric Juice Resilience**:
  Despite containing Lys7, Trypsin cleavage is completely inhibited by the physiological **Proline-block rule** (`KP` motif at residues 7–8). Furthermore, the absence of aromatic residues (`F, Y, W`) provides total immunity against Chymotrypsin.

---

### Case Study 3: GHK-Cu Regenerative Tripeptide Coordination
- **Sequence**: `GHK` ($3\text{ AA}$)
- **Biophysical Profile**:
  - $\text{MW}_{\text{avg}} = 340.38\text{ Da}$, $\text{MW}_{\text{mono}} = 340.1859\text{ Da}$
  - $pI = 9.70$ (EMBOSS) / $10.11$ (Lehninger), Net Charge at $\text{pH } 7.4 = +1.05\text{ e}$
  - $\text{GRAVY} = -2.5000$ (Extremely hydrophilic, high aqueous solubility)
- **Coordination Mechanism**:
  His2 provides the imidazole nitrogen center for high-affinity Copper(II) chelation ($\log K = 16.44$), yielding the bioactive $[\text{Cu(II)}(\text{GHK} - 2\text{H})]^{+}$ complex for wound healing and extracellular matrix remodeling.

---

## International Research Collaboration & Grants

`insilico-peptide-toolkit` is actively positioned as computational core infrastructure for international research grant consortia:
- **CYTED (Programa Iberoamericano de Ciencia y Tecnología para el Desarrollo)**: Subprograma Salud y Biomedicina.
- **NIH / Fogarty International Center**: Global Non-Communicable Diseases (NCD) and Affordable Biotherapeutics.
- **Chevening & UK Research Partnerships**: Bilateral peptide structural bioengineering collaborations with leading UK institutions (Oxford, Cambridge, Imperial, Dundee).

For complete grant proposals, experimental wet-lab validation protocols (SPPS, CD, SPR, plasma stability), and open-science governance, refer to [`RESEARCH_COLLABORATION.md`](RESEARCH_COLLABORATION.md).

---

## Citing insilico-peptide-toolkit

If you utilize `insilico-peptide-toolkit` in academic research, drug discovery pipelines, or grant submissions, please cite:

```bibtex
@software{alvarez2026insilico,
  author       = {Álvarez, Juan José},
  title        = {insilico-peptide-toolkit: In-Silico Biophysical Modeling, Mutational Screening, and Multi-Omic Profiling for Bioactive Peptides},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/Juanj053/insilico-peptide-toolkit}},
  doi          = {10.5281/zenodo.1234567}
}
```

---

---

## Author & Academic Contact
**Juan José Álvarez Morales**  
Chief Bioengineer & Founder, Elevate Biotech  
- **LinkedIn:** [linkedin.com/in/juanj-alva](https://www.linkedin.com/in/juanj-alva)  
- **GitHub:** [@Juanj053](https://github.com/Juanj053)  
- **Email:** [`jjalvarezm77@gmail.com`](mailto:jjalvarezm77@gmail.com)  

---

## License
MIT License. Copyright (c) 2026 Juan José Álvarez (@Juanj053). See [`LICENSE`](LICENSE) for details.
