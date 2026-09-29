"""
insilico-peptide-toolkit: Biophysical modeling, in-silico mutational screening,
and multi-omic profiling toolkit for bioactive peptides.

Published by Juan José Álvarez (@Juanj053) for Elevate Biotech research collaboration.
"""

from __future__ import annotations

from insilico_peptide.connectors.alphafold import fetch_alphafold_plddt
from insilico_peptide.connectors.chembl import fetch_chembl_bioactivity
from insilico_peptide.connectors.uniprot import fetch_uniprot_data
from insilico_peptide.core.biophysics import (
    PeptideProperties,
    calculate_aliphatic_index,
    calculate_average_mw,
    calculate_gravy,
    calculate_hydrophobic_moment,
    calculate_instability_index,
    calculate_isoelectric_point,
    calculate_molecular_weight,
    calculate_monoisotopic_mw,
    calculate_net_charge,
    calculate_properties,
    calculate_sliding_hydropathy,
    predict_proteolytic_cleavage,
    titrate_charge,
    validate_sequence,
)
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

__version__ = "0.1.0"
__author__ = "Juan José Álvarez (@Juanj053)"

__all__ = [
    "__version__",
    "__author__",
    "validate_sequence",
    "calculate_properties",
    "calculate_molecular_weight",
    "calculate_average_mw",
    "calculate_monoisotopic_mw",
    "calculate_net_charge",
    "titrate_charge",
    "calculate_isoelectric_point",
    "calculate_gravy",
    "calculate_sliding_hydropathy",
    "calculate_hydrophobic_moment",
    "calculate_instability_index",
    "calculate_aliphatic_index",
    "predict_proteolytic_cleavage",
    "PeptideProperties",
    "run_alanine_scan",
    "get_conservative_substitutions",
    "run_homologous_scan",
    "run_site_mutagenesis",
    "MutantProfile",
    "predict_ddg",
    "predict_solubility_delta",
    "rank_pareto",
    "fetch_uniprot_data",
    "fetch_alphafold_plddt",
    "fetch_chembl_bioactivity",
]
