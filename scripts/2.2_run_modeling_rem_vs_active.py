#!/usr/bin/env python
"""
Command-line entry point for remission vs active modeling

This script provides functions for:
- Parsing CLI flags (`parse_arguments`)
- Running remission vs active analyses and aggregations (`main`)
"""

# Imports
import argparse
import sys
from typing import List, Optional, Dict

import pandas as pd

# Specific imports from ccf_medication modules
from ccf_medication.models.calls import (
    run_me_model_analysis,
)

# Constants
from ccf_medication.constants.pathing import (
    REMISSION_TX_DATAFRAME_PATH,
    REMISSION_PX_DATAFRAME_PATH,
    ACTIVE_TX_DATAFRAME_PATH,
    ACTIVE_PX_DATAFRAME_PATH,
    TX_GENE_COLS,
    PX_PROTEIN_COLS,
    REM_VS_ACT_PX_DIR,
    REM_VS_ACT_TX_DIR,
    AGG_REM_VS_ACT_TX_RESULTS_PATH,
    AGG_REM_VS_ACT_PX_RESULTS_PATH,
)

from ccf_medication.constants.tables import DRUG_FAMILY_COL

from ccf_medication.constants.thresholds import (
    MIN_PATIENTS_PER_SEVERITY_GROUP,
    MIN_PATIENTS_PER_SEVERITY_GROUP,
)

# Functions
def parse_arguments(argv: List[str]) -> argparse.Namespace:
    """
    Parse command line arguments.

    :param argv: Command line arguments.
    :return: Parsed command line arguments.
    """
    
    parser = argparse.ArgumentParser(
        description="Run disease severity analyses for proteomics and transcriptomics and aggregate results."
    )
    parser.add_argument(
        "--rerun-proteomics",
        action="store_true",
        default=True,
        help="Rerun proteomics analysis before aggregation (default: True)",
    )
    parser.add_argument(
        "--no-rerun-proteomics",
        dest="rerun_proteomics",
        action="store_false",
        help="Skip re-running proteomics analysis",
    )
    parser.add_argument(
        "--rerun-transcriptomics",
        action="store_true",
        default=True,
        help="Rerun transcriptomics analysis before aggregation (default: True)",
    )
    parser.add_argument(
        "--no-rerun-transcriptomics",
        dest="rerun_transcriptomics",
        action="store_false",
        help="Skip re-running transcriptomics analysis",
    )

    parser.add_argument(
        "--min-patients",
        type=int,
        default=MIN_PATIENTS_PER_SEVERITY_GROUP,
        help=f"Override minimum patients per group (for both remission and active groups) (default: {MIN_PATIENTS_PER_SEVERITY_GROUP})",
    )

    parser.add_argument(
        "--exclude-values",
        default={"drug_family": ["Combination Therapy"]},
        help='Values to exclude from analysis (default: {"drug_family": ["Combination Therapy"]})',
    )
    return parser.parse_args(argv)


def _map_to_binary_severity(
    df: pd.DataFrame,
    *,
    source_col: str = "endo_category",
    output_col: str = "disease_severity",
    severity_map: Optional[Dict[str, str]] = {
            "remission": "remission",
            "moderate": "active",
            "severe": "active",
        },
) -> pd.DataFrame:
    """
    Map the severity of the disease to a binary severity.

    :param df: DataFrame to map the severity of the disease to a binary severity.
    :param source_col: Column to map the severity of the disease to a binary severity.
    :param output_col: Column to output the binary severity.
    :param severity_map: Map of the severity of the disease to a binary severity.
    :return: DataFrame with the binary severity.
    """

    if severity_map is None:
        raise ValueError("severity_map must be provided")

    output = df.copy()
    if source_col not in output.columns:
        raise ValueError(f"Severity source column '{source_col}' not found in data")

    output[output_col] = output[source_col].map(severity_map)
    output = output.dropna(subset=[output_col])
    return output

def main(argv: List[str]) -> int:
    """
    Main function to run remission vs active modeling analyses and aggregations.

    :param argv: Command line arguments.
    :return: 0 if successful.
    """
    
    args = parse_arguments(argv)

    print("Loading input datasets…")
    remission_tx_data = pd.read_parquet(REMISSION_TX_DATAFRAME_PATH)
    remission_px_data = pd.read_parquet(REMISSION_PX_DATAFRAME_PATH)
    active_tx_data = pd.read_parquet(ACTIVE_TX_DATAFRAME_PATH)
    active_px_data = pd.read_parquet(ACTIVE_PX_DATAFRAME_PATH)

    tx_data = pd.concat([remission_tx_data, active_tx_data])
    px_data = pd.concat([remission_px_data, active_px_data])

    tx_data = _map_to_binary_severity(tx_data, source_col="endo_category", output_col="disease_severity")
    px_data = _map_to_binary_severity(px_data, source_col="endo_category", output_col="disease_severity")

    # Load feature columns
    tx_gene_cols = open(TX_GENE_COLS, "r").read().split("\n")
    px_gene_cols = open(PX_PROTEIN_COLS, "r").read().split("\n")

    # Run analyses
    if args.rerun_proteomics:
        print("Running proteomics analysis…")

        agg_px_results = run_me_model_analysis(
            data=px_data,
            gene_cols=px_gene_cols,
            min_patients=args.min_patients,
            omics_type="proteomics",
            cols_to_split_on={"diagnosis": ["uc", "cd"], DRUG_FAMILY_COL: None},
            primary_effect_col="disease_severity",
            reference_groups={"disease_severity": "remission"},
            additional_fixed_effects=[],
            random_effects=["patient_id"],
            exclude_values=args.exclude_values,
            save_path=REM_VS_ACT_PX_DIR,
        )
        pre_edit_path = AGG_REM_VS_ACT_PX_RESULTS_PATH.split('.parquet')[0] + "_pre_edit.parquet"
        agg_px_results.to_parquet(pre_edit_path)

        agg_px_results.drop(columns=['reference_values', 'non_reference_values'], inplace=True)
        agg_px_results.rename(columns={'num_patients_non_reference': 'n_patients_in_active_group',
                                'num_samples_non_reference': 'n_samples_in_active_group',
                                f'num_patients_remission': 'n_patients_in_rem_group',
                                f'num_samples_remission': 'n_samples_in_rem_group'}, inplace=True)
        agg_px_results.to_parquet(AGG_REM_VS_ACT_PX_RESULTS_PATH)


    if args.rerun_transcriptomics:
        print("Running transcriptomics analysis…")
        agg_tx_results = run_me_model_analysis(
            data=tx_data,
            gene_cols=tx_gene_cols,
            min_patients=args.min_patients,
            omics_type="transcriptomics",
            cols_to_split_on={"diagnosis": ["uc", "cd"], "simple_tissue": ["colon", "small_intestine"], DRUG_FAMILY_COL: None},
            primary_effect_col="disease_severity",
            reference_groups={"disease_severity": "remission"},
            additional_fixed_effects=[],
            random_effects=["patient_id"],
            exclude_values=args.exclude_values,
            save_path=REM_VS_ACT_TX_DIR,
        )
        pre_edit_path = AGG_REM_VS_ACT_TX_RESULTS_PATH.split('.parquet')[0] + "_pre_edit.parquet"
        agg_tx_results.to_parquet(pre_edit_path)

        agg_tx_results.drop(columns=['reference_values', 'non_reference_values'], inplace=True)
        agg_tx_results.rename(columns={'num_patients_non_reference': 'n_patients_in_active_group',
                                'num_samples_non_reference': 'n_samples_in_active_group',
                                f'num_patients_remission': 'n_patients_in_rem_group',
                                f'num_samples_remission': 'n_samples_in_rem_group'}, inplace=True)
        agg_tx_results.to_parquet(AGG_REM_VS_ACT_TX_RESULTS_PATH)


    print("All tasks completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))