#!/usr/bin/env python
"""
Command-line entry point for med vs no med modeling

This script provides functions for:
- Parsing CLI flags (`parse_arguments`)
- Running med vs no med analyses and aggregations (`main`)
"""

# Imports
import argparse
import sys
from typing import List

import pandas as pd

# Specific imports from ccf_medication modules
from ccf_medication.models.calls import (
    run_me_model_analysis,
)

# Constants
from ccf_medication.constants.pathing import (
    REMISSION_TX_DATAFRAME_PATH,
    REMISSION_PX_DATAFRAME_PATH,
    TX_GENE_COLS,
    PX_PROTEIN_COLS,
    AGG_MED_VS_NO_MED_TX_RESULTS_PATH,
    AGG_MED_VS_NO_MED_PX_RESULTS_PATH,
    MED_VS_NO_MED_PX_DIR,
    MED_VS_NO_MED_TX_DIR,
)
from ccf_medication.constants.tables import (
    DRUG_FAMILY_COL,
    NO_MEDICATION,
)
from ccf_medication.constants.thresholds import (
    MIN_PATIENTS_PER_DRUG_FAMILY,
)

# Functions
def parse_arguments(argv: List[str]) -> argparse.Namespace:
    """
    Parse command line arguments.

    :param argv: Command line arguments.
    :return: Parsed command line arguments.
    """
    
    parser = argparse.ArgumentParser(
        description="Run modeling analyses for proteomics and transcriptomics and aggregate results."
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
        default=MIN_PATIENTS_PER_DRUG_FAMILY,
        help=f"Override minimum patients per group (default: constant MIN_PATIENTS_PER_DRUG_FAMILY={MIN_PATIENTS_PER_DRUG_FAMILY})",
    )
    parser.add_argument(
        "--exclude-values",
        default={"drug_family": ["Combination Therapy"]},
        help='Values to exclude from analysis (default: {"drug_family": ["Combination Therapy"]})',
    )
    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    """
    Main function to run med vs no med modeling analyses and aggregations.

    :param argv: Command line arguments.
    :return: 0 if successful.
    """
    
    args = parse_arguments(argv)

    # Load datasets
    print("Loading input datasets…")
    tx_rem_data = pd.read_parquet(REMISSION_TX_DATAFRAME_PATH)
    px_rem_data = pd.read_parquet(REMISSION_PX_DATAFRAME_PATH)

    # Load feature columns
    tx_gene_cols = open(TX_GENE_COLS, "r").read().split("\n")
    px_gene_cols = open(PX_PROTEIN_COLS, "r").read().split("\n")

    min_patients = args.min_patients if args.min_patients is not None else MIN_PATIENTS_PER_DRUG_FAMILY

    # # Run analyses
    if args.rerun_proteomics:
        print("Running proteomics analysis…")
        agg_px_results = run_me_model_analysis(
            data=px_rem_data,
            gene_cols=px_gene_cols,
            min_patients=args.min_patients,
            omics_type="proteomics",
            cols_to_split_on={"diagnosis": ["uc", "cd"]},
            primary_effect_col=DRUG_FAMILY_COL,
            reference_groups={DRUG_FAMILY_COL: NO_MEDICATION},
            additional_fixed_effects=[],
            random_effects=["patient_id"],
            exclude_values=args.exclude_values,
            save_path=MED_VS_NO_MED_PX_DIR,
        )
        pre_edit_path = AGG_MED_VS_NO_MED_PX_RESULTS_PATH.split('.parquet')[0] + "_pre_edit.parquet"
        agg_px_results.to_parquet(pre_edit_path)

        agg_px_results.drop(columns=['reference_values'], inplace=True)
        agg_px_results.rename(columns={'non_reference_values': 'drug_family',
                                       'num_patients_non_reference': 'n_patients_in_med_group',
                                       'num_samples_non_reference': 'n_samples_in_med_group',
                                       f'num_patients_{NO_MEDICATION}': 'n_patients_in_no_med_group',
                                       f'num_samples_{NO_MEDICATION}': 'n_samples_in_no_med_group'}, inplace=True)
        agg_px_results.to_parquet(AGG_MED_VS_NO_MED_PX_RESULTS_PATH)

    if args.rerun_transcriptomics:
        print("Running transcriptomics analysis…")
        agg_tx_results = run_me_model_analysis(
            data=tx_rem_data,
            gene_cols=tx_gene_cols,
            min_patients=args.min_patients,
            omics_type="transcriptomics",
            cols_to_split_on={"diagnosis": ["uc", "cd"], "simple_tissue": ["colon", "small_intestine"]},
            primary_effect_col=DRUG_FAMILY_COL,
            reference_groups={DRUG_FAMILY_COL: NO_MEDICATION},
            additional_fixed_effects=[],
            random_effects=["patient_id"],
            exclude_values=args.exclude_values,
            save_path=MED_VS_NO_MED_TX_DIR,
        )
        pre_edit_path = AGG_MED_VS_NO_MED_TX_RESULTS_PATH.split('.parquet')[0] + "_pre_edit.parquet"
        agg_tx_results.to_parquet(pre_edit_path)

        agg_tx_results.drop(columns=['reference_values'], inplace=True)
        agg_tx_results.rename(columns={'non_reference_values': 'drug_family',
                                       'num_patients_non_reference': 'n_patients_in_med_group',
                                       'num_samples_non_reference': 'n_samples_in_med_group',
                                       f'num_patients_{NO_MEDICATION}': 'n_patients_in_no_med_group',
                                       f'num_samples_{NO_MEDICATION}': 'n_samples_in_no_med_group'}, inplace=True)
        agg_tx_results.to_parquet(AGG_MED_VS_NO_MED_TX_RESULTS_PATH)

    print("All tasks completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


