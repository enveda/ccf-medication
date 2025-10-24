"""
Data loading, preparation, and formatting helpers

This module provides functions for:
- Loading aggregated results and formatting per-analysis outputs (`load_aggregated_results`, `load_and_format_results`, `load_and_format_all_results`)
- Processing RNA-seq counts (`process_rnaseq_counts`)
- Merging datasets and detecting omics columns (`merge_datasets`, `get_proteomics_and_transcriptomics_columns`)
- Reformatting 2024 metadata to match original schema (`reformat_2024_metadata`)
- Formatting samples with disease severity (`format_samples_with_disease_severity`)
"""

import logging
from typing import List, Tuple

import pandas as pd

from ccf_medication.constants.system import (
    LOGGER_NAME,
)
from ccf_medication.constants.tables import (
    METADATA_2024_COLS_RENAME_MAP,
    PVAL_COL_MED_VS_NO_MED,
    FC_COL_MED_VS_NO_MED,
    PVAL_COL_REM_VS_ACT,
    FC_COL_REM_VS_ACT,
)
from ccf_medication.constants.pathing import (
    AGG_REM_VS_ACT_TX_RESULTS_PATH,
    AGG_REM_VS_ACT_PX_RESULTS_PATH,
    AGG_MED_VS_NO_MED_TX_RESULTS_PATH,
    AGG_MED_VS_NO_MED_PX_RESULTS_PATH,
)

from ccf_medication.constants.thresholds import (
    ADJ_PVAL_THRESH,
    TX_FC_THRESH,
    PX_FC_THRESH,
)

from ccf_medication.utils.analysis import get_response_genes, get_significant_genes

logger = logging.getLogger(LOGGER_NAME)

def load_aggregated_results(file_path: str) -> pd.DataFrame:
    """
    Load aggregated results from a Parquet file.

    :param file_path: Path to the aggregated results Parquet file.

    :return: DataFrame of aggregated results. Removes a deprecated ``effect`` column
        if present for backward compatibility.
    """
    df = pd.read_parquet(file_path)

    # Remove deprecated 'effect' column if present
    if "effect" in df.columns:
        df = df.drop(columns=["effect"])

    return df


def process_rnaseq_counts(
    input_data: pd.DataFrame,
    number_of_samples: int,
    ordered_samples: List,
    sample_id: str = "sample_id",
    min_count_sum_across_samples_dynamic: int = 0,
    min_count_sum_across_samples: int = 10,
    min_counts_in_at_least_one_sample: int = 15,
) -> pd.DataFrame:
    """Process raw RNA-seq counts by filtering low-count and zero-variance features.

    This function does not perform normalization (e.g., DESeq2); normalization is
    assumed to happen upstream.

    :param input_data: Raw counts with genes/features in columns and samples in rows.
    :param number_of_samples: Total number of samples; used when applying a dynamic
        threshold.
    :param ordered_samples: Ordered sample identifiers to assign to the output's
        ``sample_id`` column.
    :param sample_id: Name of the column to store ``ordered_samples`` in.
    :param min_count_sum_across_samples_dynamic: If > 0, dynamic threshold applied as
        ``number_of_samples * min_count_sum_across_samples_dynamic``. When set,
        overrides ``min_count_sum_across_samples``.
    :param min_count_sum_across_samples: Fixed minimum sum across samples to retain a
        gene/feature. Used only when the dynamic threshold is 0.
    :param min_counts_in_at_least_one_sample: Minimum counts required in at least one
        sample to retain a gene/feature.
    :return: DataFrame of filtered counts with an added ``sample_id`` column reflecting
        ``ordered_samples``.
    """
    if min_count_sum_across_samples_dynamic > 0:
        min_count_sum_across_samples = (
            number_of_samples * min_count_sum_across_samples_dynamic
        )  # This is a dynamic threshold based on the number of samples, thus the threshold is multiplied by the number of samples

    counts = input_data.astype(int).reset_index(drop=True)

    logger.info(f"Number of genes before filtering: {counts.shape[1]}")
    logger.info(f"Threshold: {min_count_sum_across_samples}")
    # Remove genes that don't have at least necessary_min_counts in any sample
    counts = counts.loc[:, (counts >= min_counts_in_at_least_one_sample).any()]
    logger.info(
        f"Number of genes after removing those without minimum counts ({min_counts_in_at_least_one_sample}) in at least one sample: {counts.shape[1]}"
    )

    logger.info(f"Original number of genes: {counts.shape[1]}")
    counts = counts.loc[:, counts.var() > 0]
    logger.info(f"Number of genes after removing 0 variance: {counts.shape[1]}")
    counts = counts.loc[:, counts.sum() > min_count_sum_across_samples]
    logger.info(f"Number of genes after removing low sum: {counts.shape[1]}")

    normalized_counts = counts
    normalized_counts = normalized_counts.reset_index(drop=True)

    normalized_counts[sample_id] = ordered_samples
    return normalized_counts


def merge_datasets(
    metadata,
    omics,
    omics_sample_id_col,
    categories_of_interest,
    metadata_columns_of_interest,
):
    """
    Merge omics measurements with metadata and apply standard filters.

    :param metadata: Metadata table.
    :param omics: Omics measurements that include a ``sample_id`` column.
    :param omics_sample_id_col: Column in ``metadata`` that maps to ``omics.sample_id``
        (e.g., ``"transcriptomics"`` when those values map to the omics
        ``sample_id``).
    :param categories_of_interest: Mapping from categorical column name to the
        ordered list of allowed values. Used for filtering and to set categorical order.
    :param metadata_columns_of_interest: Metadata columns to retain in the merge.

    :return: DataFrame containing the merged and filtered omics data with selected
        metadata columns.
    """

    def apply_filters(data, categories_of_interest):
        """Apply standard filters to merged data"""
        # Filter to CD and UC only
        data = data[data["diagnosis"].isin(["cd", "uc"])]

        # Apply category filters
        for category, category_order in categories_of_interest.items():
            data = data[data[category].isin(category_order)]
            # Make categorical
            data[category] = pd.Categorical(
                data[category], categories=category_order, ordered=True
            )

        return data.reset_index(drop=True)

    # Merge transcriptomics with metadata
    omics_data = (
        pd.merge(
            metadata[[f"{omics_sample_id_col}"] + metadata_columns_of_interest],
            omics,
            left_on=f"{omics_sample_id_col}",
            right_on="sample_id",
            how="inner",
        )
        .dropna(subset=[f"{omics_sample_id_col}"])
        .drop(columns=[f"{omics_sample_id_col}"])
    )
    omics_data = apply_filters(omics_data, categories_of_interest)

    return omics_data


def get_proteomics_and_transcriptomics_columns(
    merged_data, metadata_df, transcriptomics_df, proteomics_df
):
    """
    Identify numeric proteomics and transcriptomics columns in the dataset.

    :param merged_data: DataFrame containing merged omics data and metadata.
    :param metadata_df: Metadata table used to exclude metadata columns.
    :param transcriptomics_df: Original transcriptomics data used to seed candidate
        column names.
    :param proteomics_df: Original proteomics data used to seed candidate column names.

    :return: Tuple[List[str], List[str]]: ``(proteomics_cols, transcriptomics_cols)``
        lists of numeric column names excluding metadata columns and ``sample_id``.
    """

    transcriptomics_cols = transcriptomics_df.columns.tolist()
    proteomics_cols = proteomics_df.columns.tolist()

    # remove any metadata columns from the proteomics and transcriptomics dataframes
    proteomics_cols = set(proteomics_cols) - set(metadata_df.columns)
    transcriptomics_cols = set(transcriptomics_cols) - set(metadata_df.columns)

    # remove sample_id from the proteomics and transcriptomics dataframes
    proteomics_cols = set(proteomics_cols) - set(["sample_id"])
    transcriptomics_cols = set(transcriptomics_cols) - set(["sample_id"])

    # remove any columns that are not numeric from the proteomics and transcriptomics dataframes
    proteomics_cols = [
        col
        for col in proteomics_cols
        if pd.api.types.is_numeric_dtype(merged_data[col])
    ]
    transcriptomics_cols = [
        col
        for col in transcriptomics_cols
        if pd.api.types.is_numeric_dtype(merged_data[col])
    ]

    print(
        f"Using {len(proteomics_cols)} proteomics columns and {len(transcriptomics_cols)} transcriptomics columns"
    )

    return proteomics_cols, transcriptomics_cols


def reformat_2024_metadata(
    metadata_2024, metadata_2024_cols_rename_map=METADATA_2024_COLS_RENAME_MAP
):
    """
    Reformat the 2024 metadata to match the original schema.

    :param metadata_2024: Raw 2024 metadata table to be reformatted.
    :param metadata_2024_cols_rename_map: Mapping from raw column names to
        standardized names.

    :return: DataFrame with standardized column names and normalized values for
        medication, endoscopic category, and tissue labels.
    """
    metadata_2024.rename(columns=metadata_2024_cols_rename_map, inplace=True)

    # This makes sure "Combination Therapy" is the name for any combination of drugs like the original metadata
    metadata_2024["medication"] = metadata_2024["raw_medication"].apply(
        lambda x: "Combination Therapy" if ";" in x else x
    )

    # The original endo_category is all lowercase (e.g. "remission", "mild", "moderate", "severe")
    metadata_2024["endo_category"] = metadata_2024["endo_category"].str.lower()

    # This makes sure "small intestine" is "small_intestine"
    metadata_2024["simple_tissue"] = metadata_2024["simple_tissue"].str.replace(
        " ", "_"
    )

    # if there are any str "nan" values in the medication column, replace them with pd.NA
    metadata_2024 = metadata_2024.replace("nan", pd.NA)

    metadata_2024["birth_year"] = metadata_2024["birth_year"].astype(int)

    return metadata_2024

def load_and_format_results(file_path: str, is_rem_vs_act: bool = False) -> pd.DataFrame:
    """
    Load and format the results.

    :param file_path: The path to the results file.
    :param is_rem_vs_act: Whether the results are for the remission vs active (true) 
                            or the medication vs no medication analysis (false).


    :return: The formatted results.
    """
    results = load_aggregated_results(file_path)
    if is_rem_vs_act:
        cols_to_rename = {"adjusted_p_value": PVAL_COL_REM_VS_ACT, 
                         "coefficient": FC_COL_REM_VS_ACT,
                         "p_value": "pval_rem_vs_act",
                         "standard_error": "std_error_rem_vs_act"}
        
    else:
        cols_to_rename = {"adjusted_p_value": PVAL_COL_MED_VS_NO_MED, 
                          "coefficient": FC_COL_MED_VS_NO_MED,
                          "p_value": "pval_med_vs_no_med",
                          "standard_error": "std_error_med_vs_no_med"}
    

    results = results.rename(columns=cols_to_rename)
    return results

def load_and_format_all_results(*, 
                             rem_vs_act_tx_path: str = AGG_REM_VS_ACT_TX_RESULTS_PATH, 
                             rem_vs_act_px_path: str = AGG_REM_VS_ACT_PX_RESULTS_PATH, 
                             med_vs_no_med_tx_path: str = AGG_MED_VS_NO_MED_TX_RESULTS_PATH, 
                             med_vs_no_med_px_path: str = AGG_MED_VS_NO_MED_PX_RESULTS_PATH,
                             pval_thresh_rem_vs_act: float = ADJ_PVAL_THRESH,
                             tx_fc_thresh_rem_vs_act: float = TX_FC_THRESH,
                             px_fc_thresh_rem_vs_act: float = PX_FC_THRESH,
                             pval_thresh_med_vs_no_med: float = ADJ_PVAL_THRESH,
                             tx_fc_thresh_med_vs_no_med: float = TX_FC_THRESH, 
                             px_fc_thresh_med_vs_no_med: float = PX_FC_THRESH, 
                             fc_col_rem_vs_act: str = FC_COL_REM_VS_ACT,
                             pval_col_rem_vs_act: str = PVAL_COL_REM_VS_ACT,
                             fc_col_med_vs_no_med: str = FC_COL_MED_VS_NO_MED,
                             pval_col_med_vs_no_med: str = PVAL_COL_MED_VS_NO_MED) -> Tuple[Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame], Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]]:
    """
    Load and format the all results.

    :param rem_vs_act_tx_path: The path to the remission vs active transcriptomics results.
    :param rem_vs_act_px_path: The path to the remission vs active proteomics results.
    :param med_vs_no_med_tx_path: The path to the medication vs no medication transcriptomics results.
    :param med_vs_no_med_px_path: The path to the medication vs no medication proteomics results.
    :param pval_thresh_rem_vs_act: The p-value threshold for the remission vs active analysis.
    :param tx_fc_thresh_rem_vs_act: The fold change threshold for the remission vs active analysis.
    :param px_fc_thresh_rem_vs_act: The fold change threshold for the remission vs active analysis.
    :param pval_thresh_med_vs_no_med: The p-value threshold for the medication vs no medication analysis.
    :param tx_fc_thresh_med_vs_no_med: The fold change threshold for the medication vs no medication analysis.
    :param px_fc_thresh_med_vs_no_med: The fold change threshold for the medication vs no medication analysis.
    
    :return: (tx_response, px_response, tx_signif_rem_vs_act, px_signif_rem_vs_act, tx_signif_med_vs_no_med, px_signif_med_vs_no_med), 
            (tx_rem_vs_act, px_rem_vs_act, tx_med_vs_no_med, px_med_vs_no_med)
    """

    # load and format the results
    tx_rem_vs_act = load_and_format_results(rem_vs_act_tx_path, is_rem_vs_act=True)
    px_rem_vs_act = load_and_format_results(rem_vs_act_px_path, is_rem_vs_act=True)
    tx_med_vs_no_med = load_and_format_results(med_vs_no_med_tx_path, is_rem_vs_act=False)
    px_med_vs_no_med = load_and_format_results(med_vs_no_med_px_path, is_rem_vs_act=False)


    tx_signif_rem_vs_act = get_significant_genes(tx_rem_vs_act, 
                                                 pval_threshold=pval_thresh_rem_vs_act, 
                                                 fc_threshold=tx_fc_thresh_rem_vs_act,
                                                 pval_column=pval_col_rem_vs_act,
                                                 fc_column=fc_col_rem_vs_act)
    px_signif_rem_vs_act = get_significant_genes(px_rem_vs_act, 
                                                 pval_threshold=pval_thresh_rem_vs_act, 
                                                 fc_threshold=px_fc_thresh_rem_vs_act,
                                                 pval_column=pval_col_rem_vs_act,
                                                 fc_column=fc_col_rem_vs_act)
    tx_signif_med_vs_no_med = get_significant_genes(tx_med_vs_no_med, 
                                                    pval_threshold=pval_thresh_med_vs_no_med, 
                                                    fc_threshold=tx_fc_thresh_med_vs_no_med,
                                                    pval_column=pval_col_med_vs_no_med,
                                                    fc_column=fc_col_med_vs_no_med)
    px_signif_med_vs_no_med = get_significant_genes(px_med_vs_no_med, 
                                                    pval_threshold=pval_thresh_med_vs_no_med, 
                                                    fc_threshold=px_fc_thresh_med_vs_no_med,
                                                    pval_column=pval_col_med_vs_no_med,
                                                    fc_column=fc_col_med_vs_no_med)

    tx_response = get_response_genes(tx_med_vs_no_med, 
                        tx_rem_vs_act, 
                        omics_type='tx',
                        pval_thresh_med_vs_no_med=pval_thresh_med_vs_no_med,
                        fc_thresh_med_vs_no_med=tx_fc_thresh_med_vs_no_med,
                        pval_thresh_rem_vs_act=pval_thresh_rem_vs_act,
                        fc_thresh_rem_vs_act=tx_fc_thresh_rem_vs_act,
                        pval_col_med_vs_no_med=pval_col_med_vs_no_med,
                        fc_col_med_vs_no_med=fc_col_med_vs_no_med,
                        pval_col_rem_vs_act=pval_col_rem_vs_act,
                        fc_col_rem_vs_act=fc_col_rem_vs_act,
                        )
    if 'effect' in tx_response.columns:
        tx_response.drop(columns=['effect'], inplace=True)

    px_response = get_response_genes(px_med_vs_no_med, 
                    px_rem_vs_act, 
                    omics_type='px',
                    pval_thresh_med_vs_no_med=pval_thresh_med_vs_no_med,
                    fc_thresh_med_vs_no_med=px_fc_thresh_med_vs_no_med,
                    pval_thresh_rem_vs_act=pval_thresh_rem_vs_act,
                    fc_thresh_rem_vs_act=px_fc_thresh_rem_vs_act,
                    pval_col_med_vs_no_med=pval_col_med_vs_no_med,
                    fc_col_med_vs_no_med=fc_col_med_vs_no_med,
                    pval_col_rem_vs_act=pval_col_rem_vs_act,
                    fc_col_rem_vs_act=fc_col_rem_vs_act,
                    )
    if 'effect' in px_response.columns:
        px_response.drop(columns=['effect'], inplace=True)


    return (tx_response, px_response, tx_signif_rem_vs_act, px_signif_rem_vs_act, tx_signif_med_vs_no_med, px_signif_med_vs_no_med), (tx_rem_vs_act, px_rem_vs_act, tx_med_vs_no_med, px_med_vs_no_med)

def format_samples_with_disease_severity(*, tx_samples: pd.DataFrame, px_samples: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Format the samples with disease severity.

    :param tx_samples: The transcriptomics samples.
    :param px_samples: The proteomics samples.
    """
    disease_severity_map = {"remission": 0, "severe": 1, "moderate": 1}

    tx_samples_subset = tx_samples.copy()
    px_samples_subset = px_samples.copy()

    tx_samples_subset['disease_severity'] = tx_samples_subset['endo_category'].map(disease_severity_map)
    px_samples_subset['disease_severity'] = px_samples_subset['endo_category'].map(disease_severity_map)

    tx_samples_subset = tx_samples_subset[tx_samples_subset['simple_tissue'].isin(['colon', 'small_intestine'])]

    tx_samples_subset = tx_samples_subset.dropna(subset=['disease_severity'])
    px_samples_subset = px_samples_subset.dropna(subset=['disease_severity'])

    return tx_samples_subset, px_samples_subset