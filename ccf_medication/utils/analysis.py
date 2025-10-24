"""
Functions used to analyze the proteomics and transcriptomics data

This module provides functions for:
- Getting the response genes (special intersection between the significant results for 'rem vs act' and 'med vs no med')
- Getting the significant genes (based on p-value and fold change thresholds)
- Saving genes to organized directory structures
- Getting the top N significant genes
- Getting the number of biomarkers
"""


from typing import Iterable, Mapping, Optional, Union

import pandas as pd
import re

from ccf_medication.constants.tables import (
    FC_COL_MED_VS_NO_MED,
    FC_COL_REM_VS_ACT,
    FOLD_CHANGE_COLUMN,
    P_VALUE_COLUMN,
)
from ccf_medication.constants.thresholds import ADJ_PVAL_THRESH
from ccf_medication.utils.file_management import GeneFileSaver


def get_response_genes(
    med_vs_no_med_results: pd.DataFrame,
    rem_vs_act_results: pd.DataFrame,
    omics_type: str,
    fc_col_rem_vs_act: str = FC_COL_REM_VS_ACT,
    fc_col_med_vs_no_med: str = FC_COL_MED_VS_NO_MED,
    pval_col_rem_vs_act: str = None,
    pval_col_med_vs_no_med: str = None,
    pval_thresh_med_vs_no_med: float = None,
    pval_thresh_rem_vs_act: float = None,
    fc_thresh_med_vs_no_med: float = None,
    fc_thresh_rem_vs_act: float = None,
) -> pd.DataFrame:
    """
    Get responder genes.

    :param med_vs_no_med_results: DataFrame of med vs no med results from the mixed effects model.
    :param rem_vs_act_results: DataFrame of rem vs act results from the mixed effects model.
    :param omics_type: Type of omics data.
    :param fc_col_rem_vs_act: Fold change column for rem vs act results.
    :param fc_col_med_vs_no_med: Fold change column for med vs no med results.
    :param pval_col_rem_vs_act: P-value column for rem vs act results.
    :param pval_col_med_vs_no_med: P-value column for med vs no med results.
    :param pval_thresh_med_vs_no_med: P-value threshold for med vs no med results.
    :param pval_thresh_rem_vs_act: P-value threshold for rem vs act results.
    :param fc_thresh_med_vs_no_med: Threshold for fold change for med vs no med results.
    :param fc_thresh_rem_vs_act: Threshold for fold change for rem vs act results.
    :return: DataFrame of responder genes.
    """
    if pval_thresh_med_vs_no_med is None:
        raise ValueError("pval_thresh_med_vs_no_med is required")
    if fc_thresh_med_vs_no_med is None:
        raise ValueError("fc_thresh_med_vs_no_med is required")
    if fc_thresh_rem_vs_act is None:
        raise ValueError("fc_thresh_rem_vs_act is required")
    if pval_thresh_rem_vs_act is None:
        raise ValueError("pval_thresh_rem_vs_act is required")

    signif_med_vs_no_med = get_significant_genes(
        med_vs_no_med_results, pval_thresh_med_vs_no_med, fc_thresh_med_vs_no_med, pval_col_med_vs_no_med, fc_col_med_vs_no_med
    )
    signif_rem_vs_act = get_significant_genes(
        rem_vs_act_results,
        pval_thresh_rem_vs_act,
        fc_thresh_rem_vs_act,
        pval_col_rem_vs_act,
        fc_col_rem_vs_act,
    )

    if omics_type.lower() in ["tx", "transcriptomics"]:
        intersection_cols = ["diagnosis", "simple_tissue", "drug_family", "feature"]
    elif omics_type.lower() in ["px", "proteomics"]:
        intersection_cols = ["diagnosis", "drug_family", "feature"]
    else:
        raise ValueError(
            f"Invalid omics type: {omics_type}. Valid options are 'tx' or 'transcriptomics' or 'px' or 'proteomics'."
        )

    intersection = signif_rem_vs_act.merge(
        signif_med_vs_no_med,
        on=intersection_cols,  # or left_on=[...], right_on=[...] if names differ
        how="inner",  # 'left', 'right', 'outer' also work
        validate="1:1",  # try '1:1', '1:m', 'm:1' to catch duplicates
        suffixes=("_x", "_y"),
    )

    return intersection[
        intersection[fc_col_rem_vs_act] * intersection[fc_col_med_vs_no_med] < 0
    ]


def get_significant_genes(
    results_df: pd.DataFrame,
    pval_threshold: float = ADJ_PVAL_THRESH,
    fc_threshold: float = None,
    pval_column: str = P_VALUE_COLUMN,
    fc_column: str = FOLD_CHANGE_COLUMN,
) -> pd.DataFrame:
    """
    Filter results dataframe for significant genes based on p-value and fold change thresholds.

    :param results_df: pd.DataFrame - The results dataframe to filter
    :param pval_threshold: float - The p-value threshold
    :param fc_threshold: float - The fold change threshold
    :param pval_column: str - The column name for the p-values
    :param fc_column: str - The column name for the fold changes
    :return: pd.DataFrame - The filtered results dataframe
    """
    return results_df[
        (results_df[pval_column] < pval_threshold)
        & (results_df[fc_column].abs() > fc_threshold)
    ]


def _validate_sort_direction(direction: Union[bool, str]) -> bool:
    """
    Convert sort direction to boolean.

    :param direction: Union[bool, str] - The direction to convert
    :return: bool - The converted direction
    """
    if isinstance(direction, bool):
        return direction

    direction_str = str(direction).strip().lower()
    if direction_str in {"asc", "ascending", "true", "1"}:
        return True
    if direction_str in {"desc", "descending", "false", "0"}:
        return False

    raise ValueError(
        f"Invalid sort direction: {direction!r} (use True/False or 'asc'/'desc')"
    )


def top_N_genes(
    df: pd.DataFrame,
    *,
    subpop_cols: Iterable[str],
    drug_family_col: str,
    n: int,
    sort_by: Mapping[str, Union[bool, str]],
    gene_col: str = "feature",
    keep_cols: Optional[Iterable[str]] = None,
    dropna_subpops: bool = False,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compute the top-N features per (drug family, subpopulation) group.

    :param df: pd.DataFrame - The dataframe to filter
    :param subpop_cols: Iterable[str] - The columns to group by
    :param drug_family_col: str - The column name for the drug family
    :param n: int - The number of features to keep
    :param sort_by: Mapping[str, Union[bool, str]] - The columns to sort by
    :param gene_col: str - The column name for the gene
    :param keep_cols: Optional[Iterable[str]] - The columns to keep
    :param dropna_subpops: bool - Whether to drop NA in subpopulations
    :param save_path: Optional[str] - The path to save the results
    :return: pd.DataFrame - The top-N features
    """
    subpop_cols = list(subpop_cols)
    required = set(subpop_cols) | {gene_col, drug_family_col}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    split_cols = [drug_family_col] + subpop_cols
    dfx = df.copy()

    if dropna_subpops:
        dfx = dfx.dropna(subset=split_cols)

    # Validate and process sort columns
    sort_cols = []
    sort_asc = []

    for col, direction in sort_by.items():
        if col not in dfx.columns:
            raise KeyError(f"Sort column '{col}' not found in DataFrame.")
        sort_cols.append(col)
        sort_asc.append(_validate_sort_direction(direction))

    # Sort by subpopulation keys then requested sort columns
    order_by = split_cols + sort_cols
    asc_flags = [True] * len(split_cols) + sort_asc
    dfx = dfx.sort_values(order_by, ascending=asc_flags, na_position="last")

    # Get top-N per combined subpopulation
    output = dfx.groupby(split_cols, sort=False, group_keys=False).head(n).copy()
    output["rank"] = output.groupby(split_cols, sort=False).cumcount() + 1

    # Organize columns
    extras = list(keep_cols) if keep_cols else []
    present_sort_cols = [c for c in sort_cols if c in output.columns and c != gene_col]
    ordered_cols = (
        split_cols
        + ["rank", gene_col]
        + [c for c in present_sort_cols if c != gene_col]
        + [c for c in extras if c in output.columns]
    )

    result = output[ordered_cols].reset_index(drop=True)

    # Save to files if path provided
    if save_path:
        saver = GeneFileSaver(save_path)
        saver.save_genes(
            result,
            dir_level_cols=subpop_cols,
            file_level_col=drug_family_col,
            col_to_save=gene_col,
            file_suffix=f"top_{n}_genes",
        )

    return result


def _get_bare_features(df: pd.DataFrame, feat_col: str = 'feature', output_feat_col: str = 'bare_feature'):
    """
    Get the bare features from the dataframe. These are the features without the Olink panel information.

    :param df: DataFrame containing the features.
    :param feat_col: Column name for the features.
    :param output_feat_col: Column name for the output features.
    :return: DataFrame with the bare features.
    """
    panels = ['Cardiometabolic', 'Inflammation', 'Neurology', 'Oncology']
    subpanels = ['_II', '']
    mask= '|'.join([p+s for s in subpanels for p in panels])

    df['Olink_panels'] = df[feat_col].str.extract(rf"({mask})$", re.IGNORECASE)
    df[output_feat_col] = [og_feat.replace("_"+panel, '') for og_feat, panel in df[[feat_col, 'Olink_panels']].to_numpy()]
    return df


def get_num_biomarkers(final_signif_df: pd.DataFrame):
    """
    Get the number of biomarkers in the dataframe. 
    This accounts for the fact that the proteomics data has multiple panels, 
    so we need to remove the panel information to get the bare features.


    :param final_signif_df: DataFrame containing the significant features.
    :return: Number of biomarkers.
    """


    tx_len = len(final_signif_df[final_signif_df['omics_type'] == 'transcriptomics'])
    px_df = final_signif_df[final_signif_df['omics_type'] == 'proteomics'].copy()

    px_df = _get_bare_features(px_df)

    px_df_no_panels = px_df[['bare_feature','drug_family', 'diagnosis']].drop_duplicates()
    
    px_len = len(px_df_no_panels)

    return tx_len + px_len
