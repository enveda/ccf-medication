"""
Miscellaneous helpers used across analysis and plotting

This module provides functions for:
- Combining proteomics and transcriptomics tables for heatmap-style outputs (`combine_px_and_tx_tables`)
- Adding a disease severity column to a table (`add_disease_severity_column`)
- Cleaning column names by normalizing separators (`clean_columns`)
"""

from typing import Dict
import pandas as pd


def combine_px_and_tx_tables(
    px_table: pd.DataFrame, tx_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Combine the transcriptomics and proteomics tables.
    This function is particularly useful for combining the outputs of the count_heatmap function.

    :param px_table: The proteomics table.
    :param tx_table: The transcriptomics table.
    :return: The combined table.
    """
    px_df = px_table.copy()
    tx_df = tx_table.copy()
    # adjust the proteomics table so it combines nicely with the transcriptomics table
    px_existing_names = px_table.columns.names  # keep original column names
    px_df.columns = pd.MultiIndex.from_product([["Proteomics"], px_df.columns, ["all"]])
    px_df.columns.names = ["omics_type"] + px_existing_names + ["simple_tissue"]

    # adjust the transcriptomics table so it combines nicely with the proteomics table
    # add omics_type as top level, keep existing column levels (incl. simple_tissue)
    tx_existing_names = tx_table.columns.names
    tx_df = pd.concat({"Transcriptomics": tx_df}, axis=1)
    tx_df.columns.names = ["omics_type"] + tx_existing_names

    # combine the tables
    return pd.concat([tx_df, px_df], axis=1)


def add_disease_severity_column(
    table: pd.DataFrame,
    disease_severity_col: str = "disease_severity",
    source_col: str = "endo_category",
    disease_severity_map: Dict[str, str] = None,
) -> pd.DataFrame:
    """
    Add a disease severity column to the table.

    Default mapping: remission -> 0; {moderate, severe} -> 1.
    :param table: The table to add the disease severity column to.
    :param disease_severity_col: The name of the disease severity column.
    :param source_col: The name of the source column.
    :param disease_severity_map: The mapping from source values to disease severity values.
    :return: The table with the disease severity column added.
    """

    if disease_severity_map is None:
        disease_severity_map = {"remission": 0, "severe": 1, "moderate": 1}
    table = table.copy()
    table[disease_severity_col] = table[source_col].map(disease_severity_map)
    return table

def clean_columns(input_list):
    """
    Normalize column names by replacing problematic characters.

    :param input_list: List of column names to clean.
    :return: List of cleaned column names.
    """
    return [x.replace("-", "_").replace("/", "_") for x in list(input_list)]