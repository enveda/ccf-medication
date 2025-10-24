"""
File and results management utilities

This module provides functions for:
- Saving gene lists to organized directory structures (`GeneFileSaver`)
- Saving combined results to Excel with organized sheets (`save_combined_df_in_excel`)
- Sanitizing strings for safe filenames (`sanitize_filename`)
- Generating filenames (`generate_filename`)
"""

import os
from pathlib import Path
from typing import List, Optional, Tuple
import re

import pandas as pd



class GeneFileSaver:
    """
    Handles saving genes to organized directory structures.
    """

    def __init__(self, root_dir: str):
        """
        :param root_dir: str - The root directory to save the files
        """
        self.root_dir = Path(root_dir)

    def save_genes(
        self,
        df: pd.DataFrame,
        dir_level_cols: List[str],
        file_level_col: str,
        col_to_save: str,
        file_suffix: Optional[str] = None,
    ) -> None:
        """Save genes from a drug family to organized directory structure.

        :param df: pd.DataFrame - The dataframe to save
        :param dir_level_cols: List[str] - The columns to use to create the directory structure
        :param file_level_col: str - The column to use to create the file name
        :param col_to_save: str - The column values to save
        :param file_suffix: Optional[str] - The suffix to add to the file name
        """
        group_cols = dir_level_cols + [file_level_col]

        for group_keys, group_data in df.groupby(group_cols, group_keys=True):
            self._save_single_group(
                group_data,
                group_keys,
                dir_level_cols,
                file_level_col,
                col_to_save,
                file_suffix,
            )

    def _save_single_group(
        self,
        group_data: pd.DataFrame,
        group_keys: Tuple,
        dir_level_cols: List[str],
        file_level_col: str,
        col_to_save: str,
        file_suffix: Optional[str],
    ) -> None:
        """
        Save a single group's data to file.

        :param group_data: pd.DataFrame - The dataframe to save
        :param group_keys: Tuple - The keys to use to create the directory structure
        :param dir_level_cols: List[str] - The columns to use to create the directory structure
        :param file_level_col: str - The column to use to create the file name
        :param col_to_save: str - The column values to save
        :param file_suffix: Optional[str] - The suffix to add to the file name
        """
        # Handle both single and multiple grouping columns
        if len(dir_level_cols) + 1 == 1:  # Single column
            group_keys = (group_keys,)

        # Create directory path
        dir_values = [str(key) for key in group_keys[:-1]]  # All except file level
        dir_path = self.root_dir / Path(*dir_values)
        dir_path.mkdir(parents=True, exist_ok=True)

        # Create file name
        base_name = sanitize_filename(str(group_keys[-1]))  # File level value
        if file_suffix:
            file_name = f"{base_name}_{file_suffix}.txt"
        else:
            file_name = f"{base_name}.txt"

        file_path = dir_path / file_name

        # Save values to file
        values = group_data[col_to_save].tolist()
        with open(file_path, "w") as f:
            f.write("\n".join(map(str, values)))


def save_combined_df_in_excel(df: pd.DataFrame, output_path: str):
    """
    Save combined dfs in excel.
    :param df: DataFrame to save.
    :param output_path: Path to save the excel file.
    """
    with pd.ExcelWriter(output_path) as writer:
        for _, row in (
            df[["omics_type", "diagnosis", "simple_tissue"]]
            .drop_duplicates()
            .iterrows()
        ):
            omics_type = row["omics_type"]
            diagnosis = row["diagnosis"]
            simple_tissue = row["simple_tissue"]
            df_subset = df[
                (df["omics_type"] == omics_type)
                & (df["diagnosis"] == diagnosis)
                & (df["simple_tissue"] == simple_tissue)
            ]
            if omics_type == "proteomics":
                sheet_name = f"{omics_type}_{diagnosis}"
                df_subset = df[
                    (df["omics_type"] == omics_type) & (df["diagnosis"] == diagnosis)
                ].copy()
                df_subset.drop(columns=["diagnosis", "simple_tissue", "omics_type"], inplace=True)
            else:
                sheet_name = f"{omics_type}_{diagnosis}_{simple_tissue}"
                df_subset = df[
                    (df["omics_type"] == omics_type)
                    & (df["diagnosis"] == diagnosis)
                    & (df["simple_tissue"] == simple_tissue)
                ].copy()
                df_subset.drop(columns=["diagnosis", "simple_tissue", "omics_type"], inplace=True)
            
            # put the drug family and feature column at the beginning of the dataframe, followed by "coef" columns and then adjusted pval columns
            first_cols = ["drug_family", "feature"]
            coef_substring = ["coef"]
            adj_pval_substrings = ["adj_p_val", "adjusted_p_value"]
            coef_pattern = re.compile("|".join(map(re.escape, coef_substring)), flags=re.I)
            adj_pval_pattern = re.compile("|".join(map(re.escape, adj_pval_substrings)), flags=re.I)

            # Column names that contain any of the substrings (case-insensitive)
            coef_cols = df.columns[df.columns.str.contains(coef_pattern)].tolist()
            adj_pval_cols = df.columns[df.columns.str.contains(adj_pval_pattern)].tolist()
            cols = first_cols + coef_cols + adj_pval_cols

            end_cols = [col for col in df_subset.columns if col not in cols]
            df_subset = df_subset[cols + end_cols]
            df_subset.to_excel(writer, sheet_name=sheet_name, index=False)


def sanitize_filename(name: str) -> str:
    """
    Replace characters that are problematic in filenames.

    :param name: str, the original filename or string
    :return: str, the sanitized string safe for use as filename
    """
    replacements = {
        "/": "_",
        "\\": "_",
        ":": "_",
        "*": "_",
        "?": "_",
        '"': "_",
        "<": "_",
        ">": "_",
        "|": "_",
        " ": "_",
    }

    sanitized = name
    for old, new in replacements.items():
        sanitized = sanitized.replace(old, new)

    return sanitized


def generate_filename(omics_type: str, split_values: List[str]) -> str:
    """
    Generate standardized filename for results.
    :param omics_type: Omics type string used as suffix (e.g., "proteomics").
    :param split_values: List of split values to include in the filename.
    :return: A sanitized filename including any split values, and the standard suffix for the provided omics type.
    """
    filename = "_".join([ sanitize_filename(split_value) for split_value in split_values])
    suffix = f"{omics_type}_fold_changes.csv"
    return filename + "_" + suffix
