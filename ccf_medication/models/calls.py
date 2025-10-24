"""
Mixed effects model analysis orchestrators.

This module provides functions for:
- Running mixed effects model analysis (`run_me_model_analysis`)

helper functions:
- Generating split combinations (`_generate_split_combinations`)
- Checking reference group size (`_check_reference_group_size`)
- Finding valid non-reference groups (`_find_valid_non_reference_groups`)

"""

import os
from typing import Any, Dict, List, Optional, Tuple
from functools import reduce

import pandas as pd

# Specific imports from ccf_medication modules
from ccf_medication.models.mixed_effects import mixed_effects_model
from ccf_medication.utils.file_management import generate_filename

# Constants
from ccf_medication.constants.tables import PATIENT_ID_COL


def _generate_split_combinations(
    data: pd.DataFrame, cols_to_split_on: Dict[str, List[str]], exclude_values: Dict[str, List[str]]
) -> List[Dict[str, Any]]:
    """Generate all combinations of split column values.

    :param data: Input dataframe containing the split columns.
    :param cols_to_split_on: Mapping from column name to a list of values to split the data on (if no list is provided, splits on all unique values of the column).
    :param exclude_values: Mapping from column name to a list of values to exclude
        from the combinations.
    :return: List of dictionaries, one per unique combination, mapping split column
        names to selected values. Returns ``[{}]`` when ``split_columns`` is
        empty.
    """
    if not cols_to_split_on:
        return [{}]

    data_subset = data.copy()
    
    # drop rows with values in exclude_values
    for col, values in exclude_values.items():
        data_subset = data_subset[~data_subset[col].isin(values)]

    # drop rows without values in cols_to_split_on if a list of values is provided 
    # no values are dropped for the given column
    for col, values in cols_to_split_on.items():
        if values is not None and len(values) > 0: 
            data_subset = data_subset[data_subset[col].isin(values)]

    # drop na values
    data_subset.dropna(axis=0, how='any', inplace=True)

    unique_splits = data_subset[list(cols_to_split_on.keys())].drop_duplicates().to_dict(orient='records')

    return unique_splits


def _check_reference_group_size(
    subset_data: pd.DataFrame,
    reference_groups: Dict[str, Any],
    min_patients: int,
    split_combo: Dict[str, Any],
) -> Tuple[bool, Dict[str, Any]]:
    """
    Check if the reference (No Medication, or remission) group is large enough.
    
    :param subset_data: Subset of the data to check the reference group size for.
    :param reference_groups: Reference groups to check the size for.
    :param min_patients: Minimum number of patients required for the reference group.
    :param split_combo: Split combination to check the reference group size for.
    :return: Tuple of True if the reference group is large enough, False otherwise, and a dictionary of the group counts.
    """
    ref_large_enough = True
    group_counts = {}
    for col, value in reference_groups.items():
        num_patients = subset_data[subset_data[col] == value][PATIENT_ID_COL].nunique()
        if not isinstance(num_patients, int):
            raise ValueError(f"when checking the reference group {col} = {value}, the number of patients is not an int")
        if num_patients < min_patients:
            print(f"Not enough patients in the {'_'.join(split_combo.values())} for the reference group: {col} = {value}")
            print(f"Number of patients: {num_patients}")
            ref_large_enough = False
            break
        else:
            num_samples = subset_data[subset_data[col] == value].shape[0]
            group_counts[value] = {
                f"num_patients_{value}": num_patients,
                f"num_samples_{value}": num_samples,
            }

    return ref_large_enough, group_counts

def _find_valid_non_reference_groups(
    subset_data: pd.DataFrame,
    primary_effect_col: str,
    min_patients: int,
    reference_groups: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Find the non-reference groups that are large enough.

    :param subset_data: Subset of the data to check the non-reference group size for.
    :param primary_effect_col: Column to check the non-reference group size for.
    :param min_patients: Minimum number of patients required for the non-reference group.
    :param reference_groups: Reference groups to check the non-reference group size for.
    :return: Dictionary of the non-reference group counts.
    """
    non_ref_values = subset_data[subset_data[primary_effect_col] != reference_groups[primary_effect_col]][primary_effect_col].unique()
    valid_non_ref_values = {}
    for value in non_ref_values:
        num_patients = subset_data[subset_data[primary_effect_col] == value][PATIENT_ID_COL].nunique()
        if not isinstance(num_patients, int):
            raise ValueError(f"when checking the non-reference group {value}, the number of patients is not an int")
        if num_patients >= min_patients:
            num_samples = subset_data[subset_data[primary_effect_col] == value].shape[0]
            valid_non_ref_values[value] = {
                f"num_patients_non_reference": num_patients,
                f"num_samples_non_reference": num_samples,
            }
    return valid_non_ref_values

def run_me_model_analysis(
    data: pd.DataFrame,
    gene_cols: List[str],
    min_patients: int = None,
    omics_type: str = None,
    cols_to_split_on: Dict[str, List[str]] = None,
    primary_effect_col: str = None,
    reference_groups: Dict[str, Any] = None,
    additional_fixed_effects: Optional[List[str]] = [],
    random_effects: Optional[List[str]] = [PATIENT_ID_COL],
    exclude_values: Optional[Dict[str, List[str]]] = None,
    save_path: str = None,
    save_filename: str = None,
) -> None:
    """
    Run a mixed effects model analysis for all subpopulations based on the columns in `cols_to_split_on`.

    :param data: Input dataframe containing the source column.
    :param gene_cols: List of gene columns to include in the analysis.
    :param min_patients: Minimum number of patients required for analysis (e.g. minimum number of patients per drug family or minimum number of patients per disease severity group).
    :param omics_type: Type of omics data (e.g."transcriptomics", "proteomics").
    :param cols_to_split_on: Columns to split the data on (e.g. {"diagnosis": ["uc", "cd"], "simple_tissue": ["colon", "small_intestine"]}).
    :param primary_effect_col: Column to analyze (e.g. "disease_severity", or "drug_family").
    :param reference_groups: Reference groups to use for the analysis (e.g. {"disease_severity": "remission"} or {"drug_family": NO_MEDICATION}).
    :param additional_fixed_effects: Additional fixed effects to include in the analysis (e.g. ["age", "gender"]).
    :param random_effects: Random effects to include in the analysis (e.g. ["patient_id"]).
    :param exclude_values: Values to exclude from the analysis (e.g. {"drug_family": ["Combination Therapy"]}).
    :param save_path: Path to save the results dataframe.
    :return: None.
    """
    # Check all inputs are valid
    if min_patients is None:
        raise ValueError("min_patients must be provided")
    if exclude_values is None:
        exclude_values = {}
    if omics_type is None:
        raise ValueError("omics_type must be provided, e.g. 'transcriptomics' or 'proteomics'")
    if primary_effect_col is None:
        raise ValueError("primary_effect_col must be provided, e.g. 'disease_severity' or 'drug_family'")
    if reference_groups is None:
        raise ValueError("reference_groups must be provided, e.g. {'disease_severity': 'remission'} or {'drug_family': NO_MEDICATION}")
    if save_path is None:
        raise ValueError("save_path must be provided")
    else:
        os.makedirs(save_path, exist_ok=True)


    split_combos = _generate_split_combinations(
        data, cols_to_split_on, exclude_values
    )

    agg_results = []

    # loop through each split combo and run the mixed effects model
    # (e.g. for disease severity transcriptomics analysis, we loop 
    # through each diagnosis, simple tissue, and drug_family combination)
    for split_combo in split_combos:
        subset_data = data.copy()
        conditions = [subset_data[col] == value for col, value in split_combo.items()]
        mask = reduce(lambda x, y: x & y, conditions)
        subset_data = subset_data[mask].copy()

        if subset_data.empty:
            print(f"Subset data is empty for split combo: {split_combo}")
            continue   

        # Checking to see if the reference group is large enough
        ref_large_enough, group_counts = _check_reference_group_size(subset_data, reference_groups, min_patients, split_combo)
        if not ref_large_enough:
            continue

        # finding the non-reference values that are large enough 
        # otherwise it is left out of the analysis
        valid_non_ref_values = _find_valid_non_reference_groups(subset_data, primary_effect_col, min_patients, reference_groups)
            
        valid_values_for_primary_effect_col = list(valid_non_ref_values.keys()) + [reference_groups[primary_effect_col]]
        filtered_subset_data = subset_data[subset_data[primary_effect_col].isin(valid_values_for_primary_effect_col)]

        if filtered_subset_data.empty:
            print(f"Filtered subset data is empty for split combo: {split_combo}")
            continue

        fixed_effects = [primary_effect_col] + additional_fixed_effects

        try:
            model_results_dict = mixed_effects_model(
                original_input_data=filtered_subset_data.reset_index(drop=True),
                input_features_list=gene_cols,
                missing_data_handling="drop",
                fixed_effects=fixed_effects,
                random_effects=random_effects,
                reference_groups_dictionary=reference_groups,
                correction="fdr_bh",
                zero_intercept=False,
                min_number_of_samples=min_patients,
            )

            results_df = model_results_dict["output_scores"]

            # Drop columns without "reference" in the effect column
            results_df = results_df[results_df["effect"].str.contains("reference")]

            # Adding the split combo to the results dataframe columns
            for col, value in split_combo.items():
                results_df[col] = value

            # add column for non-reference group values
            results_df["non_reference_values"] = results_df["effect"].apply(lambda x: x.split("[T.")[1].split("]")[0])
            results_df["reference_values"] = results_df["effect"].apply(lambda x: x.split("(reference='")[1].split("')")[0])

            # these add columns for the number of patients and samples in the reference group (for all results rows)
            for col, count_dict in group_counts.items():
                for key, value in count_dict.items():
                    results_df[f'{key}'] = int(value)

            # these add columns for the number of patients and samples in the non-reference group (only for the rows that correspond to the non-reference values)
            for non_ref_value, count_dict in valid_non_ref_values.items():
                for col_name, count_value in count_dict.items():
                    results_df.loc[results_df["non_reference_values"] == non_ref_value, f'{col_name}'] = int(count_value)

            # Save the results dataframe
            filename = generate_filename(omics_type, split_combo.values())
            results_df.to_csv(os.path.join(save_path, filename), index=False)

            agg_results.append(results_df)
            

        except Exception as e:
            print(f"Mixed effects model failed: {e}")
            return None

    # save the aggregated results dataframe as a parquet file
    aggregated_results_df = pd.concat(agg_results, ignore_index=True)

    if save_filename is not None:
        aggregated_results_df.to_parquet(os.path.join(save_path, f"{save_filename}"))

    return aggregated_results_df