"""
Model fitting utilities for mixed effects and OLS

This module provides functions for:
- Building patsy formula strings for fixed effects (`formula_builder`)
- Fitting mixed-effects/OLS models for single and multiple features (`individual_mixed_effects_model`, `mixed_effects_model`)
"""

import logging
import warnings
from multiprocessing import Pool

import pandas as pd
import statsmodels.formula.api as smf
from tqdm import tqdm

from ccf_medication.constants.system import LOGGER_NAME
from ccf_medication.models.significance import p_value_correction
from ccf_medication.utils.miscellaneous import clean_columns

warnings.filterwarnings("ignore")
logger = logging.getLogger(LOGGER_NAME)


def formula_builder(
    fixed_effects,
    continuous_fixed_effects,
    reference_groups_dictionary,
    zero_intercept,
):
    """
    Build a fixed-effects formula segment for mixed effects/OLS models.

    :param fixed_effects: List of fixed effects to include in the model.
    :param continuous_fixed_effects: List of continuous fixed effects (subset of
        ``fixed_effects``). Used for validation only.
    :param reference_groups_dictionary: Optional mapping from effect to reference
        level to apply treatment coding.
    :param zero_intercept: If True, omit the intercept (use ``"0 + ..."``).
    :return: String representing the right-hand side of a patsy formula for the fixed effects.
    """
    # Check if continuous fixed effects are in fixed_effects
    if continuous_fixed_effects is not None:
        for effect in continuous_fixed_effects:
            if effect not in fixed_effects:
                logger.warning(f"{effect} is not in fixed_effects, skipping")

    fixed_effects_terms = []
    for effect in fixed_effects:
        if (
            reference_groups_dictionary is not None
            and effect in reference_groups_dictionary
        ):
            logger.info(
                f"Reference group {reference_groups_dictionary[effect]} found for {effect}"
            )
            term = f"C({effect}, Treatment(reference='{reference_groups_dictionary[effect]}'))"
        else:
            term = effect
        fixed_effects_terms.append(term)
    if zero_intercept:
        fixed_effects_string = "0 + " + " + ".join(fixed_effects_terms)
    else:
        fixed_effects_string = " + ".join(fixed_effects_terms)
    return fixed_effects_string


def individual_mixed_effects_model(
    input_data,
    input_feature,
    random_effects,
    input_formula,
    min_number_of_samples,
    missing_data_handling,
    max_iterations=100,
):
    """
    Fit a mixed effects (or OLS if no random effects) model to one feature.

    :param input_data: DataFrame containing the data.
    :param input_feature: Name of the feature/response column to model.
    :param random_effects: Optional list of random-effects grouping columns. If not
        provided, an OLS model is fit instead.
    :param input_formula: Patsy fixed-effects formula string (right-hand side).
    :param min_number_of_samples: Minimum number of samples required to fit the
        model.
    :param missing_data_handling: Strategy for handling missing data ("drop",
        "none", or "raise").
    :param max_iterations: Maximum number of iterations for model fitting.
    :return: Tuple of ``(results, predictions)`` where ``results`` is a list of
        dictionaries per effect and ``predictions`` is a ``pd.Series`` of model
        predictions for the provided feature.
    """
    current_formula = f"{input_feature} ~ {input_formula}"

    # Raise error if NAs are present in random effects
    if random_effects:
        if input_data[random_effects].isna().any().any():
            raise ValueError("NAs present in random effects, please remove or impute")

    # Remove NAs
    feature_data = input_data.dropna(subset=input_feature)

    if len(feature_data) < min_number_of_samples:  # Minimum sample size check
        logger.warning(f"Too few samples for {input_feature}, skipping")
        return None, None

    # Fit model
    if random_effects:
        primary_group = random_effects[0]
        if len(random_effects) > 1:
            additional = " + ".join(f"0 + {re}" for re in random_effects[1:])
            re_formula = f"1 + {additional}" if additional else "1"
        else:
            re_formula = "1"

        model = smf.mixedlm(
            current_formula,
            feature_data,
            groups=feature_data[primary_group],
            re_formula=re_formula,
            missing=missing_data_handling,
        )
    else:
        model = smf.ols(current_formula, feature_data, missing=missing_data_handling)

    try:
        fitted_model = model.fit(maxiter=max_iterations)

    except Exception as e:
        logger.error(f"Error fitting model for {input_feature}: {str(e)}")
        return None, None

    params = fitted_model.params
    ses = fitted_model.bse
    pvals = fitted_model.pvalues

    results = [
        {
            "feature": input_feature,
            "effect": eff,
            "coefficient": params[eff],
            "standard_error": ses[eff],
            "p_value": pvals[eff],
        }
        for eff in params.index
    ]

    return results, fitted_model.predict(feature_data)


def _individual_model(
    input_data,
    input_feature,
    random_effects,
    formula_string,
    min_number_of_samples,
    missing_data_handling,
    max_iterations=100,
):
    """
    Wrapper around ``individual_mixed_effects_model`` adding the feature name.

    :param input_data: DataFrame containing the data.
    :param input_feature: Feature/response column to fit.
    :param random_effects: List of random-effects grouping columns.
    :param formula_string: Patsy fixed-effects formula string.
    :param min_number_of_samples: Minimum samples required to fit.
    :param missing_data_handling: Strategy for handling missing data ("drop",
        "none", or "raise").
    :param max_iterations: Maximum iterations for model fitting.
    :return: Tuple[List[dict], pd.Series, str]: Results list, predictions, and the
    """
    try:
        current_results, current_predicted_values = individual_mixed_effects_model(
            input_data,
            input_feature,
            random_effects,
            formula_string,
            min_number_of_samples,
            missing_data_handling,
            max_iterations,
        )
        if current_results is not None and current_predicted_values is not None:
            return current_results, current_predicted_values, input_feature
    except Exception as error:
        logger.error(f"Error processing {input_feature}: {str(error)}")
        return None, None, None





def mixed_effects_model(
    original_input_data,
    input_features_list,
    missing_data_handling,
    random_effects=None,
    fixed_effects=None,
    continuous_fixed_effects=None,
    reference_groups_dictionary=None,
    correction="fdr_bh",
    zero_intercept=True,
    min_number_of_samples=10,
    custom_formula=None,
    max_iterations=100,
):
    """
    Fit mixed-effects/OLS models across multiple features with p-value correction.

    :param original_input_data: DataFrame containing covariates, grouping variables,
        and feature columns.
    :param input_features_list: List of feature/response column names to model.
    :param missing_data_handling: Strategy for handling missing data ("drop",
        "none", or "raise").
    :param random_effects: Optional list of random-effects grouping columns. If not
        provided, OLS is used.
    :param fixed_effects: Optional list of fixed-effects columns used to build the
        formula when ``custom_formula`` is not provided.
    :param continuous_fixed_effects: Optional list of continuous fixed-effects
        columns (subset of ``fixed_effects``).
    :param reference_groups_dictionary: Optional mapping of effect -> reference level
        for treatment coding. When provided, the intercept is included
        regardless of ``zero_intercept``.
    :param correction: Multiple testing correction method (e.g., ``"fdr_bh"``).
    :param zero_intercept: If True, use a zero-intercept formula unless reference
        groups are specified.
    :param min_number_of_samples: Minimum number of samples required to fit each
        feature.
    :param custom_formula: Optional fixed-effects formula string to use instead of
        building from ``fixed_effects``.
    :param max_iterations: Maximum number of iterations for model fitting.
    :return: Dictionary with the following keys:
        - ``output_scores``: DataFrame of per-feature/effect coefficients,
          standard errors, p-values, and adjusted p-values.
        - ``predicted_values``: DataFrame of predictions per feature along
          with identifier/non-feature columns.
    :raises ValueError: For invalid arguments or missing required columns.
    """
    if reference_groups_dictionary is not None:
        zero_intercept = False
    # Validate that fixed_effects is not empty when reference_groups_dictionary is provided
    if reference_groups_dictionary is not None and (
        fixed_effects is None or len(fixed_effects) == 0
    ):
        raise ValueError(
            "reference_groups_dictionary provided but fixed_effects is empty"
        )
    if missing_data_handling not in ["drop", "none", "raise"]:
        raise ValueError(
            f"Invalid missing data handling method: {missing_data_handling}, must be one of 'drop', 'none' or 'raise'"
        )

    input_data = original_input_data.copy()
    input_data.columns = clean_columns(input_data.columns)
    input_features_list = clean_columns(input_features_list)
    non_feature_cols = [
        col for col in input_data.columns if col not in input_features_list
    ]
    ids = input_data[non_feature_cols].astype(str).copy()
    if continuous_fixed_effects:
        input_data[continuous_fixed_effects] = input_data[
            continuous_fixed_effects
        ].apply(pd.to_numeric, errors="coerce")
    input_data[input_features_list] = input_data[input_features_list].apply(
        pd.to_numeric, errors="coerce"
    )
    if random_effects is not None:
        # Check if random effects columns exist in the data
        for effect in random_effects:
            if effect not in input_data.columns:
                raise ValueError(f"Random effect '{effect}' not found in input data")

        # Check if we have enough groups for each random effect
        for effect in random_effects:
            group_counts = input_data[effect].value_counts()
            if len(group_counts) < 2:
                logger.warning(
                    f"Need at least 2 groups for random effect '{effect}', this effect might cause issues"
                )

    if custom_formula is not None:
        formula_string = custom_formula
        logger.info(f"Using custom formula: {formula_string}")

    elif fixed_effects is not None:
        # Create Treatment coding with reference for each fixed effect
        formula_string = formula_builder(
            fixed_effects,
            continuous_fixed_effects,
            reference_groups_dictionary,
            zero_intercept,
        )
    else:
        logger.warning("No fixed effects provided, stopping")
        return None

    covariate_cols = (random_effects or []) + (fixed_effects or [])
    input_data = input_data.dropna(subset=covariate_cols)

    logger.info(
        f"""Running mixed effects model with {len(input_features_list)} features
        - Formula: {formula_string}
        - Missing data handling: {missing_data_handling}
        - Zero intercept: {zero_intercept}
        - Reference groups dictionary: {reference_groups_dictionary}
        - Min number of samples: {min_number_of_samples}
        """
    )

    worker_args = [
        (
            input_data,
            feature,
            random_effects,
            formula_string,
            min_number_of_samples,
            missing_data_handling,
            max_iterations,
        )
        for feature in input_features_list
    ]

    results, predicted_values, processed = [], [], []
    with Pool() as pool:
        for out in tqdm(
            pool.starmap(_individual_model, worker_args), total=len(worker_args)
        ):
            try:
                res, preds, feat = out
                results.extend(res)
                predicted_values.append(preds)
                processed.append(feat)
            except Exception as e:
                logger.warning(f"Error processing {feat}: {str(e)}")
                continue

    if not results:
        logger.error("No models converged – no results to report.")
        return None

    output_df = pd.DataFrame(results)

    # Ensure the 'effect' column exists before groupby
    if "effect" not in output_df.columns:
        logger.error("'effect' column not found in results DataFrame")
        return {"output_scores": output_df, "predicted_values": pd.DataFrame()}

    # Group by effect to perform multiple testing correction
    logger.info(f"Performing multiple testing correction with {correction} method")
    output_df = output_df.groupby("effect").apply(
        lambda x: p_value_correction(x, "p_value", correction_method=correction)
    )
    output_df = output_df.reset_index(drop=True)

    if len(predicted_values) > 0:
        predicted_values = pd.concat(predicted_values, axis=1)
        predicted_values.columns = processed
        # non feature columns
        non_feature_cols = [
            col for col in input_data.columns if col not in input_features_list
        ]
        predicted_values = pd.concat(
            [ids.reset_index(drop=True), predicted_values.reset_index(drop=True)],
            axis=1,
        )
    else:
        predicted_values = pd.DataFrame()

    output_dictionary = {
        "output_scores": output_df,
        "predicted_values": predicted_values,
    }

    return output_dictionary
