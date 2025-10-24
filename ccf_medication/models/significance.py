"""
Significance testing helpers

This module provides functions for:
- Correcting p-values using multiple testing procedures (`p_value_correction`)
"""

import logging
import warnings
from statsmodels.stats.multitest import multipletests

from ccf_medication.constants.system import LOGGER_NAME

warnings.filterwarnings("ignore")
logger = logging.getLogger(LOGGER_NAME)


def p_value_correction(
    input_df,
    p_value_column,
    correction_method="fdr_bh",
    p_corr_column="adjusted_p_value",
):
    """
    Perform multiple testing correction on a p-value column.

    :param input_df: DataFrame containing the p-values to be corrected.
    :param p_value_column: Name of the column containing p-values.
    :param correction_method: Multiple testing correction method passed to
        ``statsmodels.stats.multitest.multipletests`` (e.g., ``"fdr_bh"``).
    :param p_corr_column: Name of the output column for adjusted p-values.
    :return: DataFrame with an added column ``p_corr_column`` containing adjusted
        p-values. Rows with NaN p-values are left unchanged.
    """
    output_df = input_df.copy()
    # Identify valid (non-NaN) p-values
    valid_mask = output_df[p_value_column].notna()

    # Perform correction only on valid p-values
    if valid_mask.sum() > 0:
        corrected_pvals = multipletests(
            output_df.loc[valid_mask, p_value_column], method=correction_method
        )[1]
        output_df.loc[valid_mask, p_corr_column] = corrected_pvals
        output_df[p_corr_column] = output_df[p_corr_column].astype(float)

    return output_df
