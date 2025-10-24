"""
Utility helpers for figure styling and plot composition

This module provides functions for:
- Building fixed color palettes for categories (`build_fixed_palette`)
- Finding significant drug-family pairs for a gene (`get_sig_drug_fam_pairs`)
"""
import seaborn as sns

from ccf_medication.constants.tables import NO_MEDICATION, FEATURE_COLUMN
from ccf_medication.constants.thresholds import ADJ_PVAL_THRESH

def build_fixed_palette(categories: list[str], palette_name: str = "muted", palette=None) -> dict[str, str]:
    """
    Build a fixed palette for the categories.

    :param categories: List of categories to build the palette for.
    :param palette_name: Optional name of the palette to use. (if you want to use a seaborn palette)
    :param palette: Optional list of colors to use for the palette. (if you want to use a custom palette)
    :return: Dictionary mapping categories to colors.
    """
    if palette is None:
        colors = sns.color_palette(palette_name, n_colors=len(categories)).as_hex()
    else:
        print(len(categories), len(palette))
        if len(categories) != len(palette):
            # extend the palette to the length of the categories by cycling through the palette
            colors = (palette * (len(categories) // len(palette))) + (palette[:len(categories) % len(palette)])
        else:
            colors = palette
    return {cat: color for cat, color in zip(categories, colors)}

# find the significant drug families
def get_sig_drug_fam_pairs(results_df, gene, subpop_condition_dict, 
                        pval_threshold=ADJ_PVAL_THRESH, fc_threshold=None, 
                        return_no_med_last=True, 
                        feature_col=FEATURE_COLUMN, 
                        pval_col=None, 
                        fc_col=None):

    """
    Find significant drug family pairs for a given gene within a subpopulation.

    :param results_df: DataFrame containing the results.
    :param gene: Feature/gene identifier to query.
    :param subpop_condition_dict: Mapping of subpopulation columns -> required values.
    :param pval_threshold: Maximum adjusted p-value to be considered significant.
    :param fc_threshold: Minimum absolute fold change to be considered significant.
    :param return_no_med_last: If True, place the no-medication group last in the returned order.
    :param feature_col: Column name for the feature identifier.
    :param pval_col: Column name for the (adjusted) p-value.
    :param fc_col: Column name for the fold change / coefficient.
    :return: Tuple of (significant_pairs, ordered_families) where significant_pairs is a
        list of (drug_family, NO_MEDICATION) tuples and ordered_families is the list of
        all families in the chosen order.
    :raises ValueError: If any of ``pval_col``, ``fc_col``, or ``fc_threshold`` is None.
    """
    for var in [pval_col, fc_col, fc_threshold]:
        if var is None:
            raise ValueError(f"pval_col, fc_col, or fc_threshold is None")

    sub_mask = (results_df[subpop_condition_dict.keys()] == subpop_condition_dict.values()).all(axis=1)
    subpop_results = results_df.loc[sub_mask]
    all_drug_families = subpop_results[(subpop_results[feature_col] == gene)].drug_family.unique().tolist()
    sig_drug_families = subpop_results[(subpop_results[feature_col] == gene) & (subpop_results[pval_col] < pval_threshold) & (subpop_results[fc_col].abs() > fc_threshold)].drug_family.unique().tolist()
    if return_no_med_last:
        all_drug_families = sorted(all_drug_families) + [NO_MEDICATION]
    else:
        all_drug_families = [NO_MEDICATION] + sorted(all_drug_families)
    return [(drug_fam, NO_MEDICATION) for drug_fam in sig_drug_families], all_drug_families
