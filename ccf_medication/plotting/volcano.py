"""
Functions used to create volcano plots for differential analysis

This module provides functions for:
- Creating Matplotlib volcano plots with proximity-based labeling (`plot_volcano_plot`)
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from typing import List, Optional
from pathlib import Path

from ccf_medication.constants.tables import (
    FOLD_CHANGE_COLUMN,
    P_VALUE_COLUMN,
    DRUG_FAMILY_COL,
    FEATURE_COLUMN,
)
from ccf_medication.constants.thresholds import ADJ_PVAL_THRESH

from ccf_medication.utils.file_management import sanitize_filename


def plot_volcano_plot(
    input_df,
    p_value_col=P_VALUE_COLUMN,
    fold_change_col=FOLD_CHANGE_COLUMN,
    p_value_threshold=0.05,
    p_value_floor=1e-10,
    fold_change_threshold=0.5,
    p_value_drop_threshold=None,
    fold_change_drop_threshold=None,
    save_path=None,
    ignore_significance=False,
    custom_annotation_column=None,
    ax=None,
    title=None,
    label_size=10,
    max_labels=50,
    proximity_threshold=10,
    fc_weight=10,
    label_offset_x=5,
    label_offset_y=5,
    x_axis_label="Fold Change",
    y_axis_label="-log10(p-value)",
    dpi=300,
):
    """Create a Matplotlib volcano plot with proximity-based labeling.

    :param input_df: DataFrame containing fold change and (adjusted) p-values.
    :param p_value_col: Name of the p-value column.
    :param fold_change_col: Name of the fold change column.
    :param p_value_threshold: Significance threshold for p-values.
    :param p_value_floor: Minimum p-value used to avoid log(0).
    :param fold_change_threshold: Fold change magnitude threshold for drawing
        significant points.
    :param p_value_drop_threshold: If provided, rows with p-values above this are
        dropped before plotting.
    :param fold_change_drop_threshold: If provided, rows with abs(fold change)
        below this are dropped before plotting.
    :param save_path: When provided, the figure is saved to this path and closed.
    :param ignore_significance: If True, skip significance-based filtering.
    :param custom_annotation_column: Optional column to use for point labels.
    :param ax: Optional Matplotlib axes to draw on.
    :param title: Plot title.
    :param label_size: Marker size for points.
    :param max_labels: Maximum number of text labels added using proximity rules.
    :param proximity_threshold: Proximity radius (in data units) for de-duplicating
        text labels.
    :param fc_weight: Weight applied to fold change for label priority.
    :param label_offset_x: X offset for annotation text (points).
    :param label_offset_y: Y offset for annotation text (points).
    :param x_axis_label: Label for the x-axis.
    :param y_axis_label: Label for the y-axis.
    :param dpi: DPI used when saving the figure.
    """
    input_df = input_df.copy()

    # Replace p values of 0 with the floor value
    input_df[p_value_col] = input_df[p_value_col].replace(0, p_value_floor)

    # Apply filters
    if p_value_drop_threshold is not None:
        input_df = input_df[input_df[p_value_col] < p_value_drop_threshold]
    if fold_change_drop_threshold is not None:
        input_df = input_df[abs(input_df[fold_change_col]) > fold_change_drop_threshold]

    # Handle significance
    if ignore_significance:
        significant_features = input_df[p_value_col] >= 0
    else:
        significant_features = input_df[p_value_col] < p_value_threshold

    # Log transform the p-value for volcano plot
    input_df[p_value_col] = -np.log10(input_df[p_value_col])

    # Calculate combined score: -log10(p-value) + fc_weight * |log2(fold change)|
    input_df["combined_score"] = input_df[p_value_col] + fc_weight * np.abs(
        input_df[fold_change_col]
    )

    # Create a new axis if none is provided
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # Plot all points (non-significant and significant)
    ax.scatter(
        input_df[fold_change_col],
        input_df[p_value_col],
        color="gray",
        s=2,
        alpha=0.5,
        label="Non-Significant",
    )

    # Separate the significant features above and below the fold change threshold
    above_fc_threshold = input_df[
        (input_df[fold_change_col] > fold_change_threshold) & significant_features
    ]

    below_fc_threshold = input_df[
        (input_df[fold_change_col] < -fold_change_threshold) & significant_features
    ]

    ax.scatter(
        above_fc_threshold[fold_change_col],
        above_fc_threshold[p_value_col],
        s=2,
        color="red",
        label=f"Significant - FC > {fold_change_threshold}",
    )

    ax.scatter(
        below_fc_threshold[fold_change_col],
        below_fc_threshold[p_value_col],
        s=2,
        color="blue",
        label=f"Significant - FC < {-fold_change_threshold}",
    )

    # Concatenate above and below fold change into a single DataFrame
    thresholded_df = pd.concat([above_fc_threshold, below_fc_threshold])

    # Step 1: Sort points by their combined score (p-value and fold change)
    thresholded_df = thresholded_df.sort_values(by=["combined_score"], ascending=False)

    # Step 2: Build a cKDTree for fast proximity search
    kd_tree = cKDTree(thresholded_df[[fold_change_col, p_value_col]].values)

    # Step 3: Label points based on proximity filtering
    labeled_points = []
    labeled_indices = set()

    for index, row in thresholded_df.iterrows():
        if len(labeled_points) >= max_labels:
            break

        # Find points within a proximity radius
        neighbors = kd_tree.query_ball_point(
            [row[fold_change_col], row[p_value_col]], proximity_threshold
        )

        # If no labeled neighbors within the threshold, label this point
        if not any(n in labeled_indices for n in neighbors):
            labeled_points.append(row)
            labeled_indices.add(index)

    # Step 4: Add annotations to the labeled points with offsets
    for _, row in pd.DataFrame(labeled_points).iterrows():
        ax.annotate(
            row[custom_annotation_column] if custom_annotation_column else row.name,
            xy=(row[fold_change_col], row[p_value_col]),
            textcoords="offset points",
            xytext=(label_offset_x, label_offset_y),
            ha="right",
            fontsize=label_size,
        )

    # Draw thresholds
    ax.axvline(x=fold_change_threshold, color="black", linestyle="--", linewidth=1)
    ax.axvline(x=-fold_change_threshold, color="black", linestyle="--", linewidth=1)
    ax.axhline(
        y=-np.log10(p_value_threshold), color="black", linestyle="--", linewidth=1
    )

    # Customize plot
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel(y_axis_label)
    if title is not None:
        ax.set_title(title)
    ax.legend(
        loc="center",
        bbox_to_anchor=(0.5, -0.2),
        ncol=3,
        fontsize=8,
    )

    # Save or show the plot
    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    elif ax is None:
        return fig
    else:
        return None

def plot_volcano_plots_with_labels(
    results_df: pd.DataFrame,
    subpop_cols: List[str],
    *,
    drug_family_col: str = DRUG_FAMILY_COL,
    feature_col: str = FEATURE_COLUMN,
    p_value_col: str = P_VALUE_COLUMN,
    fold_change_col: str = FOLD_CHANGE_COLUMN,
    p_value_threshold: float = ADJ_PVAL_THRESH,
    p_value_floor: float = 1e-10,
    fold_change_threshold: float = None,
    p_value_drop_threshold: Optional[float] = None,
    fold_change_drop_threshold: Optional[float] = None,
    output_path: Optional[str] = None,
    sub_dir_levels: Optional[List[str]] = None,
    file_suffix: str = "volcano",
    interactive: bool = False,
    label_size: int = 10,
    max_labels: int = 50,
    proximity_threshold: int = 10,
    fc_weight: int = 10,
    title_prefix: Optional[str] = None,
) -> None:
    """
    Create volcano plots for every subpopulation combination and drug family.

    - Filters `results_df` to each unique set of `subpop_cols`
    - Within each subpopulation, iterates over `drug_family_col`
    - Saves plots nested under `output_path/<sub_dir_levels...>/<subpop values...>/`
    - Uses Plotly when `interactive=True` (HTML), otherwise Matplotlib (PNG)
    :param results_df: pd.DataFrame - The results dataframe
    :param subpop_cols: List[str] - The columns to group by
    :param drug_family_col: str - The column name for the drug family
    :param feature_col: str - The column name for the feature
    :param p_value_col: str - The column name for the p-value
    :param fold_change_col: str - The column name for the fold change
    :param p_value_threshold: float - The p-value threshold
    :param p_value_floor: float - The p-value floor
    :param fold_change_threshold: float - The fold change threshold
    :param p_value_drop_threshold: Optional[float] - The p-value drop threshold
    :param fold_change_drop_threshold: Optional[float] - The fold change drop threshold
    :param output_path: Optional[str] - The path to save the plot
    :param sub_dir_levels: Optional[List[str]] - The subdirectory levels
    :param file_suffix: str - The suffix for the file name
    :param interactive: bool - Whether to use interactive plots
    :param label_size: int - The size for the labels
    :param max_labels: int - The maximum number of labels
    :param proximity_threshold: int - The proximity threshold
    :param fc_weight: int - The weight for the fold change
    :param title_prefix: Optional[str] - The title prefix
    :return: None
    """
    if output_path is None:
        raise ValueError("output_path is required")
    if sub_dir_levels is None:
        sub_dir_levels = []
    if not isinstance(subpop_cols, list) or len(subpop_cols) == 0:
        raise ValueError("subpop_cols must be a non-empty list of column names")

    required_cols = set(subpop_cols) | {
        drug_family_col,
        feature_col,
        p_value_col,
        fold_change_col,
    }
    missing = required_cols - set(results_df.columns)
    if missing:
        raise KeyError(f"Missing required columns in results_df: {sorted(missing)}")

    # Iterate over each unique subpopulation combination
    unique_combos = results_df[subpop_cols].drop_duplicates()
    for _, combo_row in unique_combos.iterrows():
        # Build mask for this subpopulation
        mask = (results_df[subpop_cols] == combo_row.values).all(axis=1)
        subpop_df = results_df.loc[mask].copy()
        if subpop_df.empty:
            continue

        # Directory: output_path / sub_dir_levels / <subpop values...>
        subpop_values = combo_row.values.astype(str).tolist()
        dir_levels = [str(level) for level in sub_dir_levels] + subpop_values
        dir_path = Path(output_path) / Path(*dir_levels)
        dir_path.mkdir(parents=True, exist_ok=True)

        # For each drug family in this subpopulation
        for drug_family in sorted(subpop_df[drug_family_col].dropna().unique()):
            df_family = subpop_df[subpop_df[drug_family_col] == drug_family].copy()
            if df_family.empty:
                continue

            # Title and file path
            title_parts = []
            if title_prefix:
                title_parts.append(title_prefix)
            title_parts.append("Volcano")
            title_parts.append(f"Drug: {drug_family}")
            title_parts.extend([f"{col}={combo_row[col]}" for col in subpop_cols])
            title = " | ".join(title_parts)

            base_name = f"{sanitize_filename(str(drug_family))}_{file_suffix}"
            
            save_file = dir_path / f"{base_name}.png"
            plot_volcano_plot(
                input_df=df_family,
                p_value_col=p_value_col,
                fold_change_col=fold_change_col,
                p_value_threshold=p_value_threshold,
                p_value_floor=p_value_floor,
                fold_change_threshold=fold_change_threshold,
                p_value_drop_threshold=p_value_drop_threshold,
                fold_change_drop_threshold=fold_change_drop_threshold,
                save_path=str(save_file),
                ignore_significance=False,
                custom_annotation_column=feature_col,
                title=title,
                label_size=label_size,
                max_labels=max_labels,
                proximity_threshold=proximity_threshold,
                fc_weight=fc_weight,
                x_axis_label="Fold Change",
                y_axis_label="-log10(p-value)",
            )