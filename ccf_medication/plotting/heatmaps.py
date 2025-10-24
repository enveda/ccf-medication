"""
Functions to build heatmaps and styled count tables

This module provides functions for:
- Creating count heatmaps and displaying styled tables (`count_heatmap`)
- Plotting heatmaps of significant feature counts (`significant_gene_count_heatmap`)
- Rendering grid-style figures for 3-level MultiIndex DataFrames (`render_grid`) aka the pretty formatted heatmap for the paper
"""

from typing import List, Optional

import itables
import numpy as np
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory


from ccf_medication.constants.tables import (
    DRUG_FAMILY_COL,
    FOLD_CHANGE_COLUMN,
    P_VALUE_COLUMN,
)
from ccf_medication.constants.thresholds import (
    ADJ_PVAL_THRESH,
)
from ccf_medication.utils.analysis import get_significant_genes


def count_heatmap(
    results_df: pd.DataFrame,
    subpopulation_cols: List[str],
    drug_family_col: str = DRUG_FAMILY_COL,
    value_column: str = "feature",
    aggfunc="count",
    fill_value=np.nan,
    unique_drug_families: List[str] = None,
    save_path: str = None,
    return_df: bool = False,
) -> None:
    """
    Create a heatmap of value counts per subpopulation and drug family.

    :param results_df: Input DataFrame containing grouping columns and the column
        whose values will be counted.
    :param subpopulation_cols: Columns defining subpopulation groups (used as
        column levels in the pivot table).
    :param drug_family_col: Column that encodes drug families (used as the pivot
        index).
    :param value_column: Column whose values will be aggregated using ``aggfunc``.
    :param aggfunc: Aggregation function applied to ``value_column`` (e.g.,
        ``"count"``, ``"nunique"`` or a callable).
    :param fill_value: Value used to fill missing entries in the pivot table.
    :param unique_drug_families: List of unique drug families to include in the heatmap.
    :param save_path: Path to save the heatmap image.
    :param return_df: Whether to return the dataframe.
    :return: The gene counts dataframe if ``return_df`` is True, otherwise None.
    """

    # Create gene count pivot
    gene_counts_df = results_df.pivot_table(
        index=[drug_family_col],
        columns=subpopulation_cols,
        values=value_column,
        aggfunc=aggfunc,
        fill_value=fill_value,
    )

    # add additional drug families to the heatmap if they are not in the results_df
    if unique_drug_families is not None:
        for drug_family in unique_drug_families:
            if drug_family not in gene_counts_df.index:
                gene_counts_df.loc[drug_family] = np.nan

    # Preserve original itables setting
    original_show_index = getattr(itables.options, "showIndex", True)
    itables.options.showIndex = True

    if save_path is not None:
        gene_counts_df.to_csv(save_path)

    try:
        # Use background gradient for numeric values
        styled_df = gene_counts_df.style.background_gradient(
            cmap="viridis", axis=0
        ).format(precision=0)
        display(styled_df)
    finally:
        # Restore original setting
        itables.options.showIndex = original_show_index

    if return_df:
        return gene_counts_df


def significant_gene_count_heatmap(
    results_df: pd.DataFrame,
    subpopulation_cols: List[str],
    drug_family_col: str = DRUG_FAMILY_COL,
    pval_threshold: float = ADJ_PVAL_THRESH,
    fc_threshold: float = None,
    pval_column: str = P_VALUE_COLUMN,
    fc_column: str = FOLD_CHANGE_COLUMN,
    value_column: str = "feature",
    return_df: bool = False,
) -> Optional[pd.DataFrame]:
    """
    Plot a heatmap of significant feature counts per subpopulation and drug.

    :param results_df: Aggregated results with subpopulation and statistics columns.
    :param subpopulation_cols: Columns defining subpopulation groups.
    :param drug_family_col: Column that encodes drug families.
    :param pval_threshold: Maximum adjusted p-value to be considered significant.
    :param fc_threshold: Minimum absolute fold change to be considered significant.
    :param pval_column: Name of the adjusted p-value column.
    :param fc_column: Name of the fold change/coefficients column.
    :param value_column: Column whose values will be counted per group.
    :param return_df: If True, return the heatmap dataframe.
    :return: The heatmap dataframe when ``return_df`` is True, otherwise None.
    """

    # get all unique drug families
    unique_drug_families = results_df[drug_family_col].unique()

    # Filter for significant genes
    significant_genes = get_significant_genes(
        results_df, pval_threshold, fc_threshold, pval_column, fc_column
    )

    heatmap_df = count_heatmap(
        results_df=significant_genes,
        subpopulation_cols=subpopulation_cols,
        drug_family_col=drug_family_col,
        value_column=value_column,
        aggfunc="count",
        fill_value=0,
        unique_drug_families=unique_drug_families,
        return_df=return_df,
    )

    if heatmap_df is not None:
        return heatmap_df


def render_grid(
    df: pd.DataFrame,
    title: str,
    path: str,
    *,
    show_bottom_labels: bool,
    # Set to None to auto-detect from df.columns
    omics_order=None,  # e.g. ("Transcriptomics", "Proteomics")
    disease_order=None,  # e.g. ("CD", "UC")
    sub_order_map=None,  # {"Transcriptomics": ("colon","small_intestine"), "Proteomics": ("all",)}
    gap_width=0.8,
    big_gap_width=1.6,
    norm=None,
    cmap=None,
    colorbar_label="Sample count",
):
    """
    Render a grid-style figure for a 3-level MultiIndex DataFrame. 
    (This is the really pretty formatted heatmap for the paper)

    Expects columns as (omics_type, disease, sub). If ``omics_order`` or ``disease_order``
    are None, they are inferred from the DataFrame to avoid exact-label mismatch.

    :param df: DataFrame to plot.
    :param title: Title of the plot.
    :param path: Path to save the plot.
    :param show_bottom_labels: Whether to show the bottom sub labels.
    :param omics_order: Optional order of the omics.
    :param disease_order: Optional order of the diseases.
    :param sub_order_map: Optional mapping from omics -> tuple of sub labels.
    :param gap_width: Width of the gap between groups of the same omics.
    :param big_gap_width: Width of the gap between different omics groups.
    :param norm: Optional normalization for the colormap.
    :param cmap: Optional colormap.
    :param colorbar_label: Label for the colorbar.
    :return: None.
    :raises ValueError: If ``df.columns`` is not a 3-level MultiIndex or if no matching
        columns are found for the requested layout.
    """
    # ---- validate ----
    if not isinstance(df.columns, pd.MultiIndex) or df.columns.nlevels != 3:
        raise ValueError(
            "Expected df.columns as a 3-level MultiIndex: (omics_type, disease, sub)"
        )

    df = df.copy()

    # ---- auto-detect orders if not provided ----
    lvl0, lvl1, lvl2 = (df.columns.get_level_values(i) for i in range(3))
    if omics_order is None:
        # preserve first-seen order
        omics_order = tuple(dict.fromkeys(lvl0))
    if disease_order is None:
        disease_order = tuple(dict.fromkeys(lvl1))

    # ---- defaults ----
    if cmap is None:
        cmap = plt.get_cmap("Purples")

    present = set(df.columns.to_list())

    def subs_for(omics_name, disease_name):
        """Return ordered list of (omics, disease, sub) columns present in df."""
        if sub_order_map and omics_name in sub_order_map:
            wanted = list(sub_order_map[omics_name])
            cols = [
                (omics_name, disease_name, s)
                for s in wanted
                if (omics_name, disease_name, s) in present
            ]
            # include any extra subs in native df order
            extras = [
                c
                for c in df.columns
                if c[0] == omics_name and c[1] == disease_name and c not in cols
            ]
            cols.extend(extras)
            return cols
        else:
            # whatever subs exist under (omics, disease), in df order
            return [
                c for c in df.columns if c[0] == omics_name and c[1] == disease_name
            ]

    # ---- build layout groups ----
    ordered_cols, group_sizes, group_omics, group_labels = [], [], [], []
    for omics_name in omics_order:
        for disease_name in disease_order:
            cols = subs_for(omics_name, disease_name)
            if len(cols) == 0:
                continue
            ordered_cols.extend(cols)
            group_sizes.append(len(cols))
            group_omics.append(omics_name)
            group_labels.append((omics_name, disease_name))

    if len(ordered_cols) == 0:
        # Helpful diagnostics
        raise ValueError(
            "No matching columns found for the requested layout.\n"
            f"Found omics={sorted(set(lvl0))}\n"
            f"Found disease={sorted(set(lvl1))}\n"
            f"Example columns (first 6): {list(df.columns[:6])}"
        )

    # ---- reorder & norm ----
    df = df[ordered_cols]
    values = df.to_numpy(dtype=float)
    nrows, ncols = values.shape

    if norm is None:
        vmin = float(np.nanmin(values)) if np.isfinite(values).any() else 0.0
        vmax = float(np.nanmax(values)) if np.isfinite(values).any() else 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # ---- augmented grid with gaps ----
    gap_after = []
    for gi in range(len(group_sizes) - 1):
        gap_after.append(
            big_gap_width if group_omics[gi] != group_omics[gi + 1] else gap_width
        )

    aug_cols, aug_to_real, widths = [], [], []
    real_idx = 0
    for gi, size in enumerate(group_sizes):
        for k in range(size):
            aug_cols.append(df.iloc[:, real_idx + k].to_numpy(dtype=float))
            aug_to_real.append(real_idx + k)
            widths.append(1.0)
        real_idx += size
        if gi < len(group_sizes) - 1:
            aug_cols.append(np.full(nrows, np.nan))
            aug_to_real.append(None)
            widths.append(gap_after[gi])

    aug = np.ma.masked_invalid(np.column_stack(aug_cols).astype(float))
    x_edges = np.concatenate([[0.0], np.cumsum(widths)])
    y_edges = np.arange(nrows + 1)

    spans, xc = [], 0.0
    for gi, size in enumerate(group_sizes):
        x0 = xc
        xc += size
        x1 = xc
        spans.append((x0, x1))
        if gi < len(group_sizes) - 1:
            xc += gap_after[gi]

    # super-spans for omics headers
    super_spans = []
    i = 0
    while i < len(group_sizes):
        j = i
        cur_omics = group_omics[i]
        x0 = spans[i][0]
        while j + 1 < len(group_sizes) and group_omics[j + 1] == cur_omics:
            j += 1
        x1 = spans[j][1]
        super_spans.append((x0, x1, cur_omics))
        i = j + 1

    # ---- draw ----
    fig, ax = plt.subplots(figsize=(11.0, 4.2))
    ax.grid(False)
    pm = ax.pcolormesh(
        x_edges,
        y_edges,
        aug,
        shading="flat",
        edgecolors="white",
        linewidth=1,
        norm=norm,
        cmap=cmap,
        zorder=2,
    )
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    outer_gap = 0.30
    ax.set_xlim(x_edges[0] - outer_gap, x_edges[-1] + outer_gap)
    ax.set_ylim(nrows + outer_gap, -outer_gap)

    rect_clear = 0.085
    for x0, x1 in spans:
        ax.add_patch(
            Rectangle(
                (x0 - rect_clear, -rect_clear),
                (x1 - x0) + 2 * rect_clear,
                nrows + 2 * rect_clear,
                fill=False,
                linewidth=2.2,
                zorder=6,
                color="black",
            )
        )

    # numbers with auto-contrast
    for j_aug in range(aug.shape[1]):
        rc = aug_to_real[j_aug]
        if rc is None:
            continue
        xc = 0.5 * (x_edges[j_aug] + x_edges[j_aug + 1])
        for irow in range(nrows):
            v = df.iat[irow, rc]
            if pd.isna(v):
                label, color = "NA", "#444"
            else:
                fv = float(v)
                label = f"{int(fv)}" if fv.is_integer() else f"{fv:.0f}"
                rC, gC, bC, _ = cmap(norm(fv))
                lum = 0.299 * rC + 0.587 * gC + 0.114 * bC
                color = "white" if lum < 0.52 else "black"
            ax.text(
                xc,
                irow + 0.5,
                label,
                ha="center",
                va="center",
                fontsize=7.5,
                zorder=7,
                color=color,
            )

    # headers
    trans_top = blended_transform_factory(ax.transData, ax.transAxes)
    disease_line_y = 1.008
    disease_text_y = 1.038
    ax.plot(
        [x_edges[0], x_edges[-1]],
        [disease_line_y, disease_line_y],
        transform=trans_top,
        linewidth=2,
        zorder=8,
    )

    # disease labels per box
    start_idx = 0
    for (x0, x1), (_omics, disease_name) in zip(spans, group_labels):
        ax.text(
            0.5 * (x0 + x1),
            disease_text_y,
            str(disease_name).upper(),
            transform=trans_top,
            ha="center",
            va="center",
            fontsize=11,
            zorder=9,
            color="black",
        )
        start_idx += size

    # omics super-headers
    omics_text_y = 1.125
    for x0, x1, omics_name in super_spans:
        ax.text(
            0.5 * (x0 + x1),
            omics_text_y,
            str(omics_name),
            transform=trans_top,
            ha="center",
            va="center",
            fontsize=13,
            zorder=9,
            color="black",
        )

    # bottom labels
    if show_bottom_labels:
        trans_bottom = blended_transform_factory(ax.transData, ax.transAxes)
        bottom_line_y = -0.09
        bottom_text_y = 0.0
        ax.plot(
            [x_edges[0], x_edges[-1]],
            [bottom_line_y, bottom_line_y],
            transform=trans_bottom,
            linewidth=2,
            zorder=8,
        )
        j_real = 0
        for j_aug in range(aug.shape[1]):
            rc = aug_to_real[j_aug]
            if rc is None:
                continue
            sublab = df.columns[j_real][2]
            xc = 0.5 * (x_edges[j_aug] + x_edges[j_aug + 1])
            ax.text(
                xc,
                bottom_text_y,
                str(sublab),
                transform=trans_bottom,
                ha="right",
                va="top",
                fontsize=8,
                zorder=9,
                color="black",
                rotation=45,
            )
            j_real += 1

    # row labels
    trans_left = blended_transform_factory(ax.transAxes, ax.transData)
    left_label_x = -0.010
    for irow, name in enumerate(df.index.tolist()):
        ax.text(
            left_label_x,
            irow + 0.5,
            str(name),
            transform=trans_left,
            ha="right",
            va="center",
            fontsize=8,
            clip_on=False,
            zorder=9,
            color="black",
        )

    ax.set_title(title, fontsize=16, pad=36)
    cb = fig.colorbar(pm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(colorbar_label)

    fig.subplots_adjust(
        left=0.12, right=0.96, top=0.80, bottom=0.26 if show_bottom_labels else 0.18
    )
    fig.savefig(path, dpi=400, bbox_inches="tight")
    plt.close(fig)

