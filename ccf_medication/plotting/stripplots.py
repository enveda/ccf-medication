"""
Functions for creating strip + box plots

This module provides functions for:
- Plotting strip + box plots with optional significance brackets (`plot_strip_box_plot`)
"""
from typing import Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, PathPatch
import numpy as np
import pandas as pd
import seaborn as sns


def plot_strip_box_plot(
    df: pd.DataFrame,
    x: str,
    y: str,
    # categories/palette
    all_categories: Optional[list] = None,            # order of categories to show (original x values or severity ints)
    palette_by_category: Optional[dict] = None,       # map original categories (or severity ints) -> color hex
    # mode and labels
    x_mode: str = "category",                         # "category" or "severity"
    severity_labels: dict[int, str] = {0: "Remission", 1: "Active"},
    label_map: Optional[dict] = None,                 # map original categories -> display labels (e.g., drug-family prettifier)
    # annotations
    sig_pairs: Optional[list[tuple]] = None,          # list of (left_cat, right_cat) using original category keys
    sig_text: str = "*",
    # styling
    title: str = "",
    y_label: str = "",
    legend_title: Optional[str] = None,               # default depends on mode
    legend_facecolor: bool = False,
    figsize=(9, 6),
    dpi: int = 400,
    title_fontsize: int = 14,
    xlabel_fontsize: int = 10,
    ylabel_fontsize: int = 12,
    legend_fontsize: int = 9,
    legend_title_fontsize: int = 12,
    x_label_rotation: int = 35,
    x_label_ha: str = "right",
    alpha: float = 0.95,
):
    """
    Plot a strip plot with boxplot w/optional significance brackets.

    :param df: DataFrame containing the data to plot.
    :param x: Column name for the x-axis.
    :param y: Column name for the y-axis.
    :param all_categories: List of categories to plot (original x values or severity ints).
    :param palette_by_category: Dictionary mapping categories to colors.
    :param x_mode: Mode for the x-axis ("category" or "severity").
    :param severity_labels: Dictionary mapping severity integers to labels.
    :param label_map: Dictionary mapping original categories to display labels.
    :param sig_pairs: List of (left_cat, right_cat) tuples for significance brackets. 
    :param sig_text: Text for significance brackets.
    :param title: Title of the plot.
    :param y_label: Label for the y-axis.
    :param legend_title: Title of the legend.
    :param legend_facecolor: Whether to fill the legend with the category colors.
    :param figsize: Size of the figure.
    :param dpi: DPI of the figure.
    :param title_fontsize: Font size of the title.
    :param xlabel_fontsize: Font size of the x-axis labels.
    :param ylabel_fontsize: Font size of the y-axis labels.
    :param legend_fontsize: Font size of the legend.
    :param legend_title_fontsize: Font size of the legend title.
    :param x_label_rotation: Rotation of the x-axis labels.
    :param x_label_ha: Horizontal alignment of the x-axis labels.
    :param alpha: Alpha value for the strip plot.
    :return: Figure and axis objects.
    """

    sns.set(style="ticks", context="talk")

    dfp = df.copy()
    # numeric y
    dfp[y] = pd.to_numeric(dfp[y], errors="coerce")
    dfp = dfp[np.isfinite(dfp[y]) & dfp[y].notna()]

    # set defaults
    if legend_title is None:
        legend_title = "Disease Severity" if x_mode == "severity" else "Medication Group"

    # build categories, display labels, and palette
    if x_mode == "severity":
        # coerce x to int 0/1 and map to labels
        dfp["_sev_int"] = pd.to_numeric(dfp[x], errors="coerce").astype("Int64")
        valid = dfp["_sev_int"].isin([0, 1])
        if all_categories is None:
            all_categories = [c for c in [0, 1] if c in set(dfp.loc[valid, "_sev_int"].dropna().astype(int).unique())]
        present = [c for c in all_categories if c in set(dfp.loc[valid, "_sev_int"].dropna().astype(int).unique())]
        if len(present) == 0:
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
            ax.set_visible(False)
            print("No categories with valid data to plot.")
            return fig, ax

        dfp = dfp[valid].copy()
        dfp["_x_disp"] = dfp["_sev_int"].astype(int).map(severity_labels)

        # palette (keys are ints in present)
        if palette_by_category is None:
            colors = sns.color_palette("muted", n_colors=len(present)).as_hex()
            palette_by_category = {c: col for c, col in zip(present, colors)}
        order_disp = [severity_labels[c] for c in present]
        palette_disp = {severity_labels[c]: palette_by_category[c] for c in present}

        # for significance handling, map original (ints) -> display labels
        category_to_disp = {c: severity_labels[c] for c in present}

    else:  # category mode
        # choose categories to plot and order
        uniques = dfp[x].dropna().unique().tolist()
        if all_categories is None:
            all_categories = sorted(uniques, key=lambda v: str(v))
        present = [c for c in all_categories if c in set(uniques)]
        if len(present) == 0:
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
            ax.set_visible(False)
            print("No categories with valid data to plot.")
            return fig, ax

        # display labels via label_map if provided
        def _disp(c):
            return label_map.get(c, c) if label_map else c

        dfp["_x_disp"] = dfp[x].map(lambda c: _disp(c))
        order_disp = [_disp(c) for c in present]

        # palette (keys are original categories)
        if palette_by_category is None:
            colors = sns.color_palette("muted", n_colors=len(present)).as_hex()
            palette_by_category = {c: col for c, col in zip(present, colors)}
        palette_disp = {_disp(c): palette_by_category[c] for c in present}
        category_to_disp = {c: _disp(c) for c in present}

    # plot
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # boxplot (single call) then recolor edges/whiskers per category
    sns.boxplot(
        data=dfp, x="_x_disp", y=y, order=order_disp,
        hue="_x_disp", hue_order=order_disp, dodge=False,
        palette=palette_disp,
        ax=ax, width=0.55, fliersize=0, boxprops=dict(facecolor="none"),
        whiskerprops=dict(linewidth=1.4), capprops=dict(linewidth=1.4),
        medianprops=dict(linewidth=1.6), linewidth=1.6,
    )

    # recolor edges explicitly (robust to seaborn/mpl versions)
    try:
        boxes_artists = list(getattr(ax, "artists", []))
        boxes_patches = [p for p in ax.patches if isinstance(p, PathPatch)]
        boxes = boxes_patches if len(boxes_patches) > 0 else boxes_artists
        for i, lbl in enumerate(order_disp):
            color = palette_disp[lbl]
            if i < len(boxes):
                try:
                    boxes[i].set_edgecolor(color)
                    boxes[i].set_facecolor("none")
                    boxes[i].set_linewidth(1.6)
                except Exception:
                    pass
            # 6 lines per box (2 whiskers, 2 caps, 1 median, 1 box median line variant)
            line_start = i * 6
            line_end = min(line_start + 6, len(ax.lines))
            for j in range(line_start, line_end):
                try:
                    ax.lines[j].set_color(color)
                    ax.lines[j].set_mfc(color)
                    ax.lines[j].set_mec(color)
                except Exception:
                    pass
    except Exception:
        pass

    # points
    sns.stripplot(
        data=dfp, x="_x_disp", y=y, order=order_disp,
        hue="_x_disp", hue_order=order_disp, dodge=False,
        palette=palette_disp,
        ax=ax, size=4, jitter=0.18, linewidth=0.4, edgecolor="none", alpha=alpha,
    )

    # significance brackets (if provided)
    if sig_pairs:
        y_min, y_max = ax.get_ylim()
        y_range = y_max - y_min
        current_top = y_max
        step = 0.08 * y_range

        # keep only pairs present; convert to display labels
        pairs_disp = []
        for lc, rc in sig_pairs:
            if lc in category_to_disp and rc in category_to_disp:
                ldisp, rdisp = category_to_disp[lc], category_to_disp[rc]
                if ldisp in order_disp and rdisp in order_disp:
                    x1_tmp, x2_tmp = order_disp.index(ldisp), order_disp.index(rdisp)
                    pairs_disp.append(((ldisp, rdisp), min(x1_tmp, x2_tmp)))
        pairs_disp.sort(key=lambda t: t[1])

        for i, ((ld, rd), _) in enumerate(pairs_disp):
            x1, x2 = order_disp.index(ld), order_disp.index(rd)
            yb = current_top - i * step
            h = 0.03 * y_range
            ax.plot([x1, x1, x2, x2], [yb, yb + h, yb + h, yb], lw=1.4, c="black")
            ax.text((x1 + x2) / 2, yb + h + 0.01 * y_range, sig_text,
                    ha="center", va="bottom", fontsize=12, color="black")
        ax.set_ylim(y_min, current_top + 0.12 * y_range)

    # styling
    ax.set_title(title, pad=8, fontsize=title_fontsize)
    ax.set_xlabel(""); ax.set_ylabel("")
    # set ticks and labels without FixedFormatter warning
    ax.set_xticks(range(len(order_disp)), labels=order_disp)
    if x_mode == "category":
        ax.tick_params(axis="x", labelrotation=x_label_rotation, labelsize=xlabel_fontsize)
        for lbl in ax.get_xticklabels():
            lbl.set_ha(x_label_ha)
    else:
        ax.tick_params(axis="x", labelrotation=0, labelsize=xlabel_fontsize)
        for lbl in ax.get_xticklabels():
            lbl.set_ha("center")
    ax.set_ylabel(y_label, fontsize=ylabel_fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(ylabel_fontsize)
    sns.despine(ax=ax, left=False, bottom=False, right=True, top=True)

    # remove any seaborn auto-legend before adding custom legend
    if ax.get_legend() is not None:
        try:
            ax.legend_.remove()
        except Exception:
            pass

    # legend with counts
    counts = dfp["_x_disp"].value_counts().to_dict()
    legend_labels = [f"{lbl} (n={counts.get(lbl, 0)})" for lbl in order_disp]
    if legend_facecolor:
        handles = [Patch(facecolor=palette_disp[lbl], edgecolor=palette_disp[lbl], linewidth=1.6, label=legend_labels[i])
                   for i, lbl in enumerate(order_disp)]
    else:
        handles = [Patch(facecolor="none", edgecolor=palette_disp[lbl], linewidth=1.0, label=legend_labels[i])
                   for i, lbl in enumerate(order_disp)]
    ax.legend(handles=handles, labels=legend_labels, title=legend_title, frameon=False,
              loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0,
              fontsize=legend_fontsize, title_fontsize=legend_title_fontsize)

    ax.tick_params(axis="both", which="both", direction="in")
    plt.tight_layout()
    return fig, ax
