"""
Functions for plotting categorical panels.

This module provides functions for:
- Plotting categorical panels (`plot_categorical_panels`)
"""

# Imports
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps as mpl_cmaps
from matplotlib import colors as mcolors

def plot_categorical_panels(
    df: pd.DataFrame,
    columns: list[str],
    hue: str | None = None,
    *,
    ncols: int = 3,
    figsize: tuple[float, float] | None = None,
    max_categories: int | None = None,
    sort: str = "count",    # "count" or "label"
    include_missing: bool = True,
    palette: str | list[str] = "Set2",  # string (cmap name) OR list of hex/rgb strings
    tight_layout: bool = True,
    bar_alpha: float = 0.5,
    dpi: int = 400,
):
    """
    One subplot per column in `columns`.
    - Each category within a subplot gets a color from `palette`.
      * If palette is a string, uses that Matplotlib colormap (e.g., "Set2").
      * If palette is a list of hex/rgb strings, cycles through that list.
    - If `hue` is provided, each category's bar is split into multiple bars
      with the SAME color and distinct hatches per hue level.

    :param df: DataFrame to plot.
    :param columns: List of columns to plot.
    :param hue: Column to use for the hue.
    :param ncols: Number of columns per row.
    :param figsize: Size of the figure.
    :param max_categories: Maximum number of categories to plot.
    :param sort: Sort the categories by count or label.
    :param include_missing: Whether to include missing values.
    :param palette: Palette to use for the plot.
    :param tight_layout: Whether to use tight layout.
    :param bar_alpha: Alpha for the bars.
    :param dpi: DPI for the plot.
    :return: Figure and axes.
    """
    df = df.copy()
    MISSING = "Missing"

    def _norm(s: pd.Series):
        s = s.astype("object")
        return s.fillna(MISSING) if include_missing else s.dropna()

    # hue setup
    if hue is not None and hue in df.columns:
        hue_values = _norm(df[hue])
        hue_levels = list(pd.Index(hue_values.unique()).dropna())
        if MISSING in hue_levels:
            hue_levels = [lvl for lvl in hue_levels if lvl != MISSING] + [MISSING]
        hatch_cycle = ["", "///", "\\\\\\", "xxx", "+++", "...", "***", "ooo"]
        while len(hatch_cycle) < len(hue_levels):
            hatch_cycle += hatch_cycle
    else:
        hue_levels, hatch_cycle, hue = [], [], None

    # figure layout
    n = len(columns)
    ncols = max(1, int(ncols))
    nrows = math.ceil(n / ncols)
    if figsize is None:
        figsize = (ncols * 5.0, nrows * 3.8)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, dpi=dpi)
    axes_flat = axes.ravel()
    for ax in axes_flat:
        ax.axis("off")

    # palette resolver (supports cmap name OR list of colors)
    def _category_color_map(categories: list[str]):
        if isinstance(palette, (list, tuple)):
            if len(palette) == 0:
                raise ValueError("Custom palette list is empty.")
            # normalize to RGBA
            pal = [mcolors.to_rgba(c) for c in palette]
            return {cat: pal[i % len(pal)] for i, cat in enumerate(categories)}
        else:
            cmap = mpl_cmaps.get_cmap(palette)  
            n_colors = getattr(cmap, "N", 8)
            return {cat: cmap(i % n_colors) for i, cat in enumerate(categories)}

    legend_created = False

    for i, col in enumerate(columns):
        ax = axes_flat[i]
        ax.axis("on")

        if col not in df.columns:
            ax.set_title(f"{col} (not in DataFrame)")
            continue

        col_series = _norm(df[col])
        vc = col_series.value_counts()

        if max_categories is not None and len(vc) > max_categories:
            top = vc.iloc[:max_categories]
            other_count = vc.iloc[max_categories:].sum()
            vc = pd.concat([
                top,
                pd.Series({f"Other ({len(col_series.unique())-max_categories})": other_count})
            ])

        if sort == "label":
            order = sorted(vc.index.astype(str), key=lambda x: (x == MISSING, x))
        else:
            order = list(vc.sort_values(ascending=False).index)

        x = np.arange(len(order))
        width = 0.8

        # color per category (from cmap or custom list)
        cat_colors = _category_color_map(order)

        max_bar = 0
        if hue is None:
            counts = vc.reindex(order).fillna(0)
            for j, (cat, count) in enumerate(counts.items()):
                bar = ax.bar(x[j], count, width=width, color=cat_colors[cat],
                             edgecolor="black", linewidth=0.5, alpha=bar_alpha)
                max_bar = max(max_bar, int(count))
                if count > 0:
                    ax.text(bar[0].get_x() + bar[0].get_width()/2, count, f"{int(count)}",
                            ha="center", va="bottom", fontsize=9)
        else:
            tmp = pd.DataFrame({col: _norm(df[col]), hue: _norm(df[hue])})
            ct = pd.crosstab(tmp[col], tmp[hue]).reindex(index=order, columns=hue_levels, fill_value=0)

            nh = len(hue_levels)
            inner_w = width / max(nh, 1)
            offsets = (np.arange(nh) - (nh - 1) / 2.0) * inner_w

            for j, cat in enumerate(order):
                for k, lvl in enumerate(hue_levels):
                    height = int(ct.loc[cat, lvl])
                    bars = ax.bar(x[j] + offsets[k], height, width=inner_w * 0.95,
                                  color=cat_colors[cat], edgecolor="black",
                                  linewidth=0.6, hatch=hatch_cycle[k], alpha=bar_alpha)
                    max_bar = max(max_bar, height)
                    if height > 0:
                        ax.text(bars[0].get_x() + bars[0].get_width()/2, height, f"{height}",
                                ha="center", va="bottom", fontsize=8)

        ax.set_title(_pretty(col), fontweight="bold")

        def _pretty_labels(var, col) -> pd.DataFrame:
            if col not in ["batch", "diagnosis", 'drug_family']:
                return var.title()
            elif col == "diagnosis":
                return var.upper()
            else:
                return var

        ax.set_xticks(x)
        labels = [str(_pretty_labels(v, col)) for v in order]
        rotate = any(len(lbl) > 15 for lbl in labels)
        ax.set_xticklabels(labels, rotation=(20 if rotate else 0), ha=("right" if rotate else "center"))

        # y-axis ticks: if max < 600 -> [0, 250, 500]; else -> [0, 250, 500, 750, 1000]
        if max_bar < 600:
            ticks = np.array([0, 250, 500])
            top_tick = 500
        else:
            ticks = np.array([0, 250, 500, 750, 1000])
            top_tick = 1000
        ax.set_yticks(ticks)
        # Add headroom if bars exceed the last tick
        axis_top = max(top_tick, int(math.ceil(max_bar * 1.1))) if max_bar > 0 else top_tick
        ax.set_ylim(0, axis_top)

        # Only show y-axis label on first column of subplots
        col_idx = i % ncols
        if col_idx == 0:
            ax.set_ylabel("Count")
        else:
            ax.set_ylabel("")

        ax.grid(axis="y", linestyle=":", alpha=0.5)
        ax.set_axisbelow(True)

    # hide unused axes
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")

    # Add a single, figure-level legend above the subplots to avoid covering any panel
    if hue is not None and len(hue_levels) > 0:
        neutral = (0.92, 0.92, 0.92, 1.0)
        handles, labels = [], []
        for k, lvl in enumerate(hue_levels):
            proxy = plt.Rectangle((0, 0), 1, 1, facecolor=neutral,
                                  edgecolor="black", linewidth=0.6,
                                  hatch=hatch_cycle[k])
            handles.append(proxy)
            labels.append(str(lvl))
        fig.legend(handles, [l.capitalize() for l in labels], title=hue.replace("_", " ").capitalize(), fontsize=9, title_fontsize=10,
                   frameon=True, loc="upper center", bbox_to_anchor=(0.5, 1.02),
                   ncol=min(len(hue_levels), 5))
        legend_created = True

    if tight_layout:
        if legend_created:
            fig.tight_layout(rect=[0, 0, 1, 0.92])  # leave room for legend
        else:
            fig.tight_layout()
    return fig, axes

def _pretty(s: str) -> str:
    return " ".join([w if w.isupper() else w.capitalize() for w in s.replace("_", " ").split()])