"""Interactive plots of already-calculated biological-sample qPCR results."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from qpcr_analysis import parse_sample_code


PLOT_TYPES = ("시료 점 + 기하평균", "기하평균 막대", "기하평균 선")
FACTOR_COLUMNS = ("ABA 코드", "온도", "실험군")
TEMPERATURE_LABELS = {"상온": "Ambient", "고온": "Heat"}
AXIS_LABELS = {
    "Condition": "Condition",
    "ABA 코드": "ABA dose code",
    "온도": "Temperature",
    "실험군": "Genotype",
    "ABA 농도": "ABA dose (unit unspecified)",
    "염 수준": "Salt level",
}


def prepare_plot_data(sample_results: pd.DataFrame, sample_codes: dict[str, str]) -> pd.DataFrame:
    """Attach factors only when every biological sample has a supported code."""
    plotted = sample_results.copy()
    plotted["시료 코드"] = plotted["Sample"].map(sample_codes)
    factors = plotted["시료 코드"].map(parse_sample_code)
    if factors.notna().all():
        plotted["실험군"] = factors.map(lambda item: {"0": "Col-1", "1": "ost1"}[item[0]])
        plotted["ABA 코드"] = factors.map(lambda item: int(item[1]))
        plotted["온도"] = plotted["시료 코드"].str[-1].map({"X": "상온", "H": "고온"})
    else:
        salt = plotted["시료 코드"].astype(str).str.extract(r"^\s*(\d+(?:\.\d+)?)\s+(0|LOW|HIGH)\s*$", flags=2)
        if salt.notna().all().all():
            plotted["ABA 농도"] = salt[0].astype(float)
            plotted["염 수준"] = salt[1].str.upper()
    return plotted


def ordered_values(values: pd.Series, column: str) -> list:
    unique = set(values.dropna())
    if column in ("ABA 코드", "ABA 농도"):
        return sorted(unique)
    if column == "염 수준":
        return [value for value in ("0", "LOW", "HIGH") if value in unique]
    if column == "온도":
        return [value for value in ("상온", "고온") if value in unique]
    if column == "실험군":
        return [value for value in ("Col-1", "ost1") if value in unique]
    return sorted(unique, key=str)


def make_figure(
    data: pd.DataFrame,
    gene: str,
    chart_type: str,
    x_column: str,
    series_column: str | None = None,
    facet_column: str | None = None,
    log_scale: bool = True,
):
    """Return a figure with biological N for each plotted geometric mean."""
    if chart_type not in PLOT_TYPES:
        raise ValueError("지원하지 않는 그래프 형식입니다.")
    if x_column not in data or (series_column and series_column not in data) or (facet_column and facet_column not in data):
        raise ValueError("선택한 그래프 변수가 데이터에 없습니다.")
    if chart_type == "기하평균 선" and x_column not in ("ABA 코드", "ABA 농도"):
        raise ValueError("선 그래프의 가로축은 ABA 코드만 선택할 수 있습니다.")
    subset = data.loc[data["Gene"] == gene].copy()
    if subset.empty:
        raise ValueError("선택한 유전자에 표시할 시료가 없습니다.")

    facets = ordered_values(subset[facet_column], facet_column) if facet_column else [None]
    series = ordered_values(subset[series_column], series_column) if series_column else [None]
    x_values = ordered_values(subset[x_column], x_column)
    if not x_values:
        raise ValueError("표시할 조건이 없습니다.")
    fig, axes = plt.subplots(1, len(facets), figsize=(max(7, 5 * len(facets)), 5), squeeze=False, sharey=True)
    palette = plt.get_cmap("tab10")

    for facet_index, facet in enumerate(facets):
        ax = axes[0, facet_index]
        panel = subset.loc[subset[facet_column] == facet] if facet_column else subset
        for series_index, series_value in enumerate(series):
            group = panel.loc[panel[series_column] == series_value] if series_column else panel
            color = palette(series_index % 10)
            means_x, means_y = [], []
            for x_index, x_value in enumerate(x_values):
                values = group.loc[group[x_column] == x_value, "Fold_Change"].to_numpy(dtype=float)
                if not len(values):
                    continue
                if not np.isfinite(values).all() or (values <= 0).any():
                    raise ValueError("그래프 값은 양의 유한 숫자여야 합니다.")
                center = x_index + (series_index - (len(series) - 1) / 2) * 0.22
                mean = float(np.exp(np.log(values).mean()))
                means_x.append(center)
                means_y.append(mean)
                if chart_type == "시료 점 + 기하평균":
                    offsets = np.linspace(-0.06, 0.06, len(values)) if len(values) > 1 else [0]
                    ax.scatter(center + np.asarray(offsets), values, color=color, alpha=0.7, s=28)
                    ax.plot([center - 0.10, center + 0.10], [mean, mean], color="black", linewidth=2.5)
                elif chart_type == "기하평균 막대":
                    bottom = 1 if log_scale else 0
                    ax.bar(center, mean - bottom, bottom=bottom, width=0.18, color=color, alpha=0.8)
                ax.annotate(f"n={len(values)}", (center, mean), xytext=(0, 6),
                            textcoords="offset points", ha="center", fontsize=8)
            if chart_type == "기하평균 선" and means_x:
                label = str(TEMPERATURE_LABELS.get(series_value, series_value)) if series_value is not None else None
                ax.plot(means_x, means_y, marker="o", color=color, label=label)
            elif series_value is not None and means_x:
                ax.plot([], [], color=color, marker="o", linestyle="None",
                        label=str(TEMPERATURE_LABELS.get(series_value, series_value)))
        ax.set_xticks(range(len(x_values)), [str(TEMPERATURE_LABELS.get(x, x)) for x in x_values])
        ax.set_xlabel(AXIS_LABELS.get(x_column, x_column))
        if facet is not None:
            ax.set_title(f"{AXIS_LABELS.get(facet_column, facet_column)} = {TEMPERATURE_LABELS.get(facet, facet)}")
        ax.axhline(1, linestyle="--", color="gray", linewidth=1)
        if log_scale:
            ax.set_yscale("log", base=2)
        if facet_index == 0:
            ax.set_ylabel("Fold change (2^-DeltaDeltaCt)")
        if series_column:
            ax.legend(title=AXIS_LABELS.get(series_column, series_column), fontsize=8)
    fig.suptitle(f"{gene} relative expression")
    fig.tight_layout()
    return fig
