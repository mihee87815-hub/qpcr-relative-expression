"""Evidence-linked, descriptive report from calculated qPCR results."""

from __future__ import annotations

import math

import pandas as pd

from qpcr_analysis import AnalysisResult, parse_sample_code


def _fmt(value: float) -> str:
    return f"{value:.3g}"


def _direction(fc: float) -> str:
    return "높게" if fc > 1 else "낮게" if fc < 1 else "같게"


def write_report(
    result: AnalysisResult,
    reference_gene: str,
    control_condition: str,
    title: str = "qPCR 상대 발현량 결과 보고서",
    experiment_note: str = "",
    sample_codes: dict[str, str] | None = None,
) -> str:
    """Write observations and bounded interpretations without significance claims."""
    summary = result.summary.sort_values(["Gene", "Condition"])
    qc = result.technical_qc
    lines = [f"# {title.strip() or 'qPCR 상대 발현량 결과 보고서'}", ""]
    if experiment_note.strip():
        lines += ["## 실험 정보", experiment_note.strip(), ""]
    lines += [
        "## 분석 방법과 해석 범위",
        f"기준 유전자는 {reference_gene}, 보정 대조군은 {control_condition}입니다. "
        "시료·유전자별 기술 반복 Ct 평균으로 ΔCt(표적 Ct − 기준 Ct)를 구하고, "
        "대조군의 평균 ΔCt를 빼서 ΔΔCt를 계산했습니다. 상대 발현량은 2^-ΔΔCt이며, "
        "조건별 값은 독립 생물학적 시료의 기하평균입니다.",
        "아래 수치는 대조군에 대한 기술적·탐색적 비교입니다. 통계 검정은 수행하지 않았으며, "
        "수치 차이만으로 유의성, 원인 또는 분자 기작을 확정할 수 없습니다.",
        "",
        "## 주요 관찰과 근거",
    ]
    noncontrol = summary[summary["Condition"] != control_condition]
    if noncontrol.empty:
        lines.append("대조군 이외의 분석 가능 조건이 없어 처리군 비교를 작성할 수 없습니다.")
    for gene, frame in noncontrol.groupby("Gene", sort=True):
        lines.append(f"### {gene}")
        for row in frame.itertuples(index=False):
            fc = float(row.Geometric_Mean_FC)
            n = int(row.Biological_N)
            lines.append(
                f"- {row.Condition}: 대조군 대비 {_fmt(fc)}배로 관찰되었습니다 "
                f"(생물학적 시료 n={n}, 평균 ΔCt={_fmt(float(row.Mean_Delta_Ct))}). "
                f"이는 이 데이터의 정규화 발현량이 대조군보다 {_direction(fc)} 측정됐다는 뜻입니다."
            )
        lines.append("")

    # Compare only matched genotype/dose or temperature combinations. The ratio
    # of two fold changes is a descriptive within-gene contrast.
    if sample_codes:
        coded = result.sample_results[["Sample", "Condition"]].drop_duplicates().copy()
        coded["code"] = coded["Sample"].map(sample_codes)
        coded["parts"] = coded["code"].map(parse_sample_code)
        if coded["parts"].notna().all() and coded.groupby("Condition")["code"].nunique().eq(1).all():
            coded["genotype"] = coded["parts"].map(lambda p: {"0": "Col-1", "1": "ost1"}[p[0]])
            coded["dose"] = coded["parts"].map(lambda p: int(p[1]))
            coded["temperature"] = coded["code"].str[-1].map({"X": "상온", "H": "고온"})
            metadata = coded.groupby("Condition")[["genotype", "dose", "temperature"]].first().reset_index()
            matched = summary.merge(metadata, on="Condition", how="inner")
            comparisons = []
            for (gene, genotype, dose), group in matched.groupby(["Gene", "genotype", "dose"]):
                if set(group["temperature"]) == {"상온", "고온"}:
                    ambient = group[group["temperature"] == "상온"].iloc[0]
                    heat = group[group["temperature"] == "고온"].iloc[0]
                    ratio = heat.Geometric_Mean_FC / ambient.Geometric_Mean_FC
                    if math.isfinite(ratio):
                        comparisons.append(
                            f"- {gene}, {genotype}, ABA 코드 {dose}: 고온/상온 상대 발현량 비는 {_fmt(ratio)} "
                            f"({heat.Condition} {_fmt(heat.Geometric_Mean_FC)}배, "
                            f"{ambient.Condition} {_fmt(ambient.Geometric_Mean_FC)}배; "
                            f"각 n={int(heat.Biological_N)}, {int(ambient.Biological_N)}). "
                            "동일 유전자·실험군·ABA 코드 내 조건부 관찰이며, 고온의 인과 효과를 입증하지는 않습니다."
                        )
            if comparisons:
                lines += ["## 같은 실험군·ABA 코드에서 온도 비교", *comparisons, ""]
            lines += ["시료 코드의 X는 상온, H는 고온, 첫 숫자 0은 Col-1, 1은 ost1로 해석했습니다. "
                      "ABA 숫자는 처리 코드이며 실제 농도와 단위는 실험 기록에서 확인해야 합니다.", ""]

    one = int((qc["Technical_N"] == 1).sum())
    high = int((qc["QC_Status"] == "SD 기준 초과").sum())
    min_n = int(summary["Biological_N"].min()) if not summary.empty else 0
    lines += [
        "## 품질 관리와 제한 사항",
        f"- 기술 반복 1회 묶음: {one}개. 이 묶음에서는 기술 반복 SD를 추정할 수 없습니다.",
        f"- 설정된 SD 기준 초과 묶음: {high}개. 해당 웰의 증폭곡선과 원자료를 확인해야 합니다.",
        f"- 조건·유전자별 최소 생물학적 시료 수: n={min_n}. "
        + ("독립 생물학적 반복이 부족해 조건 차이의 재현성이나 통계적 유의성을 평가할 수 없습니다." if min_n < 3
           else "생물학적 반복이 있어도 이 보고서는 통계 검정을 수행하지 않았습니다."),
        "- 기준 유전자 안정성, 프라이머 증폭 효율, 용융곡선, 제외 웰의 원인과 배치 효과는 "
        "이 계산만으로 검증할 수 없습니다. 해석 전 원본 실험 기록을 확인해야 합니다.",
        "",
        "## 결론",
        "위 조건별 발현량 차이와 일치하는 변화가 이 자료에서 관찰되었습니다. "
        "현재 근거는 정규화된 Ct에 기반한 탐색적 기술 결과이므로, "
        "생물학적 반복과 독립 검증을 확보한 뒤 효과의 재현성 및 기작을 판단해야 합니다.",
        "",
        "## 조건별 근거표",
        "| 유전자 | 조건 | 생물학적 n | 평균 ΔCt | 상대 발현량 기하평균 |",
        "|---|---|---:|---:|---:|",
    ]
    for row in summary.itertuples(index=False):
        lines.append(f"| {row.Gene} | {row.Condition} | {int(row.Biological_N)} | "
                     f"{_fmt(float(row.Mean_Delta_Ct))} | {_fmt(float(row.Geometric_Mean_FC))} |")
    return "\n".join(lines) + "\n"
