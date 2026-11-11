"""Deterministic relative-expression analysis for technical-replicate qPCR CSVs with outlier exclusion and statistical tests."""

from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
import math

import numpy as np
import pandas as pd
from scipy import stats


REQUIRED = {"Sample", "Condition", "Gene", "Replicate", "Ct"}


class QPCRInputError(ValueError):
    pass


@dataclass
class AnalysisResult:
    technical_qc: pd.DataFrame
    sample_results: pd.DataFrame
    summary: pd.DataFrame
    reference_summary: pd.DataFrame
    statistical_tests: pd.DataFrame
    warnings: list[str]


def read_qpcr_csv(content: bytes) -> pd.DataFrame:
    try:
        return pd.read_csv(BytesIO(content), encoding="utf-8-sig")
    except (UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise QPCRInputError(f"CSV 파일을 읽을 수 없습니다: {exc}") from exc


def analyze(
    data: pd.DataFrame,
    reference_gene: str,
    control_condition: str,
    max_technical_sd: float = 0.3,
    min_technical_replicates: int = 2,
    auto_exclude_outliers: bool = True
) -> AnalysisResult:
    missing = REQUIRED - set(data.columns)
    if missing:
        raise QPCRInputError(f"필수 열이 없습니다: {', '.join(sorted(missing))}")
    if not reference_gene or not control_condition:
        raise QPCRInputError("기준 유전자와 대조군을 선택하세요.")
    if not math.isfinite(max_technical_sd) or max_technical_sd <= 0:
        raise QPCRInputError("기술 반복 SD 기준은 0보다 커야 합니다.")
    if min_technical_replicates < 2:
        raise QPCRInputError("기술 반복 최소 횟수는 2 이상이어야 합니다.")

    df = data[list(REQUIRED)].copy()
    for col in ("Sample", "Condition", "Gene", "Replicate"):
        if df[col].isna().any() or df[col].astype(str).str.strip().eq("").any():
            raise QPCRInputError(f"{col} 열에 빈 값이 있습니다.")
        df[col] = df[col].astype(str).str.strip()
    df["Ct"] = pd.to_numeric(df["Ct"], errors="coerce")
    if not np.isfinite(df["Ct"].to_numpy(dtype=float)).all() or (df["Ct"] <= 0).any():
        raise QPCRInputError("Ct는 비어 있지 않은 양의 유한 숫자여야 합니다. Undetermined 값은 먼저 검토하세요.")
    if df.duplicated(["Sample", "Condition", "Gene", "Replicate"]).any():
        raise QPCRInputError("동일 시료·조건·유전자에서 Replicate 값이 중복됩니다.")
    if df.groupby("Sample")["Condition"].nunique().gt(1).any():
        raise QPCRInputError("한 Sample이 여러 Condition에 속합니다. 생물학적 시료 ID를 구분하세요.")
    if reference_gene not in set(df["Gene"]):
        raise QPCRInputError(f"기준 유전자 {reference_gene}가 데이터에 없습니다.")
    if control_condition not in set(df["Condition"]):
        raise QPCRInputError(f"대조군 {control_condition}가 데이터에 없습니다.")
    targets = sorted(set(df["Gene"]) - {reference_gene})
    if not targets:
        raise QPCRInputError("표적 유전자가 없습니다.")

    # Check sample/gene completeness
    samples = df[["Sample", "Condition"]].drop_duplicates()
    expected = pd.MultiIndex.from_product(
        [samples["Sample"], [reference_gene, *targets]], names=["Sample", "Gene"]
    )
    observed = pd.MultiIndex.from_frame(df[["Sample", "Gene"]].drop_duplicates())
    absent = expected.difference(observed)
    if len(absent):
        preview = ", ".join(f"{s}/{g}" for s, g in absent[:8])
        raise QPCRInputError(f"시료별 기준/표적 유전자 측정이 누락되었습니다: {preview}")
    if not set(df.loc[df["Condition"] == control_condition, "Gene"]).issuperset(targets):
        raise QPCRInputError("대조군에 표적 유전자 측정이 없습니다.")

    # Process Technical Replicates (with Outlier QC)
    qc_rows = []
    outlier_warnings_count = 0

    grouped = df.groupby(["Sample", "Condition", "Gene"], sort=True)
    for (sample, condition, gene), group in grouped:
        ct_vals = group["Ct"].to_numpy()
        n_init = len(ct_vals)
        sd_init = np.std(ct_vals, ddof=1) if n_init > 1 else 0.0
        mean_init = np.mean(ct_vals)

        excluded_val = None
        used_vals = ct_vals.copy()
        status = "Pass"

        if n_init < min_technical_replicates:
            status = "Fail (Low Reps)"
        elif n_init > 1 and sd_init > max_technical_sd:
            if auto_exclude_outliers and n_init >= 3:
                # Exclude the furthest value from median
                median_val = np.median(ct_vals)
                distances = np.abs(ct_vals - median_val)
                outlier_idx = np.argmax(distances)
                excluded_val = ct_vals[outlier_idx]
                used_vals = np.delete(ct_vals, outlier_idx)
                
                sd_new = np.std(used_vals, ddof=1) if len(used_vals) > 1 else 0.0
                if len(used_vals) >= min_technical_replicates and sd_new <= max_technical_sd:
                    status = "Pass (Outlier Excluded)"
                    outlier_warnings_count += 1
                else:
                    status = "Fail (High SD after exclusion)"
            else:
                status = "Fail (High SD)"

        qc_rows.append({
            "Sample": sample,
            "Condition": condition,
            "Gene": gene,
            "Technical_N_Init": n_init,
            "SD_Init": sd_init,
            "Technical_N_Used": len(used_vals),
            "Mean_Ct": np.mean(used_vals),
            "Technical_SD": np.std(used_vals, ddof=1) if len(used_vals) > 1 else 0.0,
            "Excluded_Value": excluded_val,
            "Status": status
        })

    qc_df = pd.DataFrame(qc_rows)
    failed_qc = qc_df[qc_df["Status"].str.startswith("Fail")]
    if not failed_qc.empty:
        labels = ", ".join(f"{r.Sample}/{r.Gene}" for r in failed_qc.head(8).itertuples())
        raise QPCRInputError(f"기술 반복 QC 실패 ({labels}). 반복 수 또는 SD를 확인하세요. 필요시 허용 기준을 완화하세요.")

    # Reference Gene matching
    ref = (qc_df[qc_df["Gene"] == reference_gene][["Sample", "Condition", "Mean_Ct"]]
           .rename(columns={"Mean_Ct": "Reference_Ct"}))
    target = qc_df[qc_df["Gene"] != reference_gene].copy()
    sample_results = target.merge(ref, on=["Sample", "Condition"], validate="many_to_one")
    
    # Calculate Delta_Ct
    sample_results["Delta_Ct"] = sample_results["Mean_Ct"] - sample_results["Reference_Ct"]
    
    # Calibrator: Mean Delta_Ct of Control Condition
    calibrator = (sample_results[sample_results["Condition"] == control_condition]
                  .groupby("Gene")["Delta_Ct"].mean())
    sample_results["Control_Mean_Delta_Ct"] = sample_results["Gene"].map(calibrator)
    
    # Calculate Delta_Delta_Ct (Now controls have variation around 0, thus Fold_Change around 1)
    sample_results["Delta_Delta_Ct"] = sample_results["Delta_Ct"] - sample_results["Control_Mean_Delta_Ct"]
    sample_results["Fold_Change"] = np.exp2(-sample_results["Delta_Delta_Ct"])
    if not np.isfinite(sample_results["Fold_Change"]).all():
        raise QPCRInputError("상대 발현량 계산 결과가 유한하지 않습니다. Ct 값을 확인하세요.")

    # Summary
    summary = (sample_results.groupby(["Condition", "Gene"], sort=True)
               .agg(Biological_N=("Sample", "nunique"),
                    Mean_Delta_Ct=("Delta_Ct", "mean"),
                    SD_Delta_Ct=("Delta_Ct", "std"),
                    Geometric_Mean_FC=("Fold_Change", lambda x: float(np.exp(np.log(x).mean()))),
                    SD_Fold_Change=("Fold_Change", "std"))
               .reset_index())
    
    ref_summary = (ref.groupby("Condition", sort=True)["Reference_Ct"]
                   .agg(Biological_N="size", Mean_Reference_Ct="mean", SD_Reference_Ct="std")
                   .reset_index())

    # Statistical Hypothesis Testing on Delta_Ct (Welch's t-test or ANOVA)
    stat_rows = []
    unique_conditions = sorted(set(sample_results["Condition"]))
    for gene in targets:
        gene_data = sample_results[sample_results["Gene"] == gene]
        if len(unique_conditions) == 2:
            # Welch's t-test
            c1, c2 = unique_conditions
            g1 = gene_data[gene_data["Condition"] == c1]["Delta_Ct"].to_numpy()
            g2 = gene_data[gene_data["Condition"] == c2]["Delta_Ct"].to_numpy()
            if len(g1) >= 2 and len(g2) >= 2:
                t_stat, p_val = stats.ttest_ind(g1, g2, equal_var=False)
                stat_rows.append({
                    "Gene": gene,
                    "Test": f"Welch's t-test ({c1} vs {c2})",
                    "Statistic": t_stat,
                    "p-value": p_val,
                    "Significance": "*" if p_val < 0.05 else "ns"
                })
        elif len(unique_conditions) > 2:
            # One-Way ANOVA
            groups_data = [gene_data[gene_data["Condition"] == cond]["Delta_Ct"].to_numpy() for cond in unique_conditions]
            if all(len(g) >= 2 for g in groups_data):
                f_stat, p_val = stats.f_oneway(*groups_data)
                stat_rows.append({
                    "Gene": gene,
                    "Test": "One-way ANOVA",
                    "Statistic": f_stat,
                    "p-value": p_val,
                    "Significance": "*" if p_val < 0.05 else "ns"
                })
    
    stat_df = pd.DataFrame(stat_rows) if stat_rows else pd.DataFrame(columns=["Gene", "Test", "Statistic", "p-value", "Significance"])

    # Generate Warnings
    warnings = [
        "2^-ΔΔCt는 표적과 기준 유전자의 증폭 효율이 약 100%(E=2)로 비슷하다는 전제가 필요합니다.",
        "대조군 내부의 생물학적 편차(Biological variation)를 그래프 및 통계에 정상적으로 보존하여 반영하였습니다.",
    ]
    if outlier_warnings_count > 0:
        warnings.append(f"총 {outlier_warnings_count}개의 웰 측정군에서 편차가 큰 이상치(Outlier)를 1개씩 자동 제외하고 계산했습니다.")
    if (summary["Biological_N"] < 3).any():
        warnings.append("일부 조건·유전자의 생물학적 시료 수가 3 미만입니다. 통계 검정 및 변동성 해석에 주의하세요.")
    
    return AnalysisResult(qc_df, sample_results, summary, ref_summary, stat_df, warnings)
