"""Deterministic relative-expression analysis for technical-replicate qPCR CSVs."""

from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd


REQUIRED = {"Sample", "Condition", "Gene", "Replicate", "Ct"}


def parse_sample_code(value: str) -> tuple[str, str, str] | None:
    """Parse user-defined genotype/ABA/stress codes without assigning units."""
    match = re.fullmatch(r"([01])(0|5|10)([XH])", str(value).strip())
    if not match:
        return None
    genotype, aba, stress = match.groups()
    return genotype, aba, "스트레스 없음" if stress == "X" else "스트레스 있음"


class QPCRInputError(ValueError):
    pass


@dataclass
class AnalysisResult:
    technical_qc: pd.DataFrame
    sample_results: pd.DataFrame
    summary: pd.DataFrame
    reference_summary: pd.DataFrame
    warnings: list[str]


def read_qpcr_csv(content: bytes) -> pd.DataFrame:
    try:
        return pd.read_csv(BytesIO(content), encoding="utf-8-sig")
    except (UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise QPCRInputError(f"CSV 파일을 읽을 수 없습니다: {exc}") from exc


def read_qpcr_upload(content: bytes, filename: str) -> tuple[pd.DataFrame, bool]:
    """Read either the app template or a QuantStudio Results export.

    Instrument exports need an explicit sample-to-condition mapping in the UI.
    """
    suffix = Path(filename).suffix.lower()
    if suffix == ".csv":
        try:
            text = content.decode("utf-8-sig")
        except UnicodeDecodeError as exc:
            raise QPCRInputError("CSV 인코딩을 읽을 수 없습니다. UTF-8 CSV로 다시 내보내세요.") from exc
        lines = text.splitlines()
        header = next((i for i, line in enumerate(lines) if line.startswith("Well,Well Position,")), None)
        if header is not None:
            if "CT" not in lines[header].split(","):
                raise QPCRInputError("이 CSV는 웰 배치(Sample Setup)만 포함하고 Ct가 없습니다. 같은 실험의 Results 시트가 있는 .xls/.xlsx 파일을 올리세요.")
            table = pd.read_csv(BytesIO(content), skiprows=header)
            return _read_instrument_results(table), True
        return read_qpcr_csv(content), False
    if suffix in {".xls", ".xlsx"}:
        try:
            workbook = pd.ExcelFile(BytesIO(content))
            if "Results" not in workbook.sheet_names:
                raise QPCRInputError("Excel 파일에 Results 시트가 없습니다. Ct 결과가 포함된 실험 결과 파일을 올리세요.")
            raw = pd.read_excel(workbook, sheet_name="Results", header=None)
        except QPCRInputError:
            raise
        except (ValueError, ImportError, OSError) as exc:
            raise QPCRInputError(f"Excel 결과 파일을 읽을 수 없습니다: {exc}") from exc
        header = next((i for i, row in raw.iterrows()
                       if str(row.iloc[0]).strip() == "Well"
                       and "CT" in set(row.astype(str).str.strip())), None)
        if header is None:
            raise QPCRInputError("Results 시트에서 Well과 CT가 있는 헤더를 찾지 못했습니다.")
        table = raw.iloc[header + 1:].copy()
        table.columns = raw.iloc[header].astype(str).str.strip()
        return _read_instrument_results(table), True
    raise QPCRInputError("CSV, XLS, XLSX 파일만 올릴 수 있습니다.")


def _read_instrument_results(table: pd.DataFrame) -> pd.DataFrame:
    required = {"Well Position", "Sample Name", "Target Name", "CT"}
    if not required.issubset(table.columns):
        raise QPCRInputError("결과 표에 Well Position, Sample Name, Target Name, CT 열이 필요합니다.")
    rows = table.dropna(subset=["Well Position", "Sample Name", "Target Name"]).copy()
    if rows.empty:
        raise QPCRInputError("분석할 웰이 없습니다.")
    ntc = (rows["Task"].astype(str).str.upper().eq("NTC")
           if "Task" in rows else pd.Series(False, index=rows.index))
    ntc_rows = rows[ntc]
    if not ntc_rows.empty:
        ntc_ct = pd.to_numeric(ntc_rows["CT"], errors="coerce")
        if ntc_ct.notna().any():
            wells = ", ".join(ntc_rows.loc[ntc_ct.notna(), "Well Position"].astype(str).head(8))
            raise QPCRInputError(f"NTC에서 Ct가 검출됐습니다 ({wells}). 오염·비특이 증폭 가능성을 확인하세요.")
        rows = rows[~ntc].copy()
    if rows.empty:
        raise QPCRInputError("NTC를 제외하면 분석할 시료 웰이 없습니다.")
    if "Omit" in rows and rows["Omit"].astype(str).str.lower().isin(["true", "1", "yes"]).any():
        raise QPCRInputError("Omit 표시된 웰이 있습니다. 제외 이유를 확인하고 결과 파일을 정리하세요.")
    invalid_ct = pd.to_numeric(rows["CT"], errors="coerce").isna()
    invalid_amp = (rows["Amp Status"].astype(str).str.strip().ne("Amp")
                   if "Amp Status" in rows else pd.Series(False, index=rows.index))
    invalid = rows[invalid_ct | invalid_amp]
    if not invalid.empty:
        examples = ", ".join(
            f"{r['Well Position']} ({r['Sample Name']}/{r['Target Name']}: {r['CT']})"
            for _, r in invalid.head(8).iterrows()
        )
        raise QPCRInputError(
            f"Ct를 계산할 수 없는 웰이 있습니다: {examples}. 원본 증폭 결과를 확인하고 재측정하거나, "
            "사전에 정한 제외 기준에 따라 처리하세요. 남은 기술 반복 수가 최소 기준 미만이면 분석할 수 없습니다."
        )
    rows["Sample"] = rows["Sample Name"].astype(str).str.strip()
    rows["Gene"] = rows["Target Name"].astype(str).str.strip()
    rows["Replicate"] = rows.groupby(["Sample", "Gene"]).cumcount() + 1
    rows["Condition"] = ""  # Never infer experimental treatment from a sample code.
    result = rows[["Sample", "Condition", "Gene", "Replicate", "CT"]].rename(columns={"CT": "Ct"}).reset_index(drop=True)
    result.attrs["ntc_wells"] = ntc_rows["Well Position"].astype(str).tolist()
    return result


def analyze(
    data: pd.DataFrame,
    reference_gene: str,
    control_condition: str,
    max_technical_sd: float = 0.3,
    min_technical_replicates: int = 2,
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

    # Each biological sample needs a reference and every target. Missing rows must
    # never disappear silently during the later join.
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

    qc = (df.groupby(["Sample", "Condition", "Gene"], sort=True)["Ct"]
          .agg(Technical_N="size", Mean_Ct="mean", Technical_SD="std")
          .reset_index())
    bad = qc[(qc["Technical_N"] < min_technical_replicates)
             | (qc["Technical_SD"] > max_technical_sd)]
    if not bad.empty:
        labels = ", ".join(f"{r.Sample}/{r.Gene}" for r in bad.head(8).itertuples())
        raise QPCRInputError(f"기술 반복 QC 실패 ({labels}). 반복 수 또는 SD를 확인하세요.")

    ref = (qc[qc["Gene"] == reference_gene][["Sample", "Condition", "Mean_Ct"]]
           .rename(columns={"Mean_Ct": "Reference_Ct"}))
    target = qc[qc["Gene"] != reference_gene].copy()
    sample = target.merge(ref, on=["Sample", "Condition"], validate="many_to_one")
    sample["Delta_Ct"] = sample["Mean_Ct"] - sample["Reference_Ct"]
    calibrator = (sample[sample["Condition"] == control_condition]
                  .groupby("Gene")["Delta_Ct"].mean())
    sample["Control_Mean_Delta_Ct"] = sample["Gene"].map(calibrator)
    sample["Delta_Delta_Ct"] = sample["Delta_Ct"] - sample["Control_Mean_Delta_Ct"]
    sample["Fold_Change"] = np.exp2(-sample["Delta_Delta_Ct"])
    if not np.isfinite(sample["Fold_Change"]).all():
        raise QPCRInputError("상대 발현량 계산 결과가 유한하지 않습니다. Ct 값을 확인하세요.")

    summary = (sample.groupby(["Condition", "Gene"], sort=True)
               .agg(Biological_N=("Sample", "nunique"),
                    Mean_Delta_Ct=("Delta_Ct", "mean"),
                    SD_Delta_Ct=("Delta_Ct", "std"),
                    Geometric_Mean_FC=("Fold_Change", lambda x: float(np.exp(np.log(x).mean()))))
               .reset_index())
    ref_summary = (ref.groupby("Condition", sort=True)["Reference_Ct"]
                   .agg(Biological_N="size", Mean_Reference_Ct="mean", SD_Reference_Ct="std")
                   .reset_index())
    warnings = [
        "2^-ΔΔCt는 표적과 기준 유전자의 증폭 효율이 적절하고 비슷하다는 전제가 필요합니다.",
        "기준 유전자의 조건 간 안정성은 이 앱이 검증하지 않습니다. 별도 실험 근거로 확인하세요.",
        "이 결과는 기술 반복을 시료별로 평균낸 기술 통계입니다. 유의성 검정이나 생물학적 결론을 제공하지 않습니다.",
    ]
    if (summary["Biological_N"] < 3).any():
        warnings.append("일부 조건·유전자의 생물학적 시료 수가 3 미만입니다. 변동성 해석에 주의하세요.")
    return AnalysisResult(qc, sample, summary, ref_summary, warnings)
