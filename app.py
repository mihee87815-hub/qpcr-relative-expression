"""Streamlit entry point: streamlit run app.py"""

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from qpcr_analysis import QPCRInputError, analyze, parse_sample_code, read_qpcr_upload


st.set_page_config(page_title="qPCR 상대 발현량 분석", layout="wide")
st.title("qPCR 상대 발현량 분석")
st.write("시료별 기술 반복 Ct를 평균낸 뒤 ΔCt, ΔΔCt, 2^-ΔΔCt를 계산합니다.")
st.caption("앱 형식 CSV 또는 QuantStudio Results 시트가 있는 XLS/XLSX 파일을 올리세요.")

uploaded = st.file_uploader("qPCR 결과 파일", type=["csv", "xls", "xlsx"])
if uploaded is not None:
    try:
        data, instrument_export = read_qpcr_upload(uploaded.getvalue(), uploaded.name)
        sample_codes = {}
        excluded_wells = data.attrs.get("excluded_wells", [])
        excluded_ct_rows = data.attrs.get("excluded_ct_rows", 0)
        if excluded_wells:
            st.warning(f"Undetermined 웰 {len(excluded_wells)}개를 계산에서 제외했습니다. 제외 내역을 아래에서 확인하세요.")
            st.dataframe(excluded_wells, hide_index=True)
        if excluded_ct_rows:
            st.warning(f"Undetermined Ct {excluded_ct_rows}행을 계산에서 제외했습니다.")
        if instrument_export:
            ntc_wells = data.attrs.get("ntc_wells", [])
            if ntc_wells:
                st.success(f"NTC {len(ntc_wells)}개 웰은 Ct가 검출되지 않아 상대 발현량 계산에서 제외했습니다: {', '.join(ntc_wells)}")
            st.info("기기 결과의 시료명에는 처리 조건 정보가 없습니다. 아래 표에서 각 시료의 조건을 반드시 지정하세요. 같은 생물학적 시료에서 나온 웰이면 생물학적 시료 ID도 동일하게 설정하세요.")
            mapping = data[["Sample"]].drop_duplicates().sort_values("Sample").reset_index(drop=True)
            mapping.columns = ["기기 시료명"]
            mapping["실험 조건"] = mapping["기기 시료명"]
            mapping["생물학적 시료 ID"] = mapping["기기 시료명"]
            st.caption("기본값은 시료명마다 별도 조건·생물학적 시료로 둡니다. ABA 농도가 같은 독립 시료들은 실험 조건을 같은 이름으로 바꾸세요. 농도 코드의 실제 값과 단위는 실험 기록을 보고 입력하세요.")
            factors = mapping["기기 시료명"].map(parse_sample_code)
            if factors.notna().all():
                factor_table = mapping[["기기 시료명"]].copy()
                factor_table[["실험군", "ABA 처리 코드", "온도"]] = [
                    ({"0": "Col-1", "1": "ost1"}[item[0]], item[1],
                     {"X": "상온", "H": "고온"}[str(name)[-1]])
                    for name, item in zip(mapping["기기 시료명"], factors)
                ]
                st.caption("시료 코드 해석: 0=Col-1, 1=ost1, X=상온, H=고온. ABA 단위는 실험 기록으로 확인하세요.")
                st.dataframe(factor_table, hide_index=True)
            edited = st.data_editor(mapping, hide_index=True, disabled=["기기 시료명"], key="sample_mapping")
            if edited["실험 조건"].isna().any() or edited["실험 조건"].astype(str).str.strip().eq("").any():
                st.warning("모든 기기 시료명에 실험 조건을 입력해야 분석할 수 있습니다.")
                st.stop()
            if edited["생물학적 시료 ID"].isna().any() or edited["생물학적 시료 ID"].astype(str).str.strip().eq("").any():
                raise QPCRInputError("생물학적 시료 ID에 빈 값이 있습니다.")
            sample_codes = dict(zip(edited["생물학적 시료 ID"].astype(str).str.strip(), edited["기기 시료명"].astype(str)))
            data = data.merge(edited, left_on="Sample", right_on="기기 시료명", validate="many_to_one")
            data["Sample"] = data["생물학적 시료 ID"].astype(str).str.strip()
            data["Condition"] = data["실험 조건"].astype(str).str.strip()
            data["Replicate"] = data.groupby(["Sample", "Condition", "Gene"]).cumcount() + 1
            data = data[["Sample", "Condition", "Gene", "Replicate", "Ct"]]
        else:
            sample_codes = {str(sample): str(sample) for sample in data["Sample"].unique()}
        if not {"Sample", "Condition", "Gene", "Replicate", "Ct"}.issubset(data.columns):
            raise QPCRInputError("필수 열이 없습니다: Sample, Condition, Gene, Replicate, Ct")
        genes = sorted(data["Gene"].dropna().astype(str).unique())
        conditions = sorted(data["Condition"].dropna().astype(str).unique())
        if not genes or not conditions:
            raise QPCRInputError("Gene 또는 Condition 값이 없습니다.")
        default_reference = genes.index("PP2A") if "PP2A" in genes else 0
        reference = st.selectbox("기준 유전자", genes, index=default_reference)
        default_control = next((i for i, name in enumerate(conditions) if name.lower() in {"00x", "mock", "control", "untreated"}), 0)
        control = st.selectbox("대조군 조건", conditions, index=default_control)
        st.info(f"현재 대조군: {control} — 그래프에서 이 조건의 기하평균이 1입니다. 분석 전에 확인하세요.")
        max_sd = st.number_input("기술 반복 Ct SD 허용 상한", min_value=0.01, value=0.30, step=0.05)
        min_reps = st.number_input("시료·유전자별 최소 기술 반복 수", min_value=1, value=1, step=1)
        allow_high_sd = st.checkbox("SD 기준을 넘는 묶음도 포함해 탐색적으로 분석", value=instrument_export)
        st.caption("기술 반복 1회인 묶음은 SD를 계산할 수 없습니다. SD 기준 초과 묶음은 결과에 표시되며 원본 증폭곡선을 확인해야 합니다.")
        st.dataframe(data.head(20), hide_index=True)

        if st.button("분석 실행", type="primary"):
            result = analyze(data, reference, control, max_sd, int(min_reps), allow_high_sd)
            st.info("상대 발현량 계산 완료. 기술 반복 QC 표와 경고를 함께 확인하세요. 이 결과는 탐색적 계산입니다.")
            for warning in result.warnings:
                st.warning(warning)
            st.subheader("조건별 요약")
            st.dataframe(result.summary, hide_index=True)
            st.subheader("기준 유전자 Ct")
            st.dataframe(result.reference_summary, hide_index=True)
            st.subheader("생물학적 시료별 결과")
            st.dataframe(result.sample_results, hide_index=True)
            st.subheader("기술 반복 QC")
            st.dataframe(result.technical_qc, hide_index=True)

            plotted = result.sample_results.copy()
            plotted["시료 코드"] = plotted["Sample"].map(sample_codes)
            factors = plotted["시료 코드"].map(parse_sample_code)
            if factors.notna().all():
                plotted["실험군"] = factors.map(lambda item: {"0": "Col-1", "1": "ost1"}[item[0]])
                plotted["ABA 코드"] = factors.map(lambda item: int(item[1]))
                plotted["온도"] = plotted["시료 코드"].str[-1].map({"X": "상온", "H": "고온"})
                st.caption("시료 코드 기준: X=상온, H=고온, 0=Col-1, 1=ost1. ABA 값은 파일에 단위가 없어 처리 코드로 표시합니다. 점=생물학적 시료, 검은 가로선=기하평균.")
                for gene, gene_data in plotted.groupby("Gene", sort=True):
                    st.subheader(f"{gene} — 온도·실험군·ABA별 발현량")
                    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
                    doses = sorted(gene_data["ABA 코드"].unique())
                    for row, temperature in enumerate(("상온", "고온")):
                        for col, genotype in enumerate(("Col-1", "ost1")):
                            ax = axes[row, col]
                            panel = gene_data[(gene_data["온도"] == temperature) & (gene_data["실험군"] == genotype)]
                            for x, dose in enumerate(doses):
                                values = panel.loc[panel["ABA 코드"] == dose, "Fold_Change"].to_numpy()
                                if len(values):
                                    offsets = np.linspace(-0.08, 0.08, len(values)) if len(values) > 1 else [0]
                                    ax.scatter(np.full(len(values), x) + offsets, values, alpha=0.75)
                                    ax.scatter(x, np.exp(np.log(values).mean()), marker="_", s=450, linewidths=3, color="black")
                            temperature_label = {"상온": "Ambient", "고온": "Heat"}[temperature]
                            ax.set_title(f"{temperature_label} · {genotype}")
                            ax.set_xticks(range(len(doses)), [str(dose) for dose in doses])
                            ax.axhline(1, linestyle="--", color="gray", linewidth=1)
                            ax.set_yscale("log", base=2)
                            ax.set_xlabel("ABA dose code")
                            if col == 0:
                                ax.set_ylabel("Fold change (2^-DeltaDeltaCt)")
                    fig.tight_layout()
                    st.pyplot(fig)
                    plt.close(fig)
            else:
                st.warning("시료 코드가 00X/05H 형식이 아니어서 온도·실험군·ABA별 그래프를 그릴 수 없습니다. 결과 표와 CSV는 확인할 수 있습니다.")

            csv = result.sample_results.to_csv(index=False).encode("utf-8-sig")
            st.download_button("시료별 결과 CSV 다운로드", csv, "qpcr_sample_results.csv", "text/csv")
            summary_csv = result.summary.to_csv(index=False).encode("utf-8-sig")
            st.download_button("조건별 요약 CSV 다운로드", summary_csv, "qpcr_summary.csv", "text/csv")
    except QPCRInputError as exc:
        st.error(str(exc))
