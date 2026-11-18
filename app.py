"""Streamlit entry point: streamlit run app.py"""

import matplotlib.pyplot as plt
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
                factor_table[["돌연변이체 코드", "ABA 처리 코드", "스트레스"]] = list(factors)
                st.caption("시료 코드 해석: 첫 숫자=돌연변이체, 가운데=ABA 처리량 코드, X/H=스트레스 없음/있음. ABA 단위는 실험 기록으로 확인하세요.")
                st.dataframe(factor_table, hide_index=True)
            edited = st.data_editor(mapping, hide_index=True, disabled=["기기 시료명"], key="sample_mapping")
            if edited["실험 조건"].isna().any() or edited["실험 조건"].astype(str).str.strip().eq("").any():
                st.warning("모든 기기 시료명에 실험 조건을 입력해야 분석할 수 있습니다.")
                st.stop()
            if edited["생물학적 시료 ID"].isna().any() or edited["생물학적 시료 ID"].astype(str).str.strip().eq("").any():
                raise QPCRInputError("생물학적 시료 ID에 빈 값이 있습니다.")
            data = data.merge(edited, left_on="Sample", right_on="기기 시료명", validate="many_to_one")
            data["Sample"] = data["생물학적 시료 ID"].astype(str).str.strip()
            data["Condition"] = data["실험 조건"].astype(str).str.strip()
            data["Replicate"] = data.groupby(["Sample", "Condition", "Gene"]).cumcount() + 1
            data = data[["Sample", "Condition", "Gene", "Replicate", "Ct"]]
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
        min_reps = st.number_input("시료·유전자별 최소 기술 반복 수", min_value=2, value=2, step=1)
        st.dataframe(data.head(20), hide_index=True)

        if st.button("분석 실행", type="primary"):
            result = analyze(data, reference, control, max_sd, int(min_reps))
            st.success("입력 및 기술 반복 QC 통과")
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

            fig, ax = plt.subplots(figsize=(9, 5))
            groups = result.sample_results.groupby(["Gene", "Condition"], sort=True)
            labels = []
            for i, ((gene, condition), group) in enumerate(groups):
                ax.scatter([i] * len(group), group["Fold_Change"], alpha=0.7)
                geo_mean = result.summary.loc[
                    (result.summary["Gene"] == gene) & (result.summary["Condition"] == condition),
                    "Geometric_Mean_FC",
                ].iloc[0]
                ax.scatter(i, geo_mean, marker="_", s=450, linewidths=3, color="black")
                labels.append(f"{gene}\n{condition}")
            ax.set_xticks(range(len(labels)), labels)
            ax.set_ylabel("Fold change (2^-ΔΔCt)")
            ax.axhline(1, linestyle="--", color="gray")
            ax.set_yscale("log", base=2)
            ax.set_title("Dots: biological samples; black bars: geometric mean")
            fig.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

            csv = result.sample_results.to_csv(index=False).encode("utf-8-sig")
            st.download_button("시료별 결과 CSV 다운로드", csv, "qpcr_sample_results.csv", "text/csv")
            summary_csv = result.summary.to_csv(index=False).encode("utf-8-sig")
            st.download_button("조건별 요약 CSV 다운로드", summary_csv, "qpcr_summary.csv", "text/csv")
    except QPCRInputError as exc:
        st.error(str(exc))
