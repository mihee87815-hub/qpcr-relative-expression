"""Streamlit entry point with CSV template and advanced analysis tools."""

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from qpcr_analysis import QPCRInputError, analyze, read_qpcr_csv


st.set_page_config(page_title="qPCR 상대 발현량 분석", layout="wide")
st.title("qPCR 상대 발현량 분석")
st.write("시료별 기술 반복 Ct 이상치 정제 및 ΔCt, ΔΔCt, 2^-ΔΔCt 계산과 통계 검정을 자동 수행합니다.")

# Template download helper
template_df = pd.DataFrame([
    {"Sample": "C1", "Condition": "Control", "Gene": "GAPDH", "Replicate": 1, "Ct": 21.1},
    {"Sample": "C1", "Condition": "Control", "Gene": "GAPDH", "Replicate": 2, "Ct": 21.2},
    {"Sample": "C1", "Condition": "Control", "Gene": "GAPDH", "Replicate": 3, "Ct": 21.0},
    {"Sample": "C1", "Condition": "Control", "Gene": "MYC", "Replicate": 1, "Ct": 25.4},
    {"Sample": "C1", "Condition": "Control", "Gene": "MYC", "Replicate": 2, "Ct": 25.5},
    {"Sample": "C1", "Condition": "Control", "Gene": "MYC", "Replicate": 3, "Ct": 26.9}, # 26.9 is intentional outlier
], columns=["Sample", "Condition", "Gene", "Replicate", "Ct"])

csv_template = template_df.to_csv(index=False).encode("utf-8-sig")
st.download_button(
    label="📄 예제/입력 템플릿 CSV 다운로드",
    data=csv_template,
    file_name="qpcr_input_template.csv",
    mime="text/csv"
)

uploaded = st.file_uploader("qPCR 결과 CSV 파일 업로드", type="csv")
if uploaded is not None:
    try:
        data = read_qpcr_csv(uploaded.getvalue())
        if not {"Sample", "Condition", "Gene", "Replicate", "Ct"}.issubset(data.columns):
            raise QPCRInputError("필수 열이 없습니다: Sample, Condition, Gene, Replicate, Ct")
        genes = sorted(data["Gene"].dropna().astype(str).unique())
        conditions = sorted(data["Condition"].dropna().astype(str).unique())
        if not genes or not conditions:
            raise QPCRInputError("Gene 또는 Condition 값이 없습니다.")
        
        col1, col2 = st.columns(2)
        with col1:
            reference = st.selectbox("기준 유전자 (Reference Gene)", genes)
            max_sd = st.number_input("기술 반복 Ct SD 허용 상한 (초과 시 이상치 1개 배제 시도)", min_value=0.01, value=0.30, step=0.05)
        with col2:
            default_control = next((i for i, name in enumerate(conditions) if name.lower() in {"mock", "control", "untreated"}), 0)
            control = st.selectbox("대조군 조건 (Control Condition)", conditions, index=default_control)
            min_reps = st.number_input("시료·유전자별 최소 기술 반복 수 (QC 통과용)", min_value=2, value=2, step=1)
        
        auto_exclude = st.checkbox("기술 반복 이상치 자동 플래그 및 배제 기능 사용 (3반복 이상에 적용)", value=True)

        st.subheader("업로드된 원본 데이터 미리보기")
        st.dataframe(data.head(15), hide_index=True)

        if st.button("분석 및 통계 검정 실행", type="primary"):
            result = analyze(data, reference, control, max_sd, int(min_reps), auto_exclude)
            st.success("🎉 입력 데이터 및 기술 반복 QC 통과!")
            
            for warning in result.warnings:
                st.warning(warning)
            
            if not result.statistical_tests.empty:
                st.subheader("🔬 통계적 유의성 검정 결과 (ΔCt 기준)")
                st.dataframe(result.statistical_tests, hide_index=True)
            
            st.subheader("📈 조건별 상대 발현량 요약 (Fold Change 기하평균)")
            st.dataframe(result.summary, hide_index=True)
            
            st.subheader("🧬 생물학적 시료별 상세 결과")
            st.dataframe(result.sample_results, hide_index=True)
            
            st.subheader("🧪 기술 반복 QC 상태")
            st.dataframe(result.technical_qc, hide_index=True)

            # Relative Expression Plot
            fig, ax = plt.subplots(figsize=(10, 6))
            groups = result.sample_results.groupby(["Gene", "Condition"], sort=True)
            labels = []
            for i, ((gene, condition), group) in enumerate(groups):
                # Plot individual samples
                ax.scatter([i] * len(group), group["Fold_Change"], alpha=0.6, label="Samples" if i==0 else "", zorder=3, edgecolors="k")
                geo_mean = result.summary.loc[
                    (result.summary["Gene"] == gene) & (result.summary["Condition"] == condition),
                    "Geometric_Mean_FC",
                ].iloc[0]
                # Plot Geometric Mean as a solid line
                ax.scatter(i, geo_mean, marker="_", s=500, linewidths=3, color="black", zorder=4)
                labels.append(f"{gene}\n({condition})")
            
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45)
            ax.set_ylabel("Fold Change (2^-ΔΔCt)")
            ax.axhline(1, linestyle="--", color="gray", alpha=0.7)
            ax.set_yscale("log", base=2)
            ax.set_title("Relative Fold Change (Individual Samples & Geometric Mean)")
            ax.grid(axis="y", linestyle=":", alpha=0.5)
            fig.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

            # Download Buttons
            st.download_button("시료별 상세 결과 CSV 다운로드", result.sample_results.to_csv(index=False).encode("utf-8-sig"), "qpcr_sample_results.csv", "text/csv")
            st.download_button("조건별 요약 결과 CSV 다운로드", result.summary.to_csv(index=False).encode("utf-8-sig"), "qpcr_summary.csv", "text/csv")
    except QPCRInputError as exc:
        st.error(str(exc))
