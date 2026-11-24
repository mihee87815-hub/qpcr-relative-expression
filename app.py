"""Streamlit entry point: streamlit run app.py"""

from io import BytesIO
import hashlib

import matplotlib.pyplot as plt
import streamlit as st

from qpcr_analysis import QPCRInputError, analyze, parse_sample_code, read_qpcr_upload
from qpcr_plot import FACTOR_COLUMNS, PLOT_TYPES, make_figure, ordered_values, prepare_plot_data
from qpcr_report import write_report
from qpcr_merge import read_files, merge_samples


st.set_page_config(page_title="qPCR 상대 발현량 분석", layout="wide")
st.title("qPCR 상대 발현량 분석")
st.write("시료별 기술 반복 Ct를 평균낸 뒤 ΔCt, ΔΔCt, 2^-ΔΔCt를 계산합니다.")
st.caption("앱 형식 CSV 또는 QuantStudio Results 시트가 있는 XLS/XLSX 파일을 올리세요.")

uploads = st.file_uploader("qPCR 결과 파일 (최대 3개)", type=["csv", "xls", "xlsx"], accept_multiple_files=True)
if uploads:
    try:
        sample_codes = {}
        if len(uploads) > 3:
            raise QPCRInputError("파일은 최대 3개까지 올릴 수 있습니다. 추가 파일을 제거하세요.")
        if len(uploads) > 1:
            combined, mapping, audit = read_files([(f.name, f.getvalue()) for f in uploads])
            st.subheader("파일 통합 설정")
            st.dataframe(audit, hide_index=True)
            st.info("기본값은 파일마다 독립 생물학적 시료로 구분합니다. 같은 실제 시료를 재측정한 파일이면 생물학적 시료 ID를 동일하게 바꾸세요. 독립 시료라도 처리 조건이 같으면 실험 조건을 같게 지정하세요.")
            mapping_key = hashlib.sha256(b''.join(f.getvalue() for f in uploads)).hexdigest()[:12]
            edited_multi = st.data_editor(mapping, hide_index=True, disabled=["입력 ID", "파일", "원본 시료명"], key=f"merge_{mapping_key}")
            st.warning("서로 다른 실험일·플레이트의 Ct를 통합하려면 측정 조건과 기준 유전자 안정성, 공통 보정 시료를 확인해야 합니다. 이 앱은 파일 간 배치 효과를 보정하지 않습니다.")
            confirmed = st.checkbox("파일 간 시료 관계와 통합 가능한 측정 조건을 확인했습니다", key=f"confirm_{mapping_key}")
            if not confirmed:
                st.stop()
            data, sample_codes = merge_samples(combined, edited_multi)
            instrument_export = False
            st.download_button("통합 시료 설정 CSV 다운로드", edited_multi.to_csv(index=False).encode("utf-8-sig"), "qpcr_file_mapping.csv", "text/csv")
        else:
            uploaded = uploads[0]
            data, instrument_export = read_qpcr_upload(uploaded.getvalue(), uploaded.name)
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
        elif not sample_codes:
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

        analysis_key = hashlib.sha256(
            (data.to_csv(index=False) + repr((reference, control, max_sd, int(min_reps), allow_high_sd))).encode("utf-8")
        ).hexdigest()
        if st.button("분석 실행", type="primary"):
            st.session_state["qpcr_result"] = analyze(data, reference, control, max_sd, int(min_reps), allow_high_sd)
            st.session_state["qpcr_result_key"] = analysis_key
        result = st.session_state.get("qpcr_result") if st.session_state.get("qpcr_result_key") == analysis_key else None
        if result is not None:
            st.info("상대 발현량 계산 완료. 기술 반복 QC 표와 경고를 함께 확인하세요. 이 결과는 탐색적 계산입니다.")
            for warning in result.warnings:
                st.warning(warning)
            st.subheader("그래프 그리기 도우미")
            plotted = prepare_plot_data(result.sample_results, sample_codes)
            coded = all(column in plotted for column in FACTOR_COLUMNS)
            st.caption("시료 점은 독립 생물학적 시료입니다. 막대와 선은 기하평균이며 n은 생물학적 시료 수입니다. n=1에는 오차막대를 그리지 않습니다.")
            if coded:
                st.caption("시료 코드: X=상온, H=고온, 0=Col-1, 1=ost1. ABA 값은 실험 기록의 단위가 없어 코드로 표시합니다.")
            elif "염 수준" not in plotted:
                st.info("처리 코드를 인식하지 못해 실험 조건별 그래프를 제공합니다.")
            graph_key = analysis_key[:12]
            gene = st.selectbox("그래프 유전자", sorted(plotted["Gene"].unique()), key=f"gene_{graph_key}")
            if "염 수준" in plotted:
                st.caption("앞 숫자는 ABA 농도, 뒤 0/LOW/HIGH는 염 처리 수준입니다. 농도 단위와 LOW/HIGH의 실제 값은 실험 기록을 확인하세요.")
                grouping = st.radio("비교 방식", ["ABA 농도별로 염 처리 비교", "염 수준별로 ABA 농도 비교"], key=f"salt_mode_{graph_key}")
                filtered = plotted.loc[plotted["Gene"] == gene].copy()
                fixed = "ABA 농도" if grouping == "ABA 농도별로 염 처리 비교" else "염 수준"
                x_column = "염 수준" if fixed == "ABA 농도" else "ABA 농도"
                chosen = st.selectbox(f"비교할 {fixed}", ordered_values(filtered[fixed], fixed), key=f"salt_fixed_{fixed}_{graph_key}")
                filtered = filtered.loc[filtered[fixed] == chosen]
                chart_options = PLOT_TYPES if x_column == "ABA 농도" else PLOT_TYPES[:2]
                chart_type = st.selectbox("그래프 형식", chart_options, key=f"salt_chart_{fixed}_{graph_key}")
                series_column = None
                facet_column = fixed
                st.caption(f"현재 표시: {fixed} = {chosen}. 상대 발현량의 보정 대조군은 {control}로 유지됩니다.")
            else:
                chart_options = PLOT_TYPES if coded else PLOT_TYPES[:2]
                chart_type = st.selectbox("그래프 형식", chart_options, key=f"chart_{graph_key}")
                dimensions = {"실험 조건": "Condition"}
                if coded:
                    dimensions = {"ABA 농도 코드": "ABA 코드", "온도": "온도", "실험군": "실험군", **dimensions}
                x_options = ["ABA 농도 코드"] if chart_type == "기하평균 선" else list(dimensions)
                x_label = st.selectbox("가로축", x_options, key=f"x_{chart_type}_{graph_key}")
                x_column = dimensions[x_label]
                group_options = {label: column for label, column in dimensions.items() if column != x_column}
                series_label = st.selectbox("색상으로 구분", ["없음", *group_options], key=f"series_{x_label}_{graph_key}")
                series_column = group_options.get(series_label)
                facet_options = {label: column for label, column in group_options.items() if column != series_column}
                facet_label = st.selectbox("패널 나누기", ["없음", *facet_options], key=f"facet_{x_label}_{series_label}_{graph_key}")
                facet_column = facet_options.get(facet_label)
                filtered = plotted.loc[plotted["Gene"] == gene].copy()
                for label, column in dimensions.items():
                    values = ordered_values(filtered[column], column)
                    selected = st.multiselect(f"{label} 선택", values, default=values, key=f"filter_{column}_{graph_key}")
                    filtered = filtered.loc[filtered[column].isin(selected)]
            log_scale = st.checkbox("Y축 로그2", value=True, key=f"log_{graph_key}")
            if filtered.empty:
                st.warning("선택한 조건에 표시할 시료가 없습니다.")
            else:
                fig = make_figure(filtered, gene, chart_type, x_column, series_column, facet_column, log_scale)
                st.pyplot(fig)
                image = BytesIO()
                fig.savefig(image, format="png", dpi=200, bbox_inches="tight")
                st.download_button("그래프 PNG 다운로드", image.getvalue(), f"qpcr_{gene}_plot.png", "image/png")
                plt.close(fig)

            with st.expander("상세 결과표 및 QC"):
                st.subheader("조건별 요약")
                st.dataframe(result.summary, hide_index=True)
                st.subheader("기준 유전자 Ct")
                st.dataframe(result.reference_summary, hide_index=True)
                st.subheader("생물학적 시료별 결과")
                st.dataframe(result.sample_results, hide_index=True)
                st.subheader("기술 반복 QC")
                st.dataframe(result.technical_qc, hide_index=True)

            st.subheader("결과 보고서 작성 도우미")
            st.caption("관찰값과 계산 근거를 연결해 설명합니다. 통계적 유의성이나 생물학적 기작은 추정하지 않습니다.")
            report_title = st.text_input("보고서 제목", value="qPCR 상대 발현량 결과 보고서")
            experiment_note = st.text_area("실험 정보·질문 (선택)", help="ABA 농도와 단위, 처리 시간, 시료 채취 시점 등을 기록하세요.")
            report_note = experiment_note
            if len(uploads) > 1:
                report_note += "\n통합 파일: " + ", ".join(f.name for f in uploads) + "\n파일 간 배치 효과는 보정하지 않았습니다. 생물학적 반복 수는 사용자가 지정한 시료 ID를 기준으로 계산했습니다."
            report = write_report(result, reference, control, report_title, report_note, sample_codes)
            with st.expander("보고서 미리보기", expanded=True):
                st.markdown(report)
            st.download_button("보고서 Markdown 다운로드", report.encode("utf-8-sig"),
                               "qpcr_report.md", "text/markdown")

            csv = result.sample_results.to_csv(index=False).encode("utf-8-sig")
            st.download_button("시료별 결과 CSV 다운로드", csv, "qpcr_sample_results.csv", "text/csv")
            summary_csv = result.summary.to_csv(index=False).encode("utf-8-sig")
            st.download_button("조건별 요약 CSV 다운로드", summary_csv, "qpcr_summary.csv", "text/csv")
    except QPCRInputError as exc:
        st.error(str(exc))
