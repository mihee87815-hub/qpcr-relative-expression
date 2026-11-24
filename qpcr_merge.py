"""Read up to three uploads while preserving file and sample provenance."""
import hashlib
import pandas as pd
from qpcr_analysis import QPCRInputError, read_qpcr_upload


def read_files(files):
    if not 1 <= len(files) <= 3:
        raise QPCRInputError("파일은 최대 3개까지 올릴 수 있습니다.")
    frames, mapping, audit, seen = [], [], [], set()
    for index, (name, content) in enumerate(files, 1):
        digest = hashlib.sha256(content).hexdigest()
        if digest in seen:
            raise QPCRInputError(f"동일한 내용의 파일이 중복되었습니다: {name}")
        seen.add(digest)
        data, instrument = read_qpcr_upload(content, name)
        audit.append({"파일": name, "유효 웰": len(data),
                      "제외 내역": repr(data.attrs)})
        data = data.copy()
        data.attrs = {}
        data["파일"] = f"{index}: {name}"
        data["원본 시료명"] = data["Sample"].astype(str)
        data["입력 ID"] = str(index) + ":" + data["Sample"].astype(str)
        for row in data[["입력 ID", "원본 시료명", "Condition", "파일"]].drop_duplicates().itertuples(index=False, name=None):
            identity, sample, condition, source = row
            mapping.append({"입력 ID": identity, "파일": source, "원본 시료명": sample,
                            "실험 조건": sample if instrument else condition,
                            "생물학적 시료 ID": identity, "시료 코드": sample})
        frames.append(data)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(mapping), pd.DataFrame(audit)


def merge_samples(data, mapping):
    for column in ["실험 조건", "생물학적 시료 ID"]:
        if mapping[column].isna().any() or mapping[column].astype(str).str.strip().eq("").any():
            raise QPCRInputError(f"{column}에 빈 값이 있습니다.")
        mapping = mapping.copy()
        mapping[column] = mapping[column].astype(str).str.strip()
    if mapping.groupby("생물학적 시료 ID")["시료 코드"].nunique().gt(1).any():
        raise QPCRInputError("같은 생물학적 시료 ID에 서로 다른 시료 코드가 지정되어 있습니다.")
    if mapping.groupby("생물학적 시료 ID")["실험 조건"].nunique().gt(1).any():
        raise QPCRInputError("같은 생물학적 시료 ID에 서로 다른 조건이 지정되어 있습니다.")
    merged = data.merge(mapping[["입력 ID", "실험 조건", "생물학적 시료 ID"]], on="입력 ID", validate="many_to_one")
    merged["Sample"] = merged["생물학적 시료 ID"]
    merged["Condition"] = merged["실험 조건"]
    merged["Replicate"] = merged.groupby(["Sample", "Gene"]).cumcount() + 1
    codes = dict(zip(mapping["생물학적 시료 ID"], mapping["시료 코드"]))
    return merged[["Sample", "Condition", "Gene", "Replicate", "Ct"]], codes
