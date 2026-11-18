import unittest

import pandas as pd

from qpcr_analysis import QPCRInputError, analyze, parse_sample_code, read_qpcr_upload


def example():
    rows = []
    for sample, condition, reference, target in [
        ("C1", "Control", 20, 25),
        ("C2", "Control", 20, 25),
        ("T1", "Treatment", 20, 24),
        ("T2", "Treatment", 20, 24),
    ]:
        for gene, ct in [("Ref", reference), ("Target", target)]:
            for replicate in (1, 2, 3):
                rows.append((sample, condition, gene, replicate, ct))
    return pd.DataFrame(rows, columns=["Sample", "Condition", "Gene", "Replicate", "Ct"])


class TestQPCR(unittest.TestCase):
    def test_user_sample_code_mapping(self):
        self.assertEqual(parse_sample_code("00X"), ("0", "0", "스트레스 없음"))
        self.assertEqual(parse_sample_code("110H"), ("1", "10", "스트레스 있음"))
        self.assertIsNone(parse_sample_code("NTC"))
    def test_known_fold_change(self):
        result = analyze(example(), "Ref", "Control")
        treated = result.summary.query("Condition == 'Treatment' and Gene == 'Target'")
        self.assertAlmostEqual(treated["Geometric_Mean_FC"].iloc[0], 2.0)
        self.assertEqual(treated["Biological_N"].iloc[0], 2)

    def test_single_replicate_fails(self):
        data = example().query("not (Sample == 'C1' and Gene == 'Target' and Replicate != 1)")
        with self.assertRaisesRegex(QPCRInputError, "QC 실패"):
            analyze(data, "Ref", "Control")

    def test_missing_reference_fails(self):
        data = example().query("not (Sample == 'T1' and Gene == 'Ref')")
        with self.assertRaisesRegex(QPCRInputError, "누락"):
            analyze(data, "Ref", "Control")

    def test_sample_in_two_conditions_fails(self):
        data = example()
        data.loc[data.Sample == "T1", "Sample"] = "C1"
        with self.assertRaisesRegex(QPCRInputError, "여러 Condition"):
            analyze(data, "Ref", "Control")

    def test_nonnumeric_ct_fails(self):
        data = example()
        data["Ct"] = data["Ct"].astype(str)
        data.loc[0, "Ct"] = "Undetermined"
        with self.assertRaisesRegex(QPCRInputError, "Ct는"):
            analyze(data, "Ref", "Control")

    def test_instrument_setup_csv_reports_missing_ct(self):
        content = ("metadata,value\nWell,Well Position,Sample Name,Target Name\n"
                   "1,A1,S1,Ref\n").encode()
        with self.assertRaisesRegex(QPCRInputError, "Ct가 없습니다"):
            read_qpcr_upload(content, "setup.csv")

    def test_instrument_result_flags_undetermined_well(self):
        content = ("metadata,value\nWell,Well Position,Sample Name,Target Name,CT,Amp Status\n"
                   "1,F3,05X,NCED3,Undetermined,No Amp\n").encode()
        with self.assertRaisesRegex(QPCRInputError, "F3"):
            read_qpcr_upload(content, "result.csv")

    def test_ntc_without_amplification_is_reported_and_excluded(self):
        content = ("Well,Well Position,Sample Name,Target Name,Task,CT,Amp Status\n"
                   "1,A1,00X,PP2A,UNKNOWN,20.1,Amp\n"
                   "2,A2,05X,PP2A,NTC,Undetermined,No Amp\n").encode()
        result, instrument = read_qpcr_upload(content, "result.csv")
        self.assertTrue(instrument)
        self.assertEqual(result["Sample"].tolist(), ["00X"])
        self.assertEqual(result.attrs["ntc_wells"], ["A2"])

    def test_ntc_with_ct_blocks_analysis(self):
        content = ("Well,Well Position,Sample Name,Target Name,Task,CT,Amp Status\n"
                   "1,A1,00X,PP2A,UNKNOWN,20.1,Amp\n"
                   "2,A2,05X,PP2A,NTC,31.2,Amp\n").encode()
        with self.assertRaisesRegex(QPCRInputError, "NTC"):
            read_qpcr_upload(content, "result.csv")


if __name__ == "__main__":
    unittest.main()
