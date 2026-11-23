import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from qpcr_plot import PLOT_TYPES, make_figure, prepare_plot_data


class PlotTests(unittest.TestCase):
    def setUp(self):
        self.results = pd.DataFrame([
            {"Sample": "00X", "Condition": "00X", "Gene": "NCED3", "Fold_Change": 1.0},
            {"Sample": "05X", "Condition": "05X", "Gene": "NCED3", "Fold_Change": 1.5},
            {"Sample": "05H", "Condition": "05H", "Gene": "NCED3", "Fold_Change": 2.0},
            {"Sample": "15H", "Condition": "15H", "Gene": "NCED3", "Fold_Change": 0.5},
        ])

    def test_coded_factors_and_chart_formats(self):
        codes = {sample: sample for sample in self.results["Sample"]}
        plotted = prepare_plot_data(self.results, codes)
        self.assertEqual(set(plotted["온도"]), {"상온", "고온"})
        self.assertEqual(set(plotted["실험군"]), {"Col-1", "ost1"})
        for chart_type in PLOT_TYPES:
            fig = make_figure(plotted, "NCED3", chart_type, "ABA 코드", "온도", "실험군")
            self.assertEqual(len(fig.axes), 2)
            plt.close(fig)

    def test_arbitrary_condition_fallback(self):
        plotted = prepare_plot_data(self.results, {sample: "custom" for sample in self.results["Sample"]})
        self.assertNotIn("ABA 코드", plotted)
        fig = make_figure(plotted, "NCED3", "시료 점 + 기하평균", "Condition")
        self.assertEqual(len(fig.axes), 1)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
