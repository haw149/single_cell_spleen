from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from expression_views import expression_frame, expression_summary


class MultiGeneTests(unittest.TestCase):
    def test_filtered_cells_keep_expression_alignment_and_gene_identity(self):
        obs = pd.DataFrame({"cell_id": ["c3", "c1"], "cluster": ["A", "A"]}, index=[3, 1])
        frame = expression_frame({"GeneA": np.array([10, 0, 30, 4]),
                                  "GeneB": np.array([50, 6, 70, 8])}, obs, "cluster")
        self.assertEqual(frame.cell_id.tolist(), ["c3", "c1", "c3", "c1"])
        self.assertEqual(frame.expression.tolist(), [4, 0, 8, 6])
        self.assertEqual(frame.index.tolist(), [3, 1, 3, 1])
        summary = expression_summary(frame).set_index("gene")
        self.assertEqual(summary.loc["GeneA", "mean"], 2)
        self.assertEqual(summary.loc["GeneB", "mean"], 7)
        self.assertEqual(summary.loc["GeneA", "percent_positive"], 50)
        self.assertEqual(summary.loc["GeneB", "cells"], 2)

    def test_empty_selection_and_missing_group(self):
        obs = pd.DataFrame({"cell_id": ["c0"]})
        self.assertTrue(expression_frame({}, obs, "cluster").empty)
        frame = expression_frame({"GeneA": [2]}, obs, "cluster")
        self.assertEqual(frame.group.tolist(), ["All cells"])
        self.assertTrue(expression_summary(frame.iloc[:0]).empty)
