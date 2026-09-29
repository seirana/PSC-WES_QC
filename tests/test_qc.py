"""Small synthetic examples; no patient-level data is needed for CI."""

import unittest

import numpy as np
import pandas as pd

from wes_qc import QCConfig, exact_hwe_pvalue, run_qc, variant_statistics
from wes_qc.core import mark_variant_pass


class StatisticsTests(unittest.TestCase):
    def test_maf_uses_less_common_allele_and_actual_group_sizes(self):
        matrix = np.array([[0, 1, 2, np.nan, 2], [np.nan] * 5])
        stats = variant_statistics(
            matrix, np.array([1, 1, 0, 0, 0], dtype=bool),
            np.array([0, 0, 1, 1, 1], dtype=bool),
        )
        self.assertEqual(stats.loc[0, "case_maf"], 0.25)
        self.assertEqual(stats.loc[0, "control_maf"], 0.0)
        self.assertEqual(stats.loc[0, "control_missing_rate"], 1 / 3)
        self.assertEqual(stats.loc[0, "all_minor_carriers"], 2)
        self.assertTrue(np.isnan(stats.loc[1, "control_maf"]))
        self.assertEqual(stats.loc[1, "case_missing_rate"], 1.0)

    def test_exact_hwe_is_symmetric_and_handles_no_calls(self):
        self.assertAlmostEqual(exact_hwe_pvalue(3, 1, 5), exact_hwe_pvalue(3, 5, 1))
        self.assertAlmostEqual(exact_hwe_pvalue(0, 2, 2), 3 / 35)
        self.assertEqual(exact_hwe_pvalue(0, 10, 0), 1.0)
        self.assertTrue(np.isnan(exact_hwe_pvalue(0, 0, 0)))
        with self.assertRaises(ValueError):
            exact_hwe_pvalue(-1, 2, 3)

    def test_invalid_genotype_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "0, 1, 2"):
            variant_statistics(np.array([[0, 3]]), np.array([1, 0]),
                               np.array([0, 1]))

    def test_duplicate_sites_detected_across_batches(self):
        metadata = pd.DataFrame({
            "chrom": ["1", "1"], "pos": [10, 10],
            "snp": ["a", "b"], "a0": ["A", "A"], "a1": ["C", "G"],
        })
        stats = variant_statistics(
            np.array([[0, 1, 2, 0], [0, 1, 2, 0]]),
            np.array([1, 1, 0, 0]), np.array([0, 0, 1, 1]),
        )
        result, _ = mark_variant_pass(stats, metadata, QCConfig(min_minor_carriers=0))
        self.assertEqual(result["duplicate_site"].tolist(), [True, True])
        self.assertFalse(result["pass_qc"].any())


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.matrix = np.array([
            [0, 1, 2, 0, 0, 1, 2, np.nan],
            [0, 1, 2, 0, 0, 1, 2, 0],
            [0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, np.nan],
            [0, 1, 2, 0, 0, 1, 2, np.nan],
            [0, 1, 2, 0, 0, 1, 2, np.nan],
        ], dtype=float)
        self.bim = pd.DataFrame({
            "chrom": ["1", "1", "1", "1", "X", "2"],
            "pos": [10, 10, 20, 30, 40, 50],
            "snp": [f"rs{i}" for i in range(6)],
            "a0": ["A"] * 6, "a1": ["G"] * 6,
        })
        self.fam = pd.DataFrame({
            "fid": ["f"] * 8, "iid": [f"s{i}" for i in range(8)],
        })
        self.groups = self.fam.copy()
        self.groups["group"] = ["case"] * 4 + ["control"] * 4

    def test_chunked_two_pass_qc_filters_missing_sample(self):
        config = QCConfig(max_variant_missing_rate=0.25,
                          min_sample_call_rate=0.75,
                          min_minor_carriers=1)
        result = run_qc(
            self.matrix, self.bim, self.fam, self.groups,
            config, batch_size=2,
        )
        self.assertEqual(result.preliminary_variant_count, 3)
        self.assertEqual(result.samples["pass_qc"].tolist(), [True] * 7 + [False])
        self.assertFalse(result.variants.loc[0, "pass_qc"])
        self.assertFalse(result.variants.loc[1, "pass_qc"])
        self.assertFalse(result.variants.loc[2, "pass_qc"])
        self.assertFalse(result.variants.loc[4, "hwe_applied"])
        self.assertEqual(result.variants.loc[5, "control_missing_rate"], 0)

    def test_missing_group_mapping_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "Every FAM sample"):
            run_qc(self.matrix, self.bim, self.fam, self.groups.iloc[:-1],
                   batch_size=2)

    def test_mismatched_metadata_fails(self):
        with self.assertRaisesRegex(ValueError, "BIM fields/rows"):
            run_qc(self.matrix, self.bim.iloc[:-1], self.fam, self.groups)


if __name__ == "__main__":
    unittest.main()
