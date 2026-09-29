"""Synthetic end-to-end command test without a real patient file or PLINK install."""

import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from wes_qc.cli import main


class CLITests(unittest.TestCase):
    def test_writes_outputs_and_input_checksums(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root / "cohort"
            for suffix in ("bed", "bim", "fam"):
                Path(f"{prefix}.{suffix}").write_bytes(suffix.encode())
            group_path = root / "groups.csv"
            groups = pd.DataFrame({
                "fid": ["f"] * 4, "iid": [f"s{i}" for i in range(4)],
                "group": ["case", "case", "control", "control"],
            })
            groups.to_csv(group_path, index=False)
            bim = pd.DataFrame({
                "chrom": ["1", "1", "1"], "pos": [1, 2, 3],
                "snp": ["v1", "v2", "v3"], "a0": ["A"] * 3, "a1": ["G"] * 3,
            })
            fam = groups[["fid", "iid"]]
            matrix = np.array([
                [0, 1, 0, 1],
                [1, 2, 1, 2],
                [0, 1, 2, np.nan],
            ])
            fake_plink = types.ModuleType("pandas_plink")
            fake_plink.read_plink = lambda *_args, **_kwargs: (bim, fam, matrix)
            output = root / "private-results"
            with patch.dict(sys.modules, {"pandas_plink": fake_plink}):
                code = main([
                    "--plink-prefix", str(prefix), "--groups", str(group_path),
                    "--output", str(output), "--min-minor-carriers", "1",
                ])
            self.assertEqual(code, 0)
            variants = pd.read_csv(output / "variants.csv")
            samples = pd.read_csv(output / "samples.csv")
            metadata = json.loads((output / "run_metadata.json").read_text())
            self.assertEqual(variants["pass_qc"].tolist(), [True, True, False])
            self.assertTrue(samples["pass_qc"].all())
            self.assertEqual(metadata["n_final_variants"], 2)
            self.assertEqual(len(metadata["input_sha256"]["bed"]), 64)
            self.assertTrue(Path(f"{prefix}.bed").exists())


if __name__ == "__main__":
    unittest.main()
