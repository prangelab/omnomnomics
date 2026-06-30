import argparse
import contextlib
import io
import unittest
from unittest import mock

from omnomnomics.cli import (
    FASTP_NEXTERA_R1,
    FASTP_NEXTERA_R2,
    resolve_fastp_adapter_settings,
)


def args(**overrides):
    values = {
        "fastp_adapter_mode": None,
        "fastp_adapter_sequence": None,
        "fastp_adapter_sequence_r2": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


class FastpAdapterSettingsTests(unittest.TestCase):
    def test_assay_mode_resolves_atac_to_nextera(self):
        _, mode, r1, r2 = resolve_fastp_adapter_settings(args(), {}, "ATAC", True)

        self.assertEqual(mode, "nextera")
        self.assertEqual(r1, FASTP_NEXTERA_R1)
        self.assertEqual(r2, FASTP_NEXTERA_R2)

    def test_assay_mode_resolves_rna_and_chip_to_overlap(self):
        for assay in ("RNA", "CHIP"):
            with self.subTest(assay=assay):
                _, mode, r1, r2 = resolve_fastp_adapter_settings(args(), {}, assay, True)
                self.assertEqual(mode, "overlap")
                self.assertEqual(r1, "")
                self.assertEqual(r2, "")

    def test_explicit_paired_mode_requires_read2_adapter(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), mock.patch("sys.exit", side_effect=SystemExit) as mocked_exit:
            with self.assertRaises(SystemExit):
                resolve_fastp_adapter_settings(
                    args(
                        fastp_adapter_mode="explicit",
                        fastp_adapter_sequence="ACGT",
                    ),
                    {},
                    "RNA",
                    True,
                )

        mocked_exit.assert_called_once_with(1)


if __name__ == "__main__":
    unittest.main()
