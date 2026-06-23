import unittest

from omnomnomics.cli import (
    describe_public_steps,
    map_public_to_internal_steps,
    resolve_public_mode_steps,
)


class StepMappingTests(unittest.TestCase):
    def test_rna_public_steps_are_internal_steps(self):
        self.assertEqual(resolve_public_mode_steps("all", 12), list(range(1, 13)))
        self.assertEqual(map_public_to_internal_steps("RNA", [10, 11, 12]), [10, 11, 12])

    def test_atac_public_followup_steps_map_to_internal_rules(self):
        public_steps = [10, 11, 12, 13, 14, 15]
        self.assertEqual(map_public_to_internal_steps("ATAC", public_steps), [10, 11, 13, 14, 15, 16])

    def test_chip_public_followup_steps_map_to_internal_rules(self):
        public_steps = [10, 11, 12, 13, 14, 15]
        self.assertEqual(map_public_to_internal_steps("CHIP", public_steps), [10, 11, 13, 14, 15, 16])

    def test_post_de_public_step_description_hides_internal_rule(self):
        description = describe_public_steps("ATAC", [15])
        self.assertIn("15 (analyze differential peaks post-DE)", description)
        self.assertNotIn("internal", description)
        self.assertNotIn("16", description)


if __name__ == "__main__":
    unittest.main()
