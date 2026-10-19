import pyBigWig
import pysam

import test_chip_inputs as fixtures


class CoverageResumeTests(fixtures.LibraryFixture):
    run_workflow = fixtures.SnakemakePreprocessingTests.run_workflow

    def prepare(self, assay):
        bam = self.bam('a.sorted.dups_marked.filtered.bam', folder='filtered_BAM')
        if assay == 'rna':
            replacement = bam.with_suffix('.stranded.bam')
            with pysam.AlignmentFile(str(bam)) as source:
                with pysam.AlignmentFile(str(replacement), 'wb', header=source.header) as dest:
                    for index, read in enumerate(source):
                        read.flag = 16 if index % 2 else 0
                        dest.write(read)
            replacement.replace(bam)
        pysam.index(str(bam))
        config = self.cli_config([fixtures.row('a', 'chip')], '8-9', assay=assay)
        tracks = self.project / 'BigWigs'
        tracks.mkdir()
        names = ['a.plus.bw', 'a.minus.bw'] if assay == 'rna' else ['a.bw']
        for name in names:
            with pyBigWig.open(str(tracks / name), 'w') as handle:
                handle.addHeader([('chr1', 1000)])
                handle.addEntries(['chr1'], [100], ends=[130], values=[1.0])
        hubs = self.project / 'merged_hubs'
        hubs.mkdir()
        (hubs / 'extra_9.tmp').write_text('complete\n')
        return config, [tracks / name for name in names]

    def test_rna_reuses_both_strands_without_marker_and_rebuilds_missing_strand(self):
        config, tracks = self.prepare('rna')
        result = self.run_workflow(config)
        self.assertIn('Nothing to be done', result)
        tracks[1].unlink()
        result = self.run_workflow(config)
        self.assertIn('rule create_wiggles:', result)
        self.assertIn('rule merge_wiggles:', result)
        self.run_workflow(config, dry_run=False)
        self.assertTrue(all(path.stat().st_size > 0 for path in tracks))
        self.assertFalse((self.project / 'BigWigs/a.extra_8.tmp').exists())
        self.assertIn('Nothing to be done', self.run_workflow(config))
        for path in tracks:
            with pyBigWig.open(str(path)) as handle:
                self.assertEqual(handle.chroms(), {'chr1': 1000})

    def test_atac_reuses_track_without_marker_and_rebuilds_missing_track(self):
        config, tracks = self.prepare('atac')
        self.assertIn('Nothing to be done', self.run_workflow(config))
        tracks[0].unlink()
        self.assertIn('rule create_wiggles:', self.run_workflow(config))
