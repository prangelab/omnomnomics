<p align="center">
	<img src="https://github.com/prangelab/omnomnomics/assets/157825254/f8f83873-7c53-495e-9c54-33eece99a8db" width="700" alt="DALL·E 2024-06-23 12 04 41 - A futuristic character resembling the Cookie Monster, eating a double-stranded DNA helix  The background is high-tech, with neon lights and circuit pa">
</p>
<div align="center">

# <center>***OMNOMNOMICS***</center>

</div>

_Omnomnomics is an A-Z processing NGS pipeline for RNA-, ChIP-, and ATAC-seq data. It trims FASTQ files, runs FastQC, aligns reads, performs assay-aware downstream analysis, and records the resulting workflow provenance._

## Quickstart:
In a rush? Once the environment is installed, you should be able to run:
```
omnomnomics rna -i path/to/your/experiment_dir -g genome_to_use
```
Other assay entrypoints use the same verb style:
```bash
omnomnomics atac -i path/to/your/experiment_dir -g genome_to_use
omnomnomics chip -i path/to/your/experiment_dir -g genome_to_use
```
For more information about how to use Omnomnomics, see the usage section down below. 

Packaged helper verbs are also available:
```bash
omnomnomics --version
omnomnomics genomes --help
omnomnomics create-track-color-table -h
omnomnomics display-track-color-table -h
```

The packaged `color_data_for_hubs` directory is intended as a read-only source of default palettes. When creating a custom palette, provide a writable target directory with `-P`.

## Installation

Clone the repository:

```bash
git clone https://github.com/prangelab/omnomnomics.git path/to/install/omnomnomics
cd path/to/install/omnomnomics
```

Create the main `omnomnomics` environment and install the package:

```bash
micromamba env create -f environment.yml
micromamba activate omnomnomics
pip install -e .
```

The main environment includes the MEME Suite tools used by post-DE peak motif analysis. IDR and SPP use small companion environments because their dependency stacks conflict with the main Python/R analysis environment.

### IDR Companion Environment

Default narrow ATAC/ChIP peak calling uses IDR. Bioconda IDR currently requires an older Python than the main `omnomnomics` environment, so IDR is installed in a small companion environment. If you want to use the default `--narrow-peak-strategy idr` mode for ATAC or narrow ChIP, run this once after creating the main environment:

```bash
bash scripts/install_idr_helper.sh
```

The helper creates or updates `omnomnomics-idr` from `environment.idr.yml`, then installs an `idr` wrapper into the main `omnomnomics` environment. After this, users should continue activating only the main environment:

```bash
micromamba activate omnomnomics
```

You do not need to activate `omnomnomics-idr` manually during normal pipeline runs. If you only use `--narrow-peak-strategy macs3`, the IDR companion environment is not required.

### SPP Companion Environment

ATAC/ChIP peak QC can report strand cross-correlation metrics (`NSC` and `RSC`) through `phantompeakqualtools` / SPP. Current Bioconda SPP packages require an older R stack than the main `omnomnomics` environment, so install it through the helper after creating the main environment:

```bash
bash scripts/install_spp_helper.sh
```

The helper creates or updates `omnomnomics-spp` from `environment.spp.yml`, then installs a `run_spp.R` wrapper into the main `omnomnomics` environment. Users should still activate only the main environment during normal runs:

```bash
micromamba activate omnomnomics
```

If the SPP helper is not installed, peak QC still runs and records the SPP fields as unavailable. Install the helper when you want `NSC`/`RSC` metrics or `--spp-gate` decisions.

Reference genomes are managed separately from the code checkout. Install them under the configured genome assembly root, or use the packaged genome helper commands. A normalized assembly layout contains at least:
- `fasta/genome.fa`
- `annotation/genes.gtf`
- `star/`
- `hisat2/`
- `aux/`

Useful genome helper commands:

```bash
# Search remote assemblies from the configured provider
omnomnomics genomes list --species human
omnomnomics genomes ls --species mouse

# Show locally installed normalized assemblies
omnomnomics genomes installed
omnomnomics genomes local --species mouse

# Install the latest matching assembly for a species
omnomnomics genomes install --species mouse

# Backfill or refresh an ENCODE blacklist BED for an installed assembly
omnomnomics genomes blacklist --assembly GRCh38

# Cache the default MEME motif database used by ATAC/ChIP post-DE motif analysis
omnomnomics genomes motifs
```

When available through genomepy, `omnomnomics genomes install` also caches the ENCODE blacklist BED into the normalized assembly `aux/` directory. Omnomnomics disables genomepy STAR/HISAT2 plugins during install and builds requested aligner indexes directly from the normalized FASTA/GTF, so index layout stays predictable. The separate `genomes blacklist` command is mainly for older genome installs or explicit refreshes.
STAR genome indexing can require substantially more memory than HISAT2 indexing. On HPC systems where memory scales with requested CPU cores, request a larger CPU slice for whole-genome installs.
Genome installs also cache the default JASPAR MEME-format motif database under the configured genome assembly root in `motif_databases/`. Post-DE motif analysis reuses this permanent cache and will try to create it on first run if it is missing. Use `omnomnomics genomes motifs --force` to refresh the cached database, or set `post_de_motif_database` to an explicit MEME-format file for strict database pinning.

## Overview of Omnomnomics
OMNOM_HOME, is the directory containing _Omnomnomics_. This is where you installed _Omnomnomics_ and where the wrapper script, main snakefile, bin folder, genomes and slurm_profile are. The default directory structure generated and used by the pipeline is the EXPERIMENT_DIR. This is the main input directory containing your experiment and all your files folders for the in- and output. 
The FASTQ folder (inside your EXPERIMENT_DIR) contains raw FASTQ files. These are used by the trimming step. The trimmed_FASTQ folder contains adapter-trimmed FASTQ files. 
The trimmed reads are then used by FastQC and the mapping jobs. The fastqc_reports folder contains the QC reports (html and zip). 
The BAM folder contains raw BAM files generated by the mapping jobs. These are then used by touchup_bam to produce filtered BAM files in filtered_BAM. 
The filtered BAM files are then used by BAM indexing, BAM statistics, optional HOMER tag directory generation, BigWig creation, peak calling, and read counting.
The HOMER_tagDirs folder contains optional HOMER tag directory tarballs for users who want to run HOMER peak calling downstream.
The BigWigs folder contains per-sample BigWig files generated directly from BAMs. For RNA runs, each sample is written as a stranded plus/minus BigWig pair so the UCSC hub can show signal above and below the zero axis. The merged_hubs folder contains grouped UCSC track hubs assembled from those BigWigs.
The DE_calling folder contains count tables and DESeq2 outputs. For RNA this is the step-12 DE branch. ATAC and ChIP use the same metadata-driven DE-style framework with assay-specific feature definitions. The DE workflow renders a reproducible run script plus a customization guide script, writes resolved metadata, and generates a DE results tree under `DE_calling/<de_out_dir>`.
Lastly, the peak_calling folder contains peak or feature BED files for ChIP-seq and ATAC-seq experiments generated from filtered BAM files. In ChIP `--broad-mode genebody`, this folder contains annotation-derived gene-body features rather than MACS3 peaks.

## Usage: 
_Omnomnomics_ is a Snakemake pipeline with a Python CLI, packaged workflow config, site config, and a set of Snakemake rules. The CLI performs safety checks, prepares a run config, and submits the Snakemake controller as a SLURM job. The main Snakefile prepares threads and memory requirements for the rules, includes all rules, and sets up the workflow. The workflow config captures pipeline defaults, while the site config captures HPC-specific defaults. The rules perform the actual tool execution and data processing.

To invoke _Omnomnomics_, run:
```bash
omnomnomics rna -i path/to/your/experiment_dir -g genome_to_use
```

## Directory Structures
```
INSTALLATION
|
|- environment.yml
|- environment.idr.yml
|- environment.spp.yml
|- scripts/
|- src/omnomnomics/
     |- cli.py
     |- workflow/ (packaged configuration, rules, templates, R code, and Slurm profile)

REFERENCE_ROOT
|
|- assemblies
     |- <assembly_name>
          |- fasta/genome.fa
          |- annotation/genes.gtf
          |- star/
          |- hisat2/
          |- aux/
|- motif_databases/



EXPERIMENT_DIR
|
|- run_configs
	|- omnomnomics.run."run_date".config.yaml
	|- any backups with random numbers if multiple runs were made on same day
|- slurm_logs
		|- controller
			|- controller.job_ID.out
		|- "rule_name1"
		|- "rule_name1".job_ID.out
		|- .... 
	|- "rule_name2"
		|- "rule_name2".job_ID.out
		|- .... 
	|- "rule_name3"
		|- "rule_name3".job_ID.out
		|- ....
	|- ....
|- run_logs
	|- omnomnomics.run."run_date".log
	|- any backups with random numbers if multiple runs were made on same day
|- MultiQC
	|- omnomnomics.run.{run_date}.multiqc_report.html
|- FASTQ
     |- {sample}_R1.fastq.gz
     |- {sample}_R2.fastq.gz (if paired)
|- trimmed_FASTQ
     |- {sample}_R1.trimmed.fastq.gz (if paired else {sample}.trimmed.fastq.gz)
     |- {sample}_R2.trimmed.fastq.gz (if paired else None)
|- fastqc_reports
     |- {sample}_R1.trimmed_fastqc.html (if paired else {sample}.trimmed_fastqc.html)
     |- {sample}_R2.trimmed_fastqc.html (if paired else None)
     |- {sample}_R1.trimmed_fastqc.zip  (if paired else {sample}.trimmed_fastqc.zip)
     |- {sample}_R2.trimmed_fastqc.zip  (if paired else None)
|- BAM
     |- {sample}.bam (from aligners)
     |- {sample}.unpaired.aligned.bam   (from hisat2 if keepunpaired = True
     |- {sample}.unpaired.unaligned.bam (from hisat2 if keepunpaired = True
     |- {sample}.bam 			(after lane numbers are merged)
     |- {sample).HISAT2_stats.txt 	(if map tool = hisat2)
     |- {sample).STAR_TE_stats.txt 	(if map tool = star_te)
     |- {sample).STAR_stats.txt 	(if map tool = star)
|- filtered_BAM
     |- {sample}.sorted.dups_marked.filtered.bam (if the type of data = ATAC or RNA, from touchup_bam)
     |- {sample}.filtered.bam (if the type of data = ChIP, from touchup_bam)
     |- {sample}.ATAC_stats.txt (if the type of data = ATAC, from touchup_bam)
     |- {sample}.sorted.dups_marked.filtered.bam.bai (if the type of data = ATAC or RNA, from index_bam)
     |- {sample}.filtered.bam.bai (if the type of data = ChIP, from index_bam)
     |- {sample}.sorted.dups_marked.filtered.bam.stats.txt (if the type of data = ATAC or RNA, from bam_stats)
     |- {sample}.filtered.bam.stats.txt (if the type of data = ChIP, from bam_stats)
|- HOMER_tagDirs
     |- {sample}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz (if the type of data = ATAC or RNA, from make_HOMER_tagDIR)
     |- {sample}.filtered.HOMER_tagDir.tar.gz (if the type of data = ChIP, from make_HOMER_tagDIR)
|- BigWigs
     |- {sample}.plus.bw (if the type of data = RNA, from create_wiggles)
     |- {sample}.minus.bw (if the type of data = RNA, from create_wiggles)
     |- {sample}.bw (if the type of data = ATAC, from create_wiggles)
     |- {sample}.bw (if the type of data = ChIP, from create_wiggles)
|- merged_hubs
     |- {group_name}.hub 
			|- "genome used"
			     |- {sample}.plus.bw (for RNA)
			     |- {sample}.minus.bw (for RNA)
			     |- {sample}.bw (for ATAC and ChIP)
			     |- trackDb.txt
			|- genomes.txt
			|- hub.txt
|- peak_calling
     |- {sample}.MACS3.q-0p01_peaks.bed (if the data type = ATAC, from call_peaks)
     |- all_groups.merged_peaks.bed (if the data type = ATAC, from call_peaks)
     |- THENAME.MACS3.q-0p05_peaks.bed (if the data type = ChIP, from call_peaks)
     |- THENAME.MACS3.q-0p01_peaks.bed (if the data type = ChIP, from call_peaks)
     |- THENAME.MACS3.q-0p001_peaks.bed (if the data type = ChIP, from call_peaks)
     |- peak_qc
          |- {chip|atac}.peak_qc_metrics.tsv
          |- {chip|atac}.sample_frip_metrics.tsv
          |- {chip|atac}.sample_qc_metrics.tsv
          |- {chip|atac}.peak_qc_summary.pdf
          |- {sample}.cross_correlation.tsv
          |- {sample}.cross_correlation.pdf

Peak-QC library complexity, SPP cross-correlation, and FRiP metrics use deterministic alignment caps by default to keep large ATAC/ChIP BAM diagnostics scalable. Defaults are 5,000,000 alignments for complexity and 10,000,000 alignments for SPP/FRiP; set `--library-complexity-max-reads 0`, `--spp-max-reads 0`, or `--frip-max-reads 0` to scan full BAMs. Library complexity prefers the merged pre-deduplication BAM and records a filtered-BAM fallback explicitly when that input is unavailable. BigWig generation does not subsample BAMs and uses CPM normalization.
|- DE_calling
     |- basename.EXPERIMENT_DIR.raw_read_quant.table.txt (from count_reads)
     |- metadata_derived.tsv (resolved metadata used by step 12)
     |- DE_analysis.rendered.R (run-specific reproducible script)
     |- DE_analysis.customization_guide.R (run-specific script with documented customization recipes)
     |- basename.EXPERIMENT_DIR.results.zip (step-12 archive with scripts + results tree)
     |- qc (shared QC outputs for step-12 runs)
     |- <de_out_dir> (default: results; override with --de-out-dir)
          |- <contrast_1>
          |- <contrast_2>
          |- contrast_summary.tsv
          |- enrichment_summary.tsv
          |- contrast_summary.sig_up_down_barplot.pdf
     |- analyze_peaks_de (ATAC/ChIP post-DE interpretation)
          |- sets
          |- signal
          |- motifs
          |- summary

```
Peak QC notes:
- `FRiP` and `peak_count` are reported per MACS3 peak set, and per-sample FRiP is reported for each selected peak set.
- `NRF`, `PBC1`, and `PBC2` are reported per BAM as library complexity metrics.
- `NSC` and `RSC` are reported from strand cross-correlation analysis via `phantompeakqualtools` / SPP when `scripts/install_spp_helper.sh` has installed the `run_spp.R` wrapper.
- `NSC` and `RSC` are most informative for `TF / narrow peaks` and should be interpreted more cautiously for broad histone marks.

## Workflow:
```
Pre-processing steps:
1:	PREP: 	Trim FASTQ files
2:	QC:   	Run FastqQC
3:	MAP:  	Align reads to the genome
4:	MERGE:	Merge flowcell lanes (if needed; i.e., if samples were spread over multiple lanes (L00n) to meet the required seq depth)

Post-processing steps:
5:	Touchup BAM
			Runs a chain of samtools commands to collate, repair mate tags, sort, mark or remove duplicates, and filter the BAM files
			Specifically: samtools collate | fixmate -m | sort | markdup | view
				RNA:	keep duplicates by default, MAPQ >= 15
				ATAC:	remove duplicates by default, MAPQ >= 30
				ChIP:	remove duplicates by default, MAPQ >= 30
				Override duplicate handling with --remove-duplicates or --keep-duplicates
				ATAC:	Get chrM stats
6:	Index BAM
7:	Generate pre/post-filter alignment QC summaries
			Writes a per-sample tabular QC summary plus per-sample PDF/SVG summary plots
8:	Create BigWigs
9:	Merge Bigwigs and trackhubs by experimental group. Optional, will even in 'auto' mode only be run if the -E flag is set. See below.

Assay dependent follow-up steps use assay-specific public numbering:

RNA public steps:
10:	Create trackhubs
11:	Create count table
12:	Call DE genes (DESeq2)

ATAC public steps:
10:	Call peaks
11:	Peak QC
12:	Analyze peaks (pre-DE)
13:	Create count table
14:	Call differential chromatin regions
15:	Analyze differential peaks (post-DE)

ChIP public steps:
10:	Call peaks
11:	Peak QC
12:	Analyze peaks (pre-DE)
13:	Create count table
14:	Call differential chromatin regions
15:	Analyze differential peaks (post-DE)

Optional export:
	--create-homer-tagdirs:	Create HOMER tag directory tarballs alongside the main numbered workflow outputs.
					This requires a working HOMER installation with the requested genome installed separately via `configureHomer.pl`.


Auto mode:
		By default runs the whole public pipeline. i.e., sets mode to 'all'
		Will detect an aborted run (e.g. cancelled by user using scancel or requeued by slurm due to resouce constraints)
		If an aborted run is detected, 'auto' mode will restart the run after the last succesfully completed step.
		Optional HOMER tag directory export is not part of the numbered pipeline. Add `--create-homer-tagdirs` if you want it.
		If MultiQC is enabled, the final reporting block also writes an experiment-level QC TSV plus PDF/SVG summary plots before running MultiQC. These aggregate raw FASTQ read counts, trimmed FASTQ read counts, mapper-reported alignment metrics, and the post-touchup BAM QC summaries.

Some job mode examples:
	-j all:	Run the whole pipeline
					Input is the 'FASTQ' folder inside your <EXPERIMENT_DIR>, which contains the .fastq.gz files you want analysed.
	-j 1-3:	Run the pipeline up to and including the mapping step.
					Input is the 'FASTQ' folder inside your <EXPERIMENT_DIR>, which contains the .fastq.gz files you want analysed.
	-j 1-9:	Run the whole pipeline excluding the assay dependent follow-up.
					Input is the 'FASTQ' folder inside your <EXPERIMENT_DIR>, which contains the .fastq.gz files you want analysed.
	-j 4-11:	Run The whole pipeline, starting after the mapping step.
					Input is the 'BAM' folder inside your <EXPERIMENT_DIR>, which contains the raw .bam files you want analysed.
	-j 7:	Only run the alignment QC summary step.
					Input is the 'BAM' and 'filtered_BAM' folders inside your <EXPERIMENT_DIR>, which contain the pre- and post-filter BAM files to compare.
	-j 1-8,10-12:	Run everything but do not create merged trackhubs.
					Input is the 'FASTQ' folder inside your <EXPERIMENT_DIR>, which contains the .fastq.gz files you want analysed.
	-j 1,3,5,6,10-12:	Trim the reads, map, touch up, create index, and run the assay dependent follow-up. Skip QC, merging, and BigWigs.
					Input is the 'FASTQ' folder inside your <EXPERIMENT_DIR>, which contains the .fastq.gz files you want analysed.
	-j 1-12 --create-homer-tagdirs:	Run the full public RNA pipeline and also export optional HOMER tag directories.
					Input is the 'FASTQ' folder inside your <EXPERIMENT_DIR>, which contains the .fastq.gz files you want analysed.
```

## Command Line Options:

    -i <EXPERIMENT_DIR>:    Main directory containing your experiment.
                                    Required argument
    ASSAY VERB:             Select assay via CLI verb: `rna`, `atac`, or `chip`.
                                    Required for workflow runs.
	  -g GENOME:				Genome version. Avaliable versions:
								UCSC/RefSeq/NCBI:	mm10, mm39, hg38.
								ENSMBL/GenBank:		GRCh38.p14, GRCm39
								Required argument
    -X:                     eXclude multiQC stats aggregator. Set if you don not wish to run multiQC.
    --create-homer-tagdirs: Create optional HOMER tag directory tarballs in addition to the main numbered workflow outputs.
                                    HOMER genomes are not installed automatically with the package. Install them separately with `configureHomer.pl` if you want to use this export.
    --rerun-selected-steps: Force recomputation of the selected workflow steps by deleting their current outputs first.
                                    Default behavior is to reuse existing outputs when Snakemake sees them as up to date.
    --site-config:         Optional path to a site-specific config YAML.
                                    Default resolution order:
                                    1. $XDG_CONFIG_HOME/omnomnomics/site.yaml
                                    2. ~/.config/omnomnomics/site.yaml
                                    3. packaged site config
    --retention-policy:   Post-run output retention policy.
                                    all: keep all pipeline outputs
                                    pruned: keep FASTQ plus reusable downstream outputs
                                    minimal: keep FASTQ plus only the requested terminal outputs
                                    Default: all
    --max-project-size:   Soft project-size cap such as 200G or 800GB.
                                    Omnomnomics may delete safe intermediates and skip BigWig or trackhub creation if the cap would otherwise be exceeded.
                                    Intermediate cleanup is stage-aware:
                                    trimmed_FASTQ is only eligible for cleanup after step 4 fully finishes.
                                    BAM is only eligible for cleanup after step 5 fully finishes.
    --remove-duplicates:    Remove duplicate reads in step 5. Default is assay-aware: keep for RNA, remove for ATAC and ChIP.
    --keep-duplicates:      Keep duplicate reads in step 5. Default is assay-aware: keep for RNA, remove for ATAC and ChIP.
    -j MODE:                Job mode. Can be 'auto', 'all' or a range of jobs. See below (-h) for some 
							examples.
                                    Default: auto
    -T TOOL:                Trimming tool choice. Can be skewer or fastp
                                    Default: skewer
    --fastp-adapter-mode:   fastp adapter handling mode. Choices: assay, overlap, auto_detect,
                            nextera, truseq, explicit, off.
                                    Default: assay
                            `assay` resolves to explicit Nextera adapters for ATAC and fastp's
                            paired-end overlap trimming for RNA/ChIP. Use `auto_detect` only when
                            you specifically want fastp's paired-end adapter auto-detection.
    --fastp-adapter-sequence:
                            Read 1 adapter sequence when --fastp-adapter-mode explicit is used.
    --fastp-adapter-sequence-r2:
                            Read 2 adapter sequence when --fastp-adapter-mode explicit is used.
                            fastp is available as an opt-in trimmer, but skewer is the default because validation on HPC test data showed intermittent fastp worker hangs before any useful fastp log output.
    -M TOOL:                Mapping tool choice. Can be HISAT2, STAR, or STAR_TE. STAR(_TE) can only be used 
							for RNA-seq data.
                                    Default: HISAT2
    --de-formula:           Explicit DESeq2 design formula for DE calling.
                                    If provided, it overrides `--de-columns` and `--de-block`.
    -I <INPUT>:             Input BAM file used for ChIP peak calling with MACS3.
                                    Default: Do not use input.
    -m <METADATA_TABLE>:    Tabular metadata file describing all samples.
                                    First column should be `filename`.
                                    Omnomnomics derives `sample_id`, `sample_type`, and `sample_color` internally.
    --sample-name:          Metadata columns used to derive `sample_id`.
                                    Accepts comma-separated column names or 1-based indices.
    --sample-type:          Metadata columns used to derive `sample_type`.
                                    Accepts comma-separated column names or 1-based indices.
                                    Default: all samples grouped together.
    --sample-color:         Metadata columns used to derive `sample_color` categories.
                                    Accepts comma-separated column names or 1-based indices.
                                    Actual colors are still assigned through track color tables.
    --de-columns:           Metadata columns of interest for auto-building DESeq2 designs.
    --de-block:             Metadata columns to include as blocking terms in auto-built designs.
    --de-interactions:      Include interaction terms when auto-building the DESeq2 design from exactly two `--de-columns`.
    --broad-mode:           ChIP broad-mark handling mode.
                                    off: narrow / TF-like path
                                    domain: MACS3 broad-domain path
                                    genebody: annotated gene-body feature path
                                    diffuse: fixed-bin tiled-genome path for very diffuse broad marks
    --chip-broad-qvalue:    Relaxed MACS3 q-value for ChIP `domain` mode pooled, replicate, and pooled-pseudoreplicate calls.
                                    Default: 0.05
    --chip-broad-cutoff:    MACS3 `--broad-cutoff` for ChIP `domain` mode.
                                    Default: 0.1
    --chip-broad-min-length:
                                    Optional MACS3 `--min-length` for ChIP `domain` mode.
    --chip-broad-max-gap:   Optional MACS3 `--max-gap` for ChIP `domain` mode.
    --chip-broad-replicate-fraction:
                                    Minimum fraction of true replicates that must support a pooled broad domain in ChIP `domain` mode.
                                    Default: 1.0
    --chip-broad-overlap-fraction:
                                    Minimum overlap fraction used when evaluating support for pooled broad domains in ChIP `domain` mode.
                                    Default: 0.5
    --chip-diffuse-bin-size:       Fixed bin size in bp for ChIP `diffuse` mode.
                                    Default: 10000
    --chip-diffuse-merge-gap:      Maximum gap in bp used to merge adjacent significant diffuse bins into domains.
                                    Default: one bin width
    --narrow-peak-strategy:
                                    Narrow-peak strategy for ATAC and narrow ChIP.
                                    `idr` uses replicate-aware consensus, `macs3` uses parameter optimization.
    --idr-mode:             IDR mode for narrow ATAC and narrow ChIP.
                                    `basic` uses true replicate IDR, `encode` adds pseudoreplicate diagnostics.
    --idr-pair-fraction:    Minimum fraction of replicate pairs that must support a peak in IDR consensus.
                                    Default: 0.5
    --idr-pairing-policy:   Replicate pairing policy for IDR when groups have more than two replicates.
                                    `all_pairs` or `anchor_vs_all`
    --idr-min-input-peaks:  Minimum narrowPeak rows required before attempting an IDR comparison.
                                    Sparse comparisons are skipped and reported; if no true-pair IDR comparison is usable, pooled MACS3 peaks are used as a flagged fallback.
                                    Default: 20
    --spp-gate:             SPP QC gate mode for ATAC and ChIP peak-QC step.
                                    `none`, `warn`, `drop`, or `strict`
    --post-de-signal-policy:
                                    Policy for ATAC/ChIP post-DE heatmaps and profiles.
                                    auto: schedule missing BigWigs through step 8 and then plot them
                                    require: use existing BigWigs and fail clearly if they are absent
                                    skip: disable post-DE signal plotting and record explicit skips
                                    Default: auto
    -a STR:                 Appendix to add to track name
                                    Default: hub
    -C <Color_Table_FILE>:  File specifying which colors to use for the tracks
                          	  Default: gray.tint.color.table from the packaged palette directory
                            Can be a *txt list file with one color table per line. Different color tables will 
							be used per hub group where applicable.
                            Can be a full (relative) path to a file or a file basename only in conjuction with 
	  -P.
                            Use `omnomnomics create-track-color-table` to build a custom palette in a writable folder.
                            Use `omnomnomics display-track-color-table` to visualize an existing color table.
    -P <DIR>:               Path to a folder with color tables.
                                    Default: packaged color_data_for_hubs directory
    -o:                     Overlay type (transparent|stacked|solid|none)
                                    Default: transparent
    -L:                     Email to use in trackhub
                                    Default: your@email.com

    -h:                     Print elaborate usage information
    -k: 		    Keep output of HISAT2 unpaired or not

By default, _Omnomnomics_ lets Snakemake reuse existing outputs that are already up to date. If you want to explicitly recompute the selected steps with the current settings, add `--rerun-selected-steps`.

The optional `--retention-policy` flag controls which large intermediate folders are retained after a successful run. `all` keeps the current behavior. `pruned` removes obvious bulk intermediates such as `trimmed_FASTQ` and `BAM` while retaining reusable downstream outputs such as `filtered_BAM`. `minimal` keeps only `FASTQ` plus the requested terminal output branches, for example `merged_hubs` and `DE_calling` for an RNA run that finishes at steps 9 and 11.

The optional `--max-project-size` flag adds a soft storage guard for space-heavy RNA outputs. When the configured cap would be exceeded, _Omnomnomics_ first attempts stage-safe cleanup of intermediates. Cleanup is guarded by step completion markers so that required upstream data are not deleted while dependent rules are still running. Before deleting `trimmed_FASTQ` or `BAM`, the pipeline caches lightweight flow-QC metric files (`*.trim_metrics.tsv` and mapper `*_stats.txt`) under `run_logs/flow_qc_cache` so downstream reporting can still read them. If the project would still exceed the cap, the pipeline logs a warning and skips BigWig or trackhub creation (steps 8 and 9) instead of failing on quota, while still allowing other requested branches such as read counting to continue.

Practical combinations:
- Full run with reusable outputs and bounded storage:
  `omnomnomics rna -i <EXPERIMENT_DIR> -g GRCh38 --retention-policy pruned --max-project-size 300G`
- Minimal final outputs for tight quota:
  `omnomnomics rna -i <EXPERIMENT_DIR> -g GRCh38 --retention-policy minimal --max-project-size 200G`

Utility commands:
- Monitor the latest run log interactively:
  `omnomnomics monitor -i <EXPERIMENT_DIR>`
- Launch the DE Shiny app for a local project copy:
  `omnomnomics de-app --project-dir /path/to/project_or_DE_calling`
- Search remote genome assemblies:
  `omnomnomics genomes list --species human`
- Show locally installed normalized assemblies:
  `omnomnomics genomes installed`
- Short aliases for the same genome helper actions:
  `omnomnomics genomes ls --species mouse`
  `omnomnomics genomes local --species mouse`
- Build a custom track color table:
  `omnomnomics create-track-color-table`
- Display an existing track color table:
  `omnomnomics display-track-color-table`

### DE quick notes
- RNA DE is public step 12 and requires RNA count-table input from step 11 plus a metadata table (`-m`).
- Use `--de-columns` and `--de-block` for automatic grouped design generation, or `--de-config <yaml>` for explicit control.
- `--de-config` can be repeated to run multiple DE analyses sequentially in one step-12 run.
- With repeated `--de-config`, each YAML must set a unique `io.out_dir`. Global `--de-out-dir` is rejected in that mode.
- Use `--de-out-dir <name>` to separate result trees when testing multiple DE settings.
- Step-12 QC outputs are shared at `DE_calling/qc/` and are not nested under `--de-out-dir`.
- Step 12 writes two scripts into `DE_calling/`:
  - `DE_analysis.rendered.R`: exact run reproduction with the resolved settings.
  - `DE_analysis.customization_guide.R`: same runnable script plus documented customization recipes.
- Explicit contrasts in `de_config` support:
  - Factor contrasts: `[factor, numerator, denominator]`
  - Coefficient contrasts: `{contrast_type: coefficient, coefficient_name: "...", label: "..."}`
- Available DESeq2 coefficient names for coefficient contrasts are exported to:
  - `results/qc/deseq2_results_names.tsv`
- Contrast-level QC tables are consolidated at top level:
  - `results/qc/contrast_testing_summary.tsv`
  - `results/differential_expression/contrast_summary.tsv`
  (Per-contrast folders keep only contrast-specific result files and plots.)
- DESeq2 runs serially inside each DE job by default. This avoids nested parallelism where Snakemake/Slurm already parallelizes across jobs, which can otherwise inflate memory use on large chromatin analyses. Advanced users can opt in with `deseq2.parallel: true` in a DE config when their scheduler and memory limits are appropriate.
- ATAC/ChIP post-DE signal heatmaps and profiles use existing per-sample BigWigs. By default, `--post-de-signal-policy auto` makes those BigWigs proper upstream dependencies of the post-DE analysis, so missing tracks are generated by step 8 as normal Snakemake jobs instead of being regenerated inside the post-DE worker. Use `require` for strict report-only reruns or `skip` to disable post-DE signal plotting during lightweight debugging runs. Signal heatmaps are rendered with deepTools `plotHeatmap`; profiles are rendered from the deepTools matrix with matplotlib so sample legends stay outside the signal axes. Equivalent post-DE signal sets are detected from the capped BED content, BigWig list, sample labels, and region label, then copied instead of recomputed. The summary table `DE_calling/analyze_peaks_de/signal/signal_runs.tsv` records `OK`, `REUSE`, and `SKIP` statuses plus any cap or reuse reason.
- ATAC/ChIP post-DE motif analysis uses MEME Suite tools. SEA is run for selected peak/region sets; AME is additionally run for top-ranked DE sets. Runs are bounded by workflow defaults (`post_de_motif_max_sets`, `post_de_motif_max_peaks`, `post_de_motif_window_bp`, `post_de_motif_timeout_seconds`, `post_de_motif_threads`, `post_de_motif_database`). When a DE-derived set exceeds the peak cap, regions are retained by adjusted P value and then absolute fold change rather than genomic input order; pre-DE unique sets preserve their source order. The default prioritizes promoter/promoter-region sets, then distal/top-ranked sets, and records skipped, failed, timed-out, or no-hit motif runs in `DE_calling/analyze_peaks_de/motifs/motif_runs.tsv` without failing the whole post-DE report. `NO_MOTIFS` means MEME completed but did not report enriched motifs passing its configured threshold. With `post_de_motif_database: auto`, omnomnomics uses a permanent MEME-format motif database cache under the configured genome assembly root and creates it on first use when possible. Run `omnomnomics genomes motifs` to prefetch it, or set `post_de_motif_database` or `OMNOMNOMICS_MEME_MOTIF_DATABASE` to pin a MEME-format motif database explicitly.

ATAC/ChIP post-DE status tables:
- `DE_calling/analyze_peaks_de/summary/set_manifest.tsv`: peak, region, gene-body, or bin sets created for downstream interpretation.
- `DE_calling/analyze_peaks_de/summary/analyze_peaks_de_report.tsv`: top-level report manifest for post-DE outputs.
- `DE_calling/analyze_peaks_de/signal/signal_runs.tsv`: per-set signal plotting status. `OK` means newly rendered, `REUSE` means an equivalent capped signal output was copied, and `SKIP` means plotting was intentionally skipped or the set was not eligible.
- `DE_calling/analyze_peaks_de/motifs/motif_runs.tsv`: per-set MEME status. `OK` means at least one report was produced, `NO_MOTIFS` means the motif tool completed with no enriched motifs, `SKIP` means the set was intentionally not tested, and `TIMEOUT` or `FAIL` are report-level motif failures.

Step 12 DE config run patterns:
- Single DE config:
  `omnomnomics rna -i <EXPERIMENT_DIR> -g GRCh38 -j 12 -m <metadata.tsv> --de-config de_grouped.yaml`
- Multiple DE configs in one sequential run:
  `omnomnomics rna -i <EXPERIMENT_DIR> -g GRCh38 -j 12 -m <metadata.tsv> --de-config de_grouped.yaml --de-config de_main_effects.yaml --de-config de_interaction.yaml`
- ATAC and ChIP differential chromatin follow the same metadata-driven CLI model, but ChIP broad-mark handling now depends strongly on `--broad-mode`.
- Multi-config mode requirements:
  - Every YAML must set `io.out_dir`.
  - `io.out_dir` values must be unique across the provided YAML files.
  - Do not use global `--de-out-dir` together with repeated `--de-config`.

### ChIP mode guide

Use the ChIP verb together with `--broad-mode` to select the feature definition that best matches the biology of the mark.

- `--broad-mode off`
  Use this for TF-like or otherwise punctate ChIP profiles. This is the narrow-peak path. It supports IDR-based consensus peak calling or MACS3-based parameter optimization.
- `--broad-mode domain`
  Use this for domain-like broad marks where MACS3 broad calling still makes biological sense. This path builds relaxed broad calls, replicate-aware concordance, pooled pseudoreplicate concordance, and then counts reads over the resulting consensus domains.
- `--broad-mode genebody`
  Use this for marks whose signal is best treated as enrichment across annotated gene bodies rather than called broad peaks. In this mode, Omnomnomics skips MACS3 peak discovery for DE features, derives gene-body features from the GTF, counts reads over those gene bodies, and runs chromatin DE on that feature matrix.
- `--broad-mode diffuse`
  Use this for very diffuse broad marks that are better modeled with a fixed-bin genome tiling than with peak calling. In this mode, Omnomnomics builds filtered genomic bins, counts reads over bins, runs chromatin DE on bins, and then merges adjacent significant bins into differential domains for downstream interpretation.

Practical mark classes:
- Typical narrow / TF-like ChIP: use `--broad-mode off`
- Broad domain-like histone marks: use `--broad-mode domain`
- Gene-body marks such as `H3K36me3` or `H3K79me2`: use `--broad-mode genebody`
- Very diffuse heterochromatin-style marks: use `--broad-mode diffuse`

Current ChIP gene-body behavior:
- The peak/feature stage creates annotation-derived gene-body BED features.
- The count-table stage counts reads over those gene bodies.
- Peak-QC annotation still provides genomic-region and nearest-gene metadata for those features.
- Pre-DE `analyze_peaks` produces gene-body metadata and all-gene-body signal summaries rather than peak intersections.
- Post-DE `analyze_peaks_de` still runs for gene-body mode and produces:
  - signal heatmaps and profiles for all gene bodies
  - signal heatmaps and profiles for DE gene-body sets
  - promoter-derived BED sets for motif analysis of DE genes, where available

Current ChIP diffuse behavior:
- The peak/feature stage derives fixed genomic bins from the BAM header chromosome sizes, filters them to standard chromosomes by default, excludes chrM by default, and subtracts the ENCODE blacklist when available.
- The count-table stage counts reads over those bins.
- Pre-DE `analyze_peaks` produces bin metadata plus global signal summaries for all bins.
- The chromatin DE stage runs on bins with diffuse-specific filtering defaults.
- Post-DE `analyze_peaks_de` produces DE bin sets and merged differential-domain BED and TSV summaries.
- Diffuse/bin FRiP is useful as a within-run consistency metric, but it should not be interpreted the same way as narrow-peak FRiP because broad tiled feature sets can cover much larger genome fractions.

Example ChIP entrypoints:

```bash
# Narrow / TF-like ChIP
omnomnomics chip -i <EXPERIMENT_DIR> -g GRCh38 --broad-mode off

# Broad domain-like histone marks
omnomnomics chip -i <EXPERIMENT_DIR> -g GRCh38 --broad-mode domain

# Gene-body marks
omnomnomics chip -i <EXPERIMENT_DIR> -g GRCh38 --broad-mode genebody

# Very diffuse broad marks
omnomnomics chip -i <EXPERIMENT_DIR> -g GRCh38 --broad-mode diffuse
```

Common step-12 config errors:
- `Global --de-out-dir cannot be combined with multiple --de-config files.`  
  Fix: remove global `--de-out-dir` and set `io.out_dir` in each YAML.
- `When multiple --de-config files are provided, each config must explicitly set io.out_dir.`  
  Fix: add `io.out_dir: <unique_name>` under `io:` in every YAML.
- `Multiple DE configs resolved to the same io.out_dir.`  
  Fix: assign unique `io.out_dir` values per YAML.
- `Explicit contrast dict contrast_type must be 'factor' or 'coefficient'.`  
  Fix: use `contrast_type: factor` for standard factor-level contrasts or `contrast_type: coefficient` for coefficient tests.

### Run the DE app locally (recommended workflow)

The DE app is intended for interactive local use after step 12 finished on HPC.

Typical workflow:
1. Run step 12 on HPC.
2. Copy or sync `DE_calling/` from HPC to your local machine.
3. Launch the app locally and load that project folder (or `DE_calling/` directly).

Minimum required local data:
- `DE_calling/metadata_derived.tsv`
- `DE_calling/qc/`
- at least one DE analysis output folder (for example `DE_calling/results/` or `DE_calling/<custom_out_dir>/`)

Launch options:

- Preferred (installed CLI):
  `omnomnomics de-app --project-dir /path/to/project_or_DE_calling`

- Optional flags:
  - `--port 3838`
  - `--host 127.0.0.1`
  - `--no-browser`

- If `omnomnomics` is not on your PATH but the source tree is available:
  `PYTHONPATH=src python -m omnomnomics.cli de-app --project-dir /path/to/project_or_DE_calling`

- Direct R fallback (advanced/manual):
  `R -e "shiny::runApp('src/omnomnomics/workflow/R/shiny_app')"`

Notes:
- Local and HPC environments can be different. You only need the step-12 outputs locally; raw FASTQ/BAM data are not required for app use.
- Running the app directly on HPC via SSH tunneling can work, but is cluster-specific and not the default supported workflow.

Step 12 YAML mini-templates:

```yaml
# de_grouped.yaml
io:
  out_dir: results_grouped

design:
  formula: "~ donor + replicate + de_group"

contrasts:
  mode: explicit
  explicit:
    items:
      - label: KD24_C_vs_KD24_L
        contrast_type: factor
        factor: de_group
        numerator: KD24_C
        denominator: KD24_L
      - label: NT_C_vs_NT_L
        contrast_type: factor
        factor: de_group
        numerator: NT_C
        denominator: NT_L
```

```yaml
# de_main_effects.yaml
io:
  out_dir: results_main_effects

design:
  formula: "~ donor + replicate + type + treatment"

contrasts:
  mode: explicit
  explicit:
    items:
      - label: KD24_vs_NT
        contrast_type: factor
        factor: type
        numerator: KD24
        denominator: NT
      - label: L_vs_C
        contrast_type: factor
        factor: treatment
        numerator: L
        denominator: C
```

```yaml
# de_interaction.yaml
io:
  out_dir: results_interaction

design:
  formula: "~ donor + replicate + type * treatment"

contrasts:
  mode: explicit
  explicit:
    items:
      - label: typeKD24.treatmentL
        contrast_type: coefficient
        coefficient_name: typeKD24.treatmentL
```

Tip: if the exact interaction coefficient name differs in your run, check `DE_calling/qc/deseq2_results_names.tsv` and update `coefficient_name` accordingly.

Additional utility command:

```bash
omnomnomics monitor -i <EXPERIMENT_DIR>
```

This watches the latest main run log in `<EXPERIMENT_DIR>/run_logs`, highlights step and sample-status lines, refreshes automatically, and exits on any key press.

## Logs
To ensure proper logging, multiple logs can be found. Inside the run_logs folder in your EXPERIMENT_DIR, a run log can be found created by _Omnomnomics_ which logs a lot of information about the current run and settings of the pipeline. In addition, slurm keeps a log of every submitted job which can be found inside the slurm_logs folder in your EXPERIMENT_DIR, and then inside its rule name folder. Here information about every submitted job can be found and potential errors while executing will be directed towards. Lastly, Snakemake also provides a log in the .snakemake folder. 

## Config files
For a run of _Omnomnomics_, three config layers are involved. The packaged workflow config contains pipeline defaults. The site config contains cluster-specific defaults such as partition and node characteristics. By default the CLI now looks for a user site config in `$XDG_CONFIG_HOME/omnomnomics/site.yaml` or `~/.config/omnomnomics/site.yaml` before falling back to the packaged site config. During a run, the CLI builds a run config that passes the resolved settings to Snakemake. This run config is written into the `run_configs` folder inside `EXPERIMENT_DIR`.

The packaged site-config template lives at [src/omnomnomics/workflow/config/site.yaml](/Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/config/site.yaml). Copy that file and edit the cluster-specific values for your own environment. Important site-level settings include the partition defaults, node layout, controller settings, worker constraints, and the Snakemake submission pacing controls `max_jobs`, `max_jobs_per_second`, and `max_status_checks_per_second`.

Set `worker_constraint` only when worker jobs must request a scheduler feature on your HPC. For example, on Snellius this can be set to `scratch-node` when jobs need reliable node-local scratch. Other clusters may use a different constraint name, or no worker constraint at all:

```yaml
worker_constraint: "scratch-node"
```

On clusters other than Snellius, the intended setup is therefore:

```bash
mkdir -p ~/.config/omnomnomics
cp /path/to/site.yaml ~/.config/omnomnomics/site.yaml
```

Use `--site-config` only when you want to override that default for a specific run.

## MultiQC
After completion of a run of _Omnomnomics_, MultiQC is used to parse and combine the results from the different bioinformatics tools. This helps to summarize the experiments that were run by giving a holistic view of all the samples and their results. In addition, it provides visual representation of the data to facilitate better interpretation. When running the full pipeline, MultiQC will generally produce a summarize html report, as well as data reports and plot reports. 

## Requirements:
- A working `omnomnomics` environment should be installed from `environment.yml`.
- The environment should provide the current workflow tools, including:
  - fastp
  - skewer
  - fastqc
  - hisat2
  - star
  - samtools
  - bedtools
  - deeptools
  - subread
  - macs3
  - idr, via `scripts/install_idr_helper.sh`, for default narrow ATAC/ChIP peak calling
  - run_spp.R, via `scripts/install_spp_helper.sh`, for ATAC/ChIP SPP cross-correlation QC
  - homer
  - multiqc
- On HPC systems, `sbatch` must be available because the Snakemake controller is submitted as a SLURM job.
- A valid sequence of steps should be included on the command line, if any. This should be 'valid' in the sense that the input/output of the included steps should 'fit' together. For every specified step, the input must already be present at the start OR it must be created by another specified steps. If intermediate steps are left out, Snakemake will attempt to run extra rules to create the necesarry inputs, but due to the complex naming schemes of input/output files, success is not guaranteed.
- Only one type of data should be present at the same time. This type must be specified on the command line. **Important**: Ensure that if intermediate inputs are given, that those are also not given on mixed type data. Example: ATAC & RNA data will yield filtered BAM files of name filtered_BAM/"sample_name".sorted.dups_marked.filtered.bam while CHIP data will yield filtered BAM files of name filtered_BAM/"sample_name".filtered.bam. If starting the pipeline from here, make sure to not give mixed time filtered BAM files.

## Troubleshooting

If a run exits immediately with a Snakemake `LockException`, first make sure no other Snakemake jobs are still active for the same experiment directory. If the lock is stale, try:

```bash
cd <EXPERIMENT_DIR>
snakemake --unlock
```

If that does not clear the problem, remove the stale Snakemake working directory manually and rerun:

```bash
cd <EXPERIMENT_DIR>
rm -r .snakemake
```

If a run fails with `Disk quota exceeded` or `No space left on device`, first stop remaining worker jobs, then clear the stale lock, then resume from a later step:

```bash
squeue -h -u "$USER" -n snakejob -o "%i" | xargs -r scancel
snakemake --unlock -s /path/to/Snakefile.smk --directory <EXPERIMENT_DIR> --profile /path/to/slurm_profile --configfile <EXPERIMENT_DIR>/run_configs/<run_config>.yaml
omnomnomics rna -i <EXPERIMENT_DIR> -g GRCh38 -j 5-11 --retention-policy pruned --max-project-size 300G
```

Adjust `-j` to match the latest successful stage in your run.

If optional HOMER tag directory export fails with messages such as `Could not find genome`, HOMER itself is installed but its genome data are missing. Install the required HOMER genome explicitly, for example:

```bash
configureHomer.pl -install hg38
configureHomer.pl -install mm39
```

On HPC systems where HOMER was installed inside a Conda or Mamba environment, `configureHomer.pl` is often not on `PATH` even though `makeTagDirectory` is. In that case, run it from the HOMER tree inside the active environment, for example:

```bash
cd "$HOME/conda/envs/omnomnomics/share/homer"
perl configureHomer.pl -install hg38
```

This can also be submitted as a batch job if you prefer not to run the download on a login node:

```bash
sbatch -p <PARTITION> -N 1 -n 1 -c 4 -t 01:00:00 --wrap='cd "$HOME/conda/envs/omnomnomics/share/homer" && perl configureHomer.pl -install hg38'
```

Adjust the environment path, partition, wall time, and genome name to match your setup.

For the optional HOMER export, `omnomnomics` translates common Ensembl/GenBank assembly names to HOMER aliases where needed:
- `GRCh38` and `GRCh38.p14` use HOMER genome `hg38`
- `GRCm39` uses HOMER genome `mm39`
