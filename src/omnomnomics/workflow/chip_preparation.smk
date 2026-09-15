# ChIP library preprocessing variants

## Omnomnomics Snake Rule ##

use rule run_fastp as run_fastp_se with:
    wildcard_constraints:
        sample=library_runtime.pattern(1, False)
    input:
        fastq1=lambda wildcards: library_runtime.raw(wildcards.sample, 1),
        fastq2=[]
    output:
        trimmed_fastq1=f"{library_manifest['folders']['trim']}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2=[],
        trim_metrics=f"{library_manifest['folders']['trim']}/{{sample}}.trim_metrics.tsv"

use rule run_skewer as run_skewer_se with:
    wildcard_constraints:
        sample=library_runtime.pattern(1, False)
    input:
        fastq1=lambda wildcards: library_runtime.raw(wildcards.sample, 1),
        fastq2=[]
    output:
        trimmed_fastq1=f"{library_manifest['folders']['trim']}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2=[],
        trim_metrics=f"{library_manifest['folders']['trim']}/{{sample}}.trim_metrics.tsv"

use rule run_fastqc as run_fastqc_se with:
    wildcard_constraints:
        sample=library_runtime.pattern(2, False)
    input:
        trimmed_fastq1=lambda wildcards: library_runtime.trimmed(wildcards.sample, 1),
        trimmed_fastq2=[]
    output:
        report1=f"{library_manifest['folders']['qc']}/{{sample}}.trimmed_fastqc.html",
        report2=[],
        report3=f"{library_manifest['folders']['qc']}/{{sample}}.trimmed_fastqc.zip",
        report4=[]

if config['THETRIMTOOL'] == "fastp":
    ruleorder: run_fastp_se > run_skewer_se
else:
    ruleorder: run_skewer_se > run_fastp_se

filtered_imports = library_runtime.filtered_imports()

rule import_filtered_bam:
    wildcard_constraints:
        sample="|".join(re.escape(sample) for sample in filtered_imports) or r"$.^"
    input:
        bam=lambda wildcards: filtered_imports[wildcards.sample]
    output:
        bam=f"{library_manifest['folders']['filtered']}/{{sample}}.filtered.bam"
    threads:
        Threads_Per_Rule['4']
    resources:
        mem_mb=Memory_Per_Rule['4'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['4']
    run:
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        shutil.copy2(input.bam, output.bam)
        log_it(logfile, f"Copied supplied filtered BAM {input.bam} to {output.bam}.", "INPUT LIBRARIES")
