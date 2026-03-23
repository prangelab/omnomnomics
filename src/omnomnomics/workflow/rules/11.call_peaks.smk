# Rule 11: Call Peaks

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import re
import glob
import csv
import shutil
import subprocess
import tempfile
from shlex import quote

def input_function(wildcards):
    input_folder1 = f"{experiment_dir}/{master_config['input_folders'][master_config['callpeaks_rule_num']-1][0]}"
    input_folder2 = f"{experiment_dir}/{master_config['input_folders'][master_config['callpeaks_rule_num']-1][1]}"
    input_files = []
    if config['THETYPE'] == "CHIP":
        for sample in samples2:
            input_files.append(f"{input_folder1}/{sample}.filtered.bam")
            input_files.append(f"{input_folder2}/{sample}.filtered.HOMER_tagDir.tar.gz")
    elif config['THETYPE'] == "ATAC":
        for sample in samples2:
            input_files.append(f"{input_folder1}/{sample}.sorted.dups_marked.filtered.bam")
    return input_files

rule call_peaks:
    input:
        input_function
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/extra_11.tmp"
    params:
        thetype= lambda wildcards: config['THETYPE'],  
        broad= lambda wildcards: config['BROAD'],
        input_sample= lambda wildcards: config['INPUT'],
        homer_input= lambda wildcards: config['HOMERINPUT'],
        genome= lambda wildcards: config['THEGENOME'],  
        experiment_dir= lambda wildcards: config['EXPERIMENT_DIR'], 
        name_fields= lambda wildcards: config['NAMEFIELDS'],  
        separator= lambda wildcards: config['THESEPARATOR'],
        thetype_field = lambda wildcards: config['THETYPEFIELD'] ,
        style= lambda wildcards: config['THESTYLE'],
        homersize= lambda wildcards: config['HOMERSIZE'],
        homermindist= lambda wildcards: config['HOMERMINDIST'],
        inputfolder1=lambda wildcards: f"{experiment_dir}/{master_config['input_folders'][master_config['callpeaks_rule_num']-1][0]}",
        inputfolder2= lambda wildcards: f"{experiment_dir}/{master_config['input_folders'][master_config['callpeaks_rule_num']-1][1]}",
        outputfolder= lambda wildcards: f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}"
    threads:
        lambda wildcards: Threads_Per_Rule['11']
    resources:
        mem_mb = lambda wildcards: Memory_Per_Rule['11'],
        partition = lambda wildcards: master_config['partition'],
        runtime = lambda wildcards: Runtime_Per_Rule['11']
    run:
        log_it(logfile, "Calling Peaks...", f"EXECUTING STEP {master_config['callpeaks_rule_num']}")
        log_it(logfile, f"Input folders: {params.inputfolder1} and {params.inputfolder2}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        def chip_style_label(style, broad):
            if broad == "1" or style == "histone":
                return "broad histone marks"
            return "TF / narrow peaks"

        def macs3_chip_presets(broad):
            presets = [
                ("q-0p05", "0.05"),
                ("q-0p01", "0.01"),
                ("q-0p001", "0.001"),
            ]
            if broad == "1":
                return [(f"{label}.broad", value) for label, value in presets]
            return presets

        def get_name_from_bam(inputfolder, group, name_fields, separator):
            group_escaped = group.replace(separator, '.*')
            cmd = f'basename "$(ls "{inputfolder}" | grep "{group_escaped}" - | grep ".bam$" | head -n1)" | cut -f{name_fields} -d{separator}'
            try:
                result = subprocess.check_output(cmd, shell=True).decode('utf-8').strip()
                return result
            except subprocess.CalledProcessError:
                return None
        

        def get_name_from_homer(inputfolder, group, name_fields, separator):
            group_escaped = group.replace(separator, '.*')
            cmd = f'basename "$(ls "{inputfolder}" | grep "{group_escaped}" - | grep ".HOMER_tagDir$" | head -n1)" | cut -f{name_fields} -d{separator}'
            try:
                result = subprocess.check_output(cmd, shell=True).decode('utf-8').strip()
                return result
            except subprocess.CalledProcessError:
                return None

        def ensure_peak_qc_dir(outputfolder):
            qc_dir = os.path.join(outputfolder, "peak_qc")
            os.makedirs(qc_dir, exist_ok=True)
            return qc_dir

        def log_peak_qc_versions():
            bedtools_version = subprocess.check_output(["bedtools", "--version"], stderr=subprocess.STDOUT)
            samtools_version = subprocess.check_output(["samtools", "--version"], stderr=subprocess.STDOUT).decode("utf-8").splitlines()[0]
            log_once(logfile, "step11.bedtools_version", "\n" + bedtools_version.decode("utf-8"), "BEDTOOLS VERSION")
            log_once(logfile, "step11.samtools_version", "\n" + samtools_version + "\n", "SAMTOOLS VERSION")
            if shutil.which("run_spp.R"):
                log_once(logfile, "step11.phantompeakqualtools", f"\nrun_spp.R: {shutil.which('run_spp.R')}\n", "PHANTOMPEAKQUALTOOLS")

        def macs3_peak_to_bed_path(peak_path):
            basename = os.path.basename(peak_path)
            if basename.endswith(".narrowPeak"):
                basename = basename[:-11]
            elif basename.endswith(".broadPeak"):
                basename = basename[:-10]
            return os.path.join(os.path.dirname(peak_path), f"{basename}.bed")

        def calculate_library_complexity_metrics(bam_path):
            with tempfile.TemporaryDirectory(prefix="omnomnomics_library_qc_") as tmpdir:
                counts_path = os.path.join(tmpdir, "fragment_counts.tsv")
                if config["PAIRED"]:
                    complexity_command = (
                        f"bedtools bamtobed -bedpe -i {quote(bam_path)} | "
                        "awk 'BEGIN{OFS=\"\\t\"} $1==$4 && $1!=\".\" {"
                        "start=($2<$5?$2:$5); end=($3>$6?$3:$6); print $1,start,end"
                        "}' | "
                        f"sort -T {quote(tmpdir)} -k1,1 -k2,2n -k3,3n | uniq -c > {quote(counts_path)}"
                    )
                else:
                    complexity_command = (
                        f"bedtools bamtobed -i {quote(bam_path)} | "
                        "awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$6}' | "
                        f"sort -T {quote(tmpdir)} -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c > {quote(counts_path)}"
                    )
                shell(complexity_command)

                total_reads = 0
                distinct_reads = 0
                one_read = 0
                two_reads = 0
                with open(counts_path) as handle:
                    for line in handle:
                        if not line.strip():
                            continue
                        count = int(line.strip().split()[0])
                        total_reads += count
                        distinct_reads += 1
                        if count == 1:
                            one_read += 1
                        elif count == 2:
                            two_reads += 1

            nrf = (distinct_reads / total_reads) if total_reads else 0.0
            pbc1 = (one_read / distinct_reads) if distinct_reads else 0.0
            pbc2 = (one_read / two_reads) if two_reads else ""
            return {
                "total_reads": total_reads,
                "distinct_reads": distinct_reads,
                "one_read_sites": one_read,
                "two_read_sites": two_reads,
                "nrf": nrf,
                "pbc1": pbc1,
                "pbc2": pbc2,
            }

        def calculate_cross_correlation_metrics(bam_path, output_dir, sample_name):
            crosscorr_prefix = os.path.join(output_dir, sample_name)
            crosscorr_table = f"{crosscorr_prefix}.cross_correlation.tsv"
            crosscorr_pdf = f"{crosscorr_prefix}.cross_correlation.pdf"
            if not shutil.which("run_spp.R"):
                return {
                    "est_frag_len": "",
                    "nsc": "",
                    "rsc": "",
                    "crosscorr_table": "",
                    "crosscorr_pdf": "",
                }

            run_spp_command = (
                f"run_spp.R -c={quote(bam_path)} -savp={quote(crosscorr_pdf)} "
                f"-out={quote(crosscorr_table)} -p={threads}"
            )
            log_it(logfile, run_spp_command, "PHANTOMPEAKQUALTOOLS COMMAND")
            shell(run_spp_command)

            metrics_line = ""
            with open(crosscorr_table) as handle:
                for line in handle:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        metrics_line = line
                        break

            parts = [part.strip() for part in metrics_line.split("\t")] if "\t" in metrics_line else [part.strip() for part in metrics_line.split(",")]
            est_frag_len = parts[2] if len(parts) > 2 else ""
            nsc = parts[8] if len(parts) > 8 else ""
            rsc = parts[9] if len(parts) > 9 else ""
            return {
                "est_frag_len": est_frag_len,
                "nsc": nsc,
                "rsc": rsc,
                "crosscorr_table": crosscorr_table,
                "crosscorr_pdf": crosscorr_pdf,
            }

        def calculate_peak_qc_metrics(peak_bed, bam_files):
            with tempfile.TemporaryDirectory(prefix="omnomnomics_peak_qc_") as tmpdir:
                merged_bed = os.path.join(tmpdir, "merged_peaks.bed")
                multicov_output = os.path.join(tmpdir, "multicov.tsv")
                shell(
                    f"sort -k1,1 -k2,2n -k3,3n {quote(peak_bed)} | "
                    f"bedtools merge -i - > {quote(merged_bed)}"
                )
                shell(
                    f"bedtools multicov -bed {quote(merged_bed)} -bams {' '.join(quote(path) for path in bam_files)} "
                    f"> {quote(multicov_output)}"
                )

                peak_count = 0
                total_peak_bp = 0
                reads_in_peaks = 0
                with open(multicov_output, newline="") as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    for row in reader:
                        peak_count += 1
                        total_peak_bp += int(row[2]) - int(row[1])
                        reads_in_peaks += sum(int(value) for value in row[3:])

                total_aligned_reads = 0
                for bam_path in bam_files:
                    total_aligned_reads += int(
                        subprocess.check_output(
                            ["samtools", "view", "-c", "-F", "260", bam_path],
                            stderr=subprocess.STDOUT,
                        ).decode("utf-8").strip()
                    )

            frip = (reads_in_peaks / total_aligned_reads) if total_aligned_reads else 0.0
            return {
                "peak_count": peak_count,
                "total_peak_bp": total_peak_bp,
                "reads_in_peaks": reads_in_peaks,
                "total_aligned_reads": total_aligned_reads,
                "frip": frip,
            }

        def write_sample_qc_outputs(outputfolder, thetype, qc_rows):
            if not qc_rows:
                return

            qc_dir = ensure_peak_qc_dir(outputfolder)
            qc_table = os.path.join(qc_dir, f"{thetype.lower()}.sample_qc_metrics.tsv")
            with open(qc_table, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "assay",
                    "sample",
                    "bam_file",
                    "total_reads",
                    "distinct_reads",
                    "one_read_sites",
                    "two_read_sites",
                    "nrf",
                    "pbc1",
                    "pbc2",
                    "est_frag_len",
                    "nsc",
                    "rsc",
                    "crosscorr_table",
                    "crosscorr_pdf",
                ])
                for row in qc_rows:
                    writer.writerow([
                        row["assay"],
                        row["sample"],
                        row["bam_file"],
                        row["total_reads"],
                        row["distinct_reads"],
                        row["one_read_sites"],
                        row["two_read_sites"],
                        f"{row['nrf']:.6f}",
                        f"{row['pbc1']:.6f}",
                        row["pbc2"] if row["pbc2"] == "" else f"{float(row['pbc2']):.6f}",
                        row["est_frag_len"],
                        row["nsc"],
                        row["rsc"],
                        row["crosscorr_table"],
                        row["crosscorr_pdf"],
                    ])

            return qc_table

        def write_peak_qc_outputs(outputfolder, thetype, qc_rows):
            if not qc_rows:
                return

            qc_dir = ensure_peak_qc_dir(outputfolder)
            qc_table = os.path.join(qc_dir, f"{thetype.lower()}.peak_qc_metrics.tsv")
            with open(qc_table, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "assay",
                    "group",
                    "peak_set",
                    "peak_file",
                    "bam_count",
                    "peak_count",
                    "total_peak_bp",
                    "reads_in_peaks",
                    "total_aligned_reads",
                    "frip",
                ])
                for row in qc_rows:
                    writer.writerow([
                        row["assay"],
                        row["group"],
                        row["peak_set"],
                        row["peak_file"],
                        row["bam_count"],
                        row["peak_count"],
                        row["total_peak_bp"],
                        row["reads_in_peaks"],
                        row["total_aligned_reads"],
                        f"{row['frip']:.6f}",
                    ])

            return qc_table

        def write_peak_qc_summary_pdf(outputfolder, thetype, sample_rows, peak_rows):
            if not sample_rows and not peak_rows:
                return

            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages

            log_once(logfile, "step11.matplotlib_version", f"\nmatplotlib {matplotlib.__version__}\n", "MATPLOTLIB VERSION")

            qc_dir = ensure_peak_qc_dir(outputfolder)
            pdf_path = os.path.join(qc_dir, f"{thetype.lower()}.peak_qc_summary.pdf")

            plt.style.use("seaborn-v0_8-whitegrid")
            with PdfPages(pdf_path) as pdf:
                if peak_rows:
                    peak_labels = [row["peak_set"] for row in peak_rows]
                    frip_values = [float(row["frip"]) for row in peak_rows]
                    peak_counts = [float(row["peak_count"]) for row in peak_rows]
                    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
                    axes[0].bar(peak_labels, frip_values, color="#4C78A8")
                    axes[0].set_title(f"{thetype} peak QC: FRiP")
                    axes[0].set_ylabel("FRiP")
                    axes[0].tick_params(axis="x", rotation=45, labelsize=9)
                    axes[1].bar(peak_labels, peak_counts, color="#D95F02")
                    axes[1].set_title(f"{thetype} peak QC: peak count")
                    axes[1].set_ylabel("Peaks")
                    axes[1].tick_params(axis="x", rotation=45, labelsize=9)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

                if sample_rows:
                    sample_labels = [row["sample"] for row in sample_rows]
                    nrf_values = [float(row["nrf"]) for row in sample_rows]
                    pbc1_values = [float(row["pbc1"]) for row in sample_rows]
                    pbc2_values = [float(row["pbc2"]) if row["pbc2"] != "" else 0.0 for row in sample_rows]
                    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
                    axes[0].bar(sample_labels, nrf_values, color="#2A9D8F")
                    axes[0].set_title(f"{thetype} library complexity: NRF")
                    axes[0].set_ylabel("NRF")
                    axes[0].tick_params(axis="x", rotation=45, labelsize=9)
                    axes[1].bar(sample_labels, pbc1_values, color="#E9C46A")
                    axes[1].set_title(f"{thetype} library complexity: PBC1")
                    axes[1].set_ylabel("PBC1")
                    axes[1].tick_params(axis="x", rotation=45, labelsize=9)
                    axes[2].bar(sample_labels, pbc2_values, color="#A44A3F")
                    axes[2].set_title(f"{thetype} library complexity: PBC2")
                    axes[2].set_ylabel("PBC2")
                    axes[2].tick_params(axis="x", rotation=45, labelsize=9)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

                    nsc_values = [float(row["nsc"]) for row in sample_rows if row["nsc"] not in ("", "NA")]
                    rsc_values = [float(row["rsc"]) for row in sample_rows if row["rsc"] not in ("", "NA")]
                    if nsc_values and rsc_values:
                        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
                        axes[0].bar(sample_labels, [float(row["nsc"]) if row["nsc"] not in ("", "NA") else 0.0 for row in sample_rows], color="#6A4C93")
                        axes[0].set_title(f"{thetype} enrichment: NSC")
                        axes[0].set_ylabel("NSC")
                        axes[0].tick_params(axis="x", rotation=45, labelsize=9)
                        axes[1].bar(sample_labels, [float(row["rsc"]) if row["rsc"] not in ("", "NA") else 0.0 for row in sample_rows], color="#F28482")
                        axes[1].set_title(f"{thetype} enrichment: RSC")
                        axes[1].set_ylabel("RSC")
                        axes[1].tick_params(axis="x", rotation=45, labelsize=9)
                        fig.tight_layout()
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)

            log_it(logfile, f"Peak QC summary PDF: {pdf_path}")

        def collect_sample_qc_rows(inputfolder, thetype):
            sample_rows = []
            qc_dir = ensure_peak_qc_dir(params.outputfolder)
            bam_suffix = ".sorted.dups_marked.filtered.bam" if thetype == "ATAC" else ".filtered.bam"
            for sample in samples2:
                bam_path = os.path.join(inputfolder, f"{sample}{bam_suffix}")
                if not os.path.exists(bam_path):
                    continue
                log_it(logfile, f"Calculating library complexity metrics for {sample}...")
                complexity_metrics = calculate_library_complexity_metrics(bam_path)
                log_it(logfile, f"Calculating cross-correlation metrics for {sample}...")
                crosscorr_metrics = calculate_cross_correlation_metrics(bam_path, qc_dir, sample)
                sample_rows.append({
                    "assay": thetype,
                    "sample": sample,
                    "bam_file": os.path.basename(bam_path),
                    **complexity_metrics,
                    **crosscorr_metrics,
                })
            return sample_rows

        def call_peaks(logfile, thetype, inputfolder1, inputfolder2, outputfolder, separator, thetype_field, name_fields, broad, input_sample, homer_input, homer_size, homer_mindist, thestyle):
            if thetype == "RNA":
                log_it(logfile, "Not a ChIP- or ATAC-seq experiment, skipping this step...")
                return
            
            log_it(logfile, "Calling peaks...", f"EXECUTING STEP {master_config['callpeaks_rule_num']}") 
            log_it(logfile, f"The type = {thetype}")
            log_peak_qc_versions()

            # Report version
            macs3_version = subprocess.check_output(["macs3", "--version"], stderr=subprocess.STDOUT)
            log_it(logfile, "\n"+macs3_version.decode("utf-8"), "MACS3 VERSION")

            if thetype == "CHIP":
                chip_qc_rows = []
                chip_sample_qc_rows = []
                #If ChIP, call peaks
                log_it(logfile, f"Finding ChIP enriched regions for {chip_style_label(thestyle, broad)}...")
                if broad == "1" or thestyle == "histone":
                    log_it(logfile, "NSC and RSC are reported for broad histone marks, but they are typically less informative there than for TF / narrow peaks.")

                log_it(logfile, "Calling peaks with MACS3...")
                
                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder1,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1][0]) 
                
                # Distribute samples over groups based on the file name pattern
                bams = [os.path.basename(bam).split(separator)[thetype_field - 1] for bam in glob.glob(f"{inputfolder1}/*.bam")]
                chip_groups = sorted(set(bams))
                
                for group in chip_groups:
                    # BAMS=( $(ls "$THEINFOLDER" | grep $(echo $THEGROUP | sed "s|$THESEPARATOR|.*|g") - | grep ".bam$") )
                    bams = [bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)] 
                    bams = [os.path.join(inputfolder1, bam) for bam in bams]
                    
                    log_it(logfile, f"Calling peaks for: {group}...")
                    log_it(logfile, f"Files in group: {', '.join(bams)}")

                    bam_basename = os.path.basename(bams[0])
                    fields = bam_basename.split(separator)

                    # Function to parse the NAMEFIELDS and return the specified field indices
                    def parse_namefields(namefields):
                        indices = []
                        for part in namefields.split(','):
                            if '-' in part:
                                start, end = map(int, part.split('-'))
                                indices.extend(range(start, end + 1))
                            else:
                                indices.append(int(part))
                        return sorted(set(indices))  # Remove duplicates and sort the indices
                    field_indices = parse_namefields(name_fields)

                    the_name = separator.join(fields[i-1] for i in field_indices if i-1 < len(fields))

                    if input_sample == "NA":#No input sample provided
                        if broad == "1": # Check if we should call broad peaks
                            log_it(logfile, "Calling broad histone marks with MACS3 using q values 0.05, 0.01, and 0.001...")
                            for ext, qvalue in macs3_chip_presets(broad):
                                shell(f"""macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -q {qvalue} --broad --broad-cutoff 0.1 --verbose 0 """)
                        else:
                            log_it(logfile, "Calling TF / narrow peaks with MACS3 using q values 0.05, 0.01, and 0.001...")
                            for ext, qvalue in macs3_chip_presets(broad):
                                shell(f"""macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -q {qvalue} --verbose 0""")
                    else:# We have input!
                        if broad == "1": # Check if we should call broad peaks
                            log_it(logfile, "Calling broad histone marks with MACS3 using q values 0.05, 0.01, and 0.001...")
                            for ext, qvalue in macs3_chip_presets(broad):
                                shell(f"""macs3 callpeak -t {' '.join(bams)} -c {input_sample} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -q {qvalue} --broad --broad-cutoff 0.1 --verbose 0""")
                        else:
                            log_it(logfile, "Calling TF / narrow peaks with MACS3 using q values 0.05, 0.01, and 0.001...")
                            for ext, qvalue in macs3_chip_presets(broad):
                                shell(f"""macs3 callpeak -t {' '.join(bams)} -c {input_sample} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -q {qvalue} --verbose 0""")
                # Clean up the output
                log_it(logfile, "Cleaning up MACS3 output (keeping only narrowPeak and broadPeak files)...")
                for group in chip_groups: 
                    the_name = get_name_from_bam(inputfolder1, group, name_fields, separator)
                    log_it(logfile, f"find {outputfolder} -type f -name {the_name}.MACS3* ! \( -name *broadPeak -o -name *narrowPeak \) -delete")
                    shell(f"find {outputfolder} -type f -name {the_name}.MACS3* ! \( -name *broadPeak -o -name *narrowPeak \) -delete ") # Delete all redundant MACS3 output

                log_it(logfile, "Converting to a clean 3 column BED format...")
                for file in glob.glob(f"{outputfolder}/*{'narrowPeak' if broad != '1' else 'broadPeak'}"):
                    sorted_bed = macs3_peak_to_bed_path(file)
                    log_it(logfile, f"cut -f1-3 {file} | sort -k1,1 -k2,2n -k3,3n > {sorted_bed}" )
                    shell(f"cut -f1-3 {file} | sort -k1,1 -k2,2n -k3,3n > {sorted_bed}")
                    os.remove(file)

                log_it(logfile, "Calculating peak QC metrics for MACS3 peak sets...")
                for group in chip_groups:
                    bams = [bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)]
                    bams = [os.path.join(inputfolder1, bam) for bam in bams]
                    the_name = get_name_from_bam(inputfolder1, group, name_fields, separator)
                    for peak_bed in sorted(glob.glob(f"{outputfolder}/{the_name}.MACS3*.bed")):
                        metrics = calculate_peak_qc_metrics(peak_bed, bams)
                        chip_qc_rows.append({
                            "assay": "CHIP",
                            "group": group,
                            "peak_set": os.path.basename(peak_bed).replace(".bed", ""),
                            "peak_file": os.path.basename(peak_bed),
                            "bam_count": len(bams),
                            **metrics,
                        })
                chip_sample_qc_rows = collect_sample_qc_rows(inputfolder1, thetype)
                peak_qc_table = write_peak_qc_outputs(outputfolder, thetype, chip_qc_rows)
                sample_qc_table = write_sample_qc_outputs(outputfolder, thetype, chip_sample_qc_rows)
                if peak_qc_table:
                    log_it(logfile, f"Peak QC metrics: {peak_qc_table}")
                if sample_qc_table:
                    log_it(logfile, f"Sample QC metrics: {sample_qc_table}")
                write_peak_qc_summary_pdf(outputfolder, thetype, chip_sample_qc_rows, chip_qc_rows)
                
                # Call peaks using HOMER
                log_it(logfile, "Calling peaks with HOMER...")
                log_it(logfile, f"Input folder: {inputfolder2}")
                
                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder2,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1][1])
                
                # See if we need to unpack the tag dirs
                if glob.glob(f"{inputfolder2}/*tagDir.tar.gz"):
                    for tagdir in glob.glob(f"{inputfolder2}/*tagDir.tar.gz"):
                        tagdir_basename = os.path.basename(tagdir)
                        log_it(logfile, f"cd {inputfolder2} && tar --strip-components=1 -xzf {tagdir_basename}")
                        shell(f"cd {inputfolder2} && tar --strip-components=1 -xzf {tagdir_basename}")
                # Fetch tagdirs
                tagdirs = glob.glob(f"{inputfolder2}/*.HOMER_tagDir")
                
                if homer_input == "NA": # No input sample was provided
                    if broad == "1": # Check if we should call broad peaks
                        log_it(logfile, "Calling broad histone marks with HOMER...") 
                        for tagdir in tagdirs:
                            log_it(logfile, f"findPeaks {tagdir} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                            shell(f"findPeaks {tagdir} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                    else:
                        log_it(logfile, "Calling peaks with HOMER...")
                        for tagdir in tagdirs:
                            log_it(logfile, f"findPeaks {tagdir} -style {thestyle} -o auto")
                            shell(f"findPeaks {tagdir} -style {thestyle} -o auto")
                else:  # We have input!
                    if broad == "1": # Check if we should call broad peaks
                        log_it(logfile, "Calling broad histone marks with HOMER...")
                        for tagdir in tagdirs:
                            log_it(logfile, f"findPeaks {tagdir} -i {homer_input} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                            shell(f"findPeaks {tagdir} -i {homer_input} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                    else:
                        log_it(logfile, "Calling peaks with HOMER...")
                        for tagdir in tagdirs:
                            log_it(logfile, f"findPeaks {tagdir} -i {homer_input} -style {thestyle} -o auto")
                            shell(f"findPeaks {tagdir} -i {homer_input} -style {thestyle} -o auto")
                
                # Convert peaks to BED format and clean up unwanted contigs (chrUn | alt | random)
                log_it(logfile, "Converting peaks to BED format and cleaning up unwanted contigs (chrUn | alt | random)...")
                for tagdir in tagdirs:
                    log_it(logfile, f"Testpp {tagdir}")
                    if thestyle == "factor":
                        shell(f"pos2bed.pl {tagdir}/peaks.txt | grep -v '#\\|alt\\|Un\\|random' | sort -k1,1 -k2,2n -k3,3n | cut -f-3 > {tagdir}/peaks.bed ")
                    elif thestyle == "histone":
                        shell(f"pos2bed.pl {tagdir}/regions.txt | grep -v '#\\|alt\\|Un\\|random' | sort -k1,1 -k2,2n -k3,3n | cut -f-3 > {tagdir}/regions.bed ")

                ## Distribute samples over groups based on the file name pattern
                # Build group list
                tagdirs2 = [os.path.basename(tagdir).split(separator)[thetype_field - 1] for tagdir in tagdirs]
                chip_groups = sorted(set(tagdirs2))

                log_it(logfile, f"chip_groups = {chip_groups}")

                # Iterate over the groups to merge the BED files
                for group in chip_groups:
                    # Fetch samples
                    beds = [os.path.join(inputfolder2, bed) for bed in os.listdir(inputfolder2) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.HOMER_tagDir$", bed)]
                    if thestyle == "factor":
                        beds = [f"{bed}/peaks.bed" for bed in beds]
                    elif thestyle == "histone":
                        beds = [f"{bed}/regions.bed" for bed in beds]
                    log_it(logfile, f"Merging peaks for: {group}...")
                    log_it(logfile, f"Samples in group: {', '.join(beds)}")
                    
                    the_name = get_name_from_homer(inputfolder2, group, name_fields, separator)
                    sorted_bed = f"{outputfolder}/{the_name}.HOMER.merged_peaks.bed"
                    log_it(logfile, f"module load bedtools && cat {' '.join(beds)} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {sorted_bed}")
                    shell(f"""module load bedtools && cat {' '.join(beds)} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {sorted_bed}""")
            else:
                atac_qc_rows = []
                atac_sample_qc_rows = []
                # If ATAC, call open regions
                log_it(logfile, "Calling ATAC open chromatin peaks...")

                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder1,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1][0]) 

                # Distribute samples over groups based on the file name pattern
                bams = [os.path.basename(bam).split(separator)[thetype_field - 1] for bam in glob.glob(f"{inputfolder1}/*.bam")]
                atac_groups = sorted(set(bams))

                for group in atac_groups: # Iterate over the groups
                    bams = [bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)]
                    bams = [os.path.join(inputfolder1, bam) for bam in bams] # Fetch BAM files
                    
                    log_it(logfile, f"Calling ATAC open chromatin peaks for: {group}...")
                    log_it(logfile, f"Files in group: {', '.join(bams)}")
                    # Set group name
                    the_name = get_name_from_bam(inputfolder1, group, name_fields, separator)

                    # Call open regions with MACS3 hmmr
                    # Turns out this MACS subroutine is till buggy and not functining prooerply, check back with it later but for now just call as ChIP peaks
                    #logIt "Calling open regions with MACS3 hmmratac..."
                    #logIt "	macs3 hmmratac -b $BAMS --outdir \"$THEOUTFOLDER\" -n $NAMEFIELDS --verbose 0 &"
                    #macs3 hmmratac -b ${BAMS[*]} --outdir "$THEOUTFOLDER" -n $THENAME --verbose 2 &

                    log_it(logfile, "Calling ATAC open chromatin peaks with MACS3, q value 0.01...")
                    log_it(logfile, f"macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.q-0p01 -q 0.01 --verbose 0" )
                    shell(f"""macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.q-0p01 -q 0.01 --verbose 0 """)

                log_it(logfile, "Cleaning up MACS3 output (keeping only narrowPeak files)...")
                for group in atac_groups: ## Clean up the output
                    # Set group name
                    the_name = get_name_from_bam(inputfolder1, group, name_fields, separator)
                    shell(f"find {outputfolder} -type f -name {the_name}.MACS3* ! -name *narrowPeak -delete")
                
                # Create clean BED files
                log_it(logfile, "Converting to a clean 3 column BED format...")
                for file in glob.glob(f"{outputfolder}/*narrowPeak"):
                    sorted_bed = macs3_peak_to_bed_path(file)
                    shell(f"cut -f1-3 {file} | sort -k1,1 -k2,2n -k3,3n > {sorted_bed}")
                    os.remove(file)

                # Concatenate the peak files
                shell(f"""cat {outputfolder}/*.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {outputfolder}/all_groups.merged_peaks.bed """)
                log_it(logfile, f"Total merged peaks: {subprocess.getoutput(f'wc -l {outputfolder}/all_groups.merged_peaks.bed').strip().split()[0]}")

                log_it(logfile, "Calculating peak QC metrics for MACS3 peak sets...")
                for group in atac_groups:
                    bams = [bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)]
                    bams = [os.path.join(inputfolder1, bam) for bam in bams]
                    the_name = get_name_from_bam(inputfolder1, group, name_fields, separator)
                    peak_bed = os.path.join(outputfolder, f"{the_name}.MACS3.q-0p01_peaks.bed")
                    if os.path.exists(peak_bed):
                        metrics = calculate_peak_qc_metrics(peak_bed, bams)
                        atac_qc_rows.append({
                            "assay": "ATAC",
                            "group": group,
                            "peak_set": os.path.basename(peak_bed).replace(".bed", ""),
                            "peak_file": os.path.basename(peak_bed),
                            "bam_count": len(bams),
                            **metrics,
                        })
                atac_sample_qc_rows = collect_sample_qc_rows(inputfolder1, thetype)
                peak_qc_table = write_peak_qc_outputs(outputfolder, thetype, atac_qc_rows)
                sample_qc_table = write_sample_qc_outputs(outputfolder, thetype, atac_sample_qc_rows)
                if peak_qc_table:
                    log_it(logfile, f"Peak QC metrics: {peak_qc_table}")
                if sample_qc_table:
                    log_it(logfile, f"Sample QC metrics: {sample_qc_table}")
                write_peak_qc_summary_pdf(outputfolder, thetype, atac_sample_qc_rows, atac_qc_rows)

        call_peaks(
            logfile,
            thetype=params.thetype,
            inputfolder1=params.inputfolder1,
            inputfolder2=params.inputfolder2,
            outputfolder=params.outputfolder,
            separator=params.separator,
            thetype_field=int(params.thetype_field),
            name_fields=params.name_fields,  
            broad=params.broad,
            input_sample=params.input_sample,
            homer_input=params.homer_input,
            homer_size=int(params.homersize),
            homer_mindist=int(params.homermindist),
            thestyle=params.style
        )
        shell(f"""echo "necessity file for callpeaks. can delete this." > {params.outputfolder}/extra_11.tmp""")
