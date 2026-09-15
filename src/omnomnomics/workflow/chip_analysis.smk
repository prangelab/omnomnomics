# ChIP analytical inputs and testing regions

## Omnomnomics Snake Rule ##
from omnomnomics.chip_analysis import (
    ControlSets, annotate_regions, bam_fragments, bedgraph_bigwig, count_libraries,
    enrichment_rows, macs_command, read_regions, region_id, region_specs_match, sha256_file,
    support_annotations, write_regions, write_tsv,
)
from omnomnomics.library_manifest import reference_lengths, inspect_bam
from shlex import quote
from pathlib import Path
import csv
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile

chip_analysis_enabled = config.get("THETYPE") == "CHIP"
chip_control_sets = ControlSets(library_manifest, config.get("INPUT", "NA"), samples2)
chip_peak_dir = f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num'] - 1]}"
chip_count_dir = f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num'] - 1]}"
chip_filtered_dir = f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num'] - 1]}"
chip_region_bed = f"{chip_peak_dir}/chip_testing/regions.bed"
chip_region_manifest = f"{chip_peak_dir}/chip_testing/regions.tsv"
chip_testing_metadata = f"{chip_count_dir}/chip_testing_metadata.tsv"
chip_count_selection_spec = f"{experiment_dir}/run_configs/chip_count_selection.json"
chip_custom_bed = str(config.get("CHIP_REGIONS_BED", "NA"))
chip_analytical_samples = list(samples2)
chip_fragment_samples = sorted(set(chip_analytical_samples) | set(chip_control_sets.ids(chip_analytical_samples)))
chip_parent_aliases = {}
chip_se_extension = int(config.get("CHIP_SE_FRAGMENT_LENGTH", 200))

def chip_bam(sample):
    if sample == chip_control_sets.legacy_id and chip_control_sets.legacy:
        return chip_control_sets.legacy
    if library_runtime and sample in library_runtime.libraries:
        return library_runtime.library(sample)["filtered_bam"]
    return f"{chip_filtered_dir}/{sample}.filtered.bam"

def chip_layout(sample):
    if library_runtime and sample in library_runtime.libraries:
        return library_runtime.paired(sample)
    return inspect_bam(chip_bam(sample), reference_lengths(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"]))

def chip_fragment_path(sample):
    return f"{chip_peak_dir}/chip_fragments/{sample}.bedpe"

def chip_control_dependencies(samples):
    if not chip_analysis_enabled:
        return []
    return [chip_controls_spec, *[chip_fragment_path(sample) for sample in chip_control_sets.ids(samples)]]

def chip_peak_dependencies(samples):
    return [chip_peak_calls_spec, *chip_control_dependencies(samples)] if chip_analysis_enabled else []

def chip_parents(bam_paths):
    parents = set()
    by_path = {os.path.abspath(chip_bam(sample)): sample for sample in chip_analytical_samples}
    for path in bam_paths:
        identity = os.path.abspath(str(path))
        if identity in chip_parent_aliases:
            parents.update(chip_parent_aliases[identity])
        elif identity in by_path:
            parents.add(by_path[identity])
        else:
            raise ValueError(f"Unresolved ChIP parent for peak calling: {path}")
    return sorted(parents)

def chip_register_alias(destination, source_bams):
    if chip_analysis_enabled:
        chip_parent_aliases[os.path.abspath(str(destination))] = chip_parents(source_bams)

def chip_macs_inputs(treatment_bams):
    parents = chip_parents(treatment_bams)
    converted = []
    for path in treatment_bams:
        identity = os.path.abspath(str(path))
        direct = next((sample for sample in parents if os.path.abspath(chip_bam(sample)) == identity), None)
        if direct:
            converted.append(chip_fragment_path(direct))
        else:
            destination = str(path) + ".bedpe"
            layouts = {chip_layout(sample) for sample in chip_parent_aliases[identity]}
            layout = next(iter(layouts)) if len(layouts) == 1 else None
            stat = os.stat(path)
            state = json.dumps({"bam": identity, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns, "paired": layout, "se_extension": chip_se_extension}, sort_keys=True)
            state_path = destination + ".source.json"
            if not os.path.isfile(destination) or not os.path.isfile(state_path) or Path(state_path).read_text() != state:
                bam_fragments(path, destination, layout, chip_se_extension)
                write_stable_text(state_path, state)
            converted.append(destination)
    return converted, [chip_fragment_path(sample) for sample in chip_control_sets.ids(parents)]

if chip_analysis_enabled:
    chip_genome_size = config.get("CHIP_EFFECTIVE_GENOME_SIZE") or sum(reference_lengths(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"]).values())
    control_spec = {
        "macs_genome_size": chip_genome_size,
        "samples": sorted(chip_analytical_samples),
        "associations": (library_manifest or {}).get("associations", []),
        "legacy_input": chip_control_sets.legacy,
        "libraries": [{"id": sample, "bam": chip_bam(sample), "layout": library_runtime.paired(sample) if library_runtime and sample in library_runtime.libraries else "inspect_bam"} for sample in chip_fragment_samples],
        "macs_fragment_representation": "BEDPE; PE observed insert; SE inferred extension",
        "se_fragment_length": chip_se_extension,
    }
    content = json.dumps(control_spec, sort_keys=True, indent=2) + "\n"
    chip_controls_spec = f"{experiment_dir}/run_configs/chip_controls.{hashlib.sha256(content.encode()).hexdigest()[:16]}.json"
    region_spec = {"controls": control_spec, "genome": config["THEGENOME"], "mode": config.get("BROAD_MODE", "off"), "q": config.get("CHIP_JOINT_PEAK_Q", 0.1), "bed": chip_custom_bed, "bed_sha256": sha256_file(chip_custom_bed) if chip_custom_bed != "NA" else None}
    region_content = json.dumps(region_spec, sort_keys=True, indent=2) + "\n"
    chip_regions_spec = f"{experiment_dir}/run_configs/chip_regions.{hashlib.sha256(region_content.encode()).hexdigest()[:16]}.json"
    caller_keys = [key for key in config if key.startswith("IDR_") or key.startswith("CHIP_BROAD_")]
    caller_keys.extend(["BROAD_MODE", "NARROW_PEAK_STRATEGY", "CHIP_PEAK_OPT_MODE"])
    caller_spec = {"controls": control_spec, "groups": {sample: sample_type_for_sample(sample) for sample in sorted(chip_analytical_samples)}, "settings": {key: config.get(key) for key in sorted(set(caller_keys))}}
    caller_content = json.dumps(caller_spec, sort_keys=True, indent=2) + "\n"
    chip_peak_calls_spec = f"{experiment_dir}/run_configs/chip_peak_calls.{hashlib.sha256(caller_content.encode()).hexdigest()[:16]}.json"
    count_content = json.dumps({"regions": chip_regions_spec, "enrichment_min_input": config.get("CHIP_ENRICHMENT_MIN_INPUT", 5), "spp_gate": config.get("SPP_GATE", "warn")}, sort_keys=True, indent=2) + "\n"
    chip_count_settings_spec = f"{experiment_dir}/run_configs/chip_count_settings.{hashlib.sha256(count_content.encode()).hexdigest()[:16]}.json"
    if not is_worker_job:
        write_stable_text(chip_controls_spec, content)
        write_stable_text(chip_regions_spec, region_content)
        write_stable_text(chip_peak_calls_spec, caller_content)
        write_stable_text(chip_count_settings_spec, count_content)

rule chip_fragments:
    wildcard_constraints:
        sample="|".join(re.escape(sample) for sample in chip_fragment_samples) if chip_analysis_enabled else r"$.^"
    input:
        bam=lambda wildcards: chip_bam(wildcards.sample),
        spec=lambda wildcards: chip_controls_spec
    output:
        fragments=f"{chip_peak_dir}/chip_fragments/{{sample}}.bedpe"
    threads: 1
    resources:
        mem_mb=Memory_Per_Rule['10'], partition=master_config['partition'], runtime=Runtime_Per_Rule['10']
    run:
        count = bam_fragments(input.bam, output.fragments, chip_layout(wildcards.sample), chip_se_extension)
        log_it(logfile, f"MACS fragments: {wildcards.sample}, {count}, {output.fragments}", "CHIP FRAGMENTS")

def chip_region_inputs(_wildcards):
    inputs = [chip_regions_spec]
    if chip_custom_bed != "NA":
        return [*inputs, chip_custom_bed]
    if str(config.get("BROAD_MODE", "off")) in {"genebody", "diffuse"}:
        return [*inputs, f"{chip_peak_dir}/all_groups.merged_peaks.bed"]
    return [*inputs, *[chip_fragment_path(sample) for sample in chip_fragment_samples]]

rule chip_testing_regions:
    input: chip_region_inputs
    output:
        bed=chip_region_bed,
        manifest=chip_region_manifest,
        provenance=chip_region_manifest + ".source.json"
    threads: Threads_Per_Rule['10']
    resources:
        mem_mb=Memory_Per_Rule['10'], partition=master_config['partition'], runtime=Runtime_Per_Rule['10']
    run:
        reference = reference_lengths(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"])
        mode = str(config.get("BROAD_MODE", "off"))
        if chip_custom_bed != "NA":
            source = chip_custom_bed
        elif mode in {"genebody", "diffuse"}:
            source = f"{chip_peak_dir}/all_groups.merged_peaks.bed"
        else:
            command = macs_command(
                [chip_fragment_path(sample) for sample in chip_analytical_samples],
                [chip_fragment_path(sample) for sample in chip_control_sets.ids(chip_analytical_samples)],
                os.path.dirname(output.bed), "joint", q=float(config.get("CHIP_JOINT_PEAK_Q", 0.1)),
                broad=mode == "domain", genome_size=chip_genome_size, format_args=["-f", "BEDPE"],
            )
            log_it(logfile, " ".join(quote(str(part)) for part in command), "JOINT MACS COMMAND")
            log_it(logfile, subprocess.check_output(["macs3", "--version"], text=True).strip(), "MACS VERSION")
            subprocess.run(command, check=True)
            source = os.path.join(os.path.dirname(output.bed), "joint_peaks.broadPeak" if mode == "domain" else "joint_peaks.narrowPeak")
        regions = read_regions(source, reference, allow_empty=chip_custom_bed == "NA", allow_overlaps=mode == "genebody")
        write_regions(regions, output.bed, output.manifest, source, config["THEGENOME"])
        for path in output:
            Path(path).touch()
        log_it(logfile, f"Testing regions: {len(regions)}; source={source}; discovery q={config.get('CHIP_JOINT_PEAK_Q', 0.1)}. Peak overlap is annotation only.", "CHIP TESTING REGIONS")

def chip_count_outputs():
    return [f"{chip_count_dir}/{name}" for name in ("chip_input_counts.tsv", "chip_regional_enrichment.tsv", "chip_library_depths.tsv", "chip_count_provenance.json", "chip_testing_regions.bed", "chip_testing_regions.tsv", "chip_testing_regions.source.json")]

def chip_count_selection():
    selected = list(chip_analytical_samples)
    dropped = set()
    drop_file = f"{chip_filtered_dir}/peak_qc/spp_qc/dropped_samples.tsv"
    if str(config.get("SPP_GATE", "warn")).strip().lower() == "drop" and os.path.isfile(drop_file):
        with open(drop_file) as handle:
            dropped = {row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")}
        selected = [sample for sample in selected if sample_id_for_sample(sample) not in dropped]
    return selected, sorted(dropped)

def chip_count_selection_inputs(_wildcards):
    inputs = [chip_controls_spec, chip_count_settings_spec]
    if master_config["peakqc_rule_num"] in themode:
        inputs.append(f"{chip_filtered_dir}/extra_{master_config['peakqc_rule_num']}.tmp")
    drop_file = f"{chip_filtered_dir}/peak_qc/spp_qc/dropped_samples.tsv"
    if str(config.get("SPP_GATE", "warn")).strip().lower() == "drop" and os.path.isfile(drop_file):
        inputs.append(drop_file)
    return inputs

rule chip_count_selection_manifest:
    input: chip_count_selection_inputs
    output: chip_count_selection_spec
    run:
        selected, dropped = chip_count_selection()
        write_stable_text(
            output[0],
            json.dumps({"samples": selected, "dropped": dropped}, sort_keys=True) + "\n",
        )

def chip_count_run(table, threads):
    reference = reference_lengths(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"])
    regions = read_regions(chip_region_bed, reference, allow_overlaps=str(config.get("BROAD_MODE", "off")) == "genebody")
    selected, dropped = chip_count_selection()
    if dropped:
        log_it(logfile, "SPP exclusions from ChIP counts: " + ", ".join(dropped))
    if not selected:
        raise ValueError("No ChIP samples remain for differential counting.")
    counted = sorted(set(selected) | set(chip_control_sets.ids(selected)))
    libraries = {sample: {"bam": chip_bam(sample), "paired": chip_layout(sample)} for sample in counted}
    def command_log(command):
        log_it(logfile, " ".join(quote(str(part)) for part in command), "CHIP COUNT COMMAND")
    version = subprocess.check_output(["featureCounts", "-v"], stderr=subprocess.STDOUT, text=True).strip()
    log_it(logfile, version, "FEATURECOUNTS VERSION")
    log_it(logfile, subprocess.check_output(["samtools", "--version"], text=True).splitlines()[0], "SAMTOOLS VERSION")
    cache = None if config.get("RERUN_SELECTED_STEPS", False) else f"{experiment_dir}/run_logs/count_cache"
    counts, depths = count_libraries(regions, libraries, chip_count_dir, threads, command_log, cache, version)
    mode = str(config.get("BROAD_MODE", "off"))
    feature_label = "Feature" if mode == "genebody" else ("Bin" if mode == "diffuse" else "Peak")
    write_tsv(table, [feature_label, *selected], [[region_id(row), *[counts[sample][position] for sample in selected]] for position, row in enumerate(regions)])
    inputs = list(chip_control_sets.ids(selected))
    write_tsv(chip_count_outputs()[0], [feature_label, *inputs], [[region_id(row), *[counts[sample][position] for sample in inputs]] for position, row in enumerate(regions)])
    write_tsv(chip_count_outputs()[1], ["underscore", "sample_id", "input_ids", "chip_count", "input_count", "chip_depth", "input_depth", "chip_cpm", "input_cpm", "fold_enrichment", "log2_fold_enrichment", "reliability"], enrichment_rows(regions, selected, chip_control_sets, counts, depths, int(config.get("CHIP_ENRICHMENT_MIN_INPUT", 5))))
    write_tsv(chip_count_outputs()[2], ["sample_id", "role", "mapped_fragments", "read_layout"], [[sample, "chip" if sample in selected else "input", depths[sample], "PE" if libraries[sample]["paired"] else "SE"] for sample in sorted(libraries)])
    provenance = {"regions_sha256": sha256_file(chip_region_bed), "region_spec": region_spec, "featurecounts_version": version, "libraries": libraries, "controls": {sample: list(chip_control_sets.ids([sample])) for sample in selected}, "legacy_input": chip_control_sets.legacy, "enrichment_min_input_count": int(config.get("CHIP_ENRICHMENT_MIN_INPUT", 5)), "enrichment_scale": "mapped-fragment CPM; depth-weighted pooled input; no pseudocount", "counting": "unique feature assignment; SE aligned reads; PE proper same-contig fragments", "de_input": "raw ChIP counts; no input subtraction or division", "eligibility_filter": "MACS joint discovery or supplied coordinates; no overlap eligibility gate"}
    provenance.update({"project_root": experiment_dir, "spp_gate": config.get("SPP_GATE", "warn"), "spp_excluded_samples": dropped})
    write_stable_text(chip_count_outputs()[3], json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    for source, destination in zip((chip_region_bed, chip_region_manifest, chip_region_manifest + ".source.json"), chip_count_outputs()[4:]):
        shutil.copy2(source, destination)

def chip_metadata_inputs(_wildcards):
    inputs = chip_count_outputs()[4:]
    gtf = os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf")
    inputs.append(gtf)
    if master_config["callpeaks_rule_num"] in themode:
        inputs.append(f"{chip_peak_dir}/extra_10.tmp")
    else:
        inputs.extend(path for path in glob.glob(f"{chip_peak_dir}/*.MACS3.optimized.bed") if os.path.isfile(path))
    inputs.extend(glob.glob(f"{chip_peak_dir}/chip_idr_peak_calling/*idr_selected_peaks.tsv"))
    return inputs

rule chip_testing_metadata:
    input: chip_metadata_inputs
    output:
        metadata=chip_testing_metadata,
        provenance=chip_testing_metadata + ".source.json"
    threads: 1
    resources:
        mem_mb=Memory_Per_Rule['15'], partition=master_config['partition'], runtime=Runtime_Per_Rule['15']
    run:
        regions = read_regions(chip_count_outputs()[4], reference_lengths(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"]), allow_overlaps=str(config.get("BROAD_MODE", "off")) == "genebody")
        catalogue = {sample_type_for_sample(sample): f"{chip_peak_dir}/{sample_type_for_sample(sample)}.MACS3.optimized.bed" for sample in chip_analytical_samples}
        support_header, support_rows = support_annotations(regions, catalogue)
        methods = {}
        summary = f"{chip_peak_dir}/chip_idr_peak_calling/idr_selected_peaks.tsv"
        if os.path.isfile(summary):
            with open(summary) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row.get("replicate_1") == "all" and row.get("replicate_2") == "all":
                        methods[row["group"]] = {"method": row["strategy"], "note": row["note"]}
        for group, path in sorted(catalogue.items()):
            assessed = os.path.isfile(path)
            method = methods.get(group, {"method": "group_peak_catalogue", "note": "IDR selection provenance unavailable"}) if assessed else {"method": "unassessed", "note": "catalogue unavailable"}
            method["catalogue_sha256"] = sha256_file(path) if assessed else None
            methods[group] = method
            support_header.extend([f"{group}.support_method", f"{group}.support_note"])
            for row in support_rows:
                row.extend([method["method"], method["note"]])
        with open(chip_count_outputs()[5]) as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        with tempfile.TemporaryDirectory(prefix="omnom_chip_annotation.") as temporary:
            annotation_header, annotation_rows = annotate_regions(regions, os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf"), temporary)
        header = list(records[0]) + annotation_header[1:] + support_header[1:]
        if str(config.get("BROAD_MODE", "off")) == "genebody":
            for annotation in annotation_rows:
                annotation[1] = "gene_body"
        rows = [[*record.values(), *annotation[1:], *support[1:]] for record, annotation, support in zip(records, annotation_rows, support_rows)]
        write_tsv(output.metadata, header, rows)
        write_stable_text(output.provenance, json.dumps({"region_sha256": sha256_file(chip_count_outputs()[4]), "peak_catalogues": methods, "overlap_rule": "any positive overlap; annotation only", "multiple_testing": "full tested family; support subsets do not refit or adjust p-values"}, sort_keys=True, indent=2) + "\n")
        for path in output:
            Path(path).touch()

chip_track_mode = str(config.get("CHIP_INPUT_TRACKS", "fold_enrichment"))

def chip_input_track_path(sample):
    return f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num'] - 1]}/input_relative/{sample}.{chip_track_mode}.bw"

rule chip_input_track:
    wildcard_constraints:
        sample="|".join(re.escape(sample) for sample in chip_analytical_samples if chip_control_sets.ids([sample])) if chip_analysis_enabled else r"$.^"
    input:
        fragments=lambda wildcards: chip_fragment_path(wildcards.sample),
        controls=lambda wildcards: chip_control_dependencies([wildcards.sample]),
        bam=lambda wildcards: chip_bam(wildcards.sample),
        control_bams=lambda wildcards: chip_control_sets.paths([wildcards.sample])
    output:
        bw=f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num'] - 1]}/input_relative/{{sample}}.{chip_track_mode}.bw",
        provenance=f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num'] - 1]}/input_relative/{{sample}}.{chip_track_mode}.json"
    threads: Threads_Per_Rule['8']
    resources:
        mem_mb=Memory_Per_Rule['8'], partition=master_config['partition'], runtime=Runtime_Per_Rule['8']
    run:
        if chip_track_mode == "none":
            raise ValueError("No input-relative track requested.")
        os.makedirs(os.path.dirname(output.bw), exist_ok=True)
        reference = reference_lengths(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"])
        controls = [chip_fragment_path(sample) for sample in chip_control_sets.ids([wildcards.sample])]
        commands = []
        def execute(command):
            commands.append(command)
            log_it(logfile, " ".join(quote(str(part)) for part in command), "CHIP INPUT TRACK COMMAND")
            subprocess.run(command, check=True)
        with tempfile.TemporaryDirectory(prefix="omnom_chip_track.") as temporary:
            if chip_track_mode == "log2_ratio":
                pooled = os.path.join(temporary, "input.bam")
                execute(["samtools", "merge", "-f", pooled, *map(str, input.control_bams)])
                # bamCompare provides a direct read-coverage ratio, distinct from MACS background FE.
                chip_copy = os.path.join(temporary, "chip.bam")
                shutil.copy2(input.bam, chip_copy)
                execute(["samtools", "index", chip_copy])
                execute(["samtools", "index", pooled])
                execute(["bamCompare", "-b1", chip_copy, "-b2", pooled, "-o", str(output.bw), "--operation", "log2", "--scaleFactorsMethod", "None", "--normalizeUsing", "CPM", "--pseudocount", "1", "--binSize", "10", "--numberOfProcessors", str(threads)])
                scale = "direct CPM read-coverage log2 ratio; pseudocount 1 CPM; 10bp bins"
            else:
                execute(macs_command([str(input.fragments)], controls, temporary, "track", q=0.1, genome_size=chip_genome_size, format_args=["-f", "BEDPE"], tracks="spmr" if chip_track_mode == "fold_enrichment" else True))
                mode = "FE" if chip_track_mode == "fold_enrichment" else "qpois"
                graph = os.path.join(temporary, "signal.bdg")
                execute(["macs3", "bdgcmp", "-t", os.path.join(temporary, "track_treat_pileup.bdg"), "-c", os.path.join(temporary, "track_control_lambda.bdg"), "-m", mode, "-o", graph])
                sorted_graph = os.path.join(temporary, "signal.sorted.bdg")
                with open(sorted_graph, "w") as handle:
                    subprocess.run(["sort", "-k1,1", "-k2,2n", graph], stdout=handle, check=True)
                sizes = os.path.join(temporary, "chrom.sizes")
                Path(sizes).write_text("".join(f"{chrom}\t{length}\n" for chrom, length in sorted(reference.items())))
                bedgraph_bigwig(sorted_graph, output.bw, reference)
                scale = "MACS depth-scaled modelled local background fold enrichment" if mode == "FE" else "MACS background -log10(q); unscaled fragment pileups"
            provenance = {"sample": wildcards.sample, "inputs": list(chip_control_sets.ids([wildcards.sample])), "mode": chip_track_mode, "scale": scale, "commands": commands, "macs_version": subprocess.check_output(["macs3", "--version"], text=True).strip(), "de_use": "visualization only; not DE counts"}
            write_stable_text(output.provenance, json.dumps(provenance, sort_keys=True, indent=2) + "\n")
