import argparse
import gzip
import shutil
import sys
from pathlib import Path

import yaml


def check_config_file_header(config_file_path, expected_header):
    with open(config_file_path, "r") as config_file:
        lines = config_file.readlines()
        if len(lines) < 3 or expected_header not in lines[2]:
            print("Config file does not contain right header. Aborting...", file=sys.stderr)
            sys.exit(1)


def load_and_validate_yaml(config_file_path, expected_header):
    check_config_file_header(config_file_path, expected_header)
    with open(config_file_path, "r") as file:
        try:
            return yaml.safe_load(file)
        except yaml.YAMLError as exc:
            print(f"Error parsing YAML file: {exc}", file=sys.stderr)
            sys.exit(1)


def merge_configs(base_config, override_config):
    merged_config = dict(base_config)
    merged_config.update(override_config)
    return merged_config


def resolve_config_path(path_value, workflow_root):
    resolved = str(path_value)
    resolved = resolved.replace("{WORKFLOW_ROOT}", str(workflow_root))
    resolved = resolved.replace("{HOME}", str(Path.home()))
    return Path(resolved).expanduser().resolve()


def load_site_settings(workflow_root, workflow_config_file, site_config_file):
    workflow_config = load_and_validate_yaml(workflow_config_file, "## Omnomnomics pipeline config ##")
    site_config = load_and_validate_yaml(site_config_file, "## Omnomnomics pipeline config ##")
    config = merge_configs(workflow_config, site_config)
    config["genome_assembly_dir"] = resolve_config_path(config["genome_assembly_dir"], workflow_root)
    config["cellranger_reference_dir"] = resolve_config_path(config["cellranger_reference_dir"], workflow_root)
    return config


def parse_genomes_arguments(argv):
    parser = argparse.ArgumentParser(
        prog="omnomnomics genomes",
        description="Manage reference genomes for omnomnomics.",
        allow_abbrev=False,
    )
    parser.add_argument(
        "--site-config",
        help="Optional path to a site-specific config YAML. Default: packaged site config",
    )

    subparsers = parser.add_subparsers(dest="genomes_command", required=True)

    list_parser = subparsers.add_parser("list", help="List available assemblies for a species")
    list_parser.add_argument("--species", default="homo sapiens", help="Species search term")
    list_parser.add_argument("--provider", default="GENCODE", help="Provider name for genomepy search")
    list_parser.add_argument("--limit", type=int, default=25, help="Maximum number of rows to show")

    install_parser = subparsers.add_parser("install", help="Download and normalize one or more assemblies")
    install_parser.add_argument("--species", default="homo sapiens", help="Species search term if no assemblies are given")
    install_parser.add_argument(
        "--assembly",
        action="append",
        default=[],
        help="Assembly name to install. Repeatable. If omitted, resolves the latest matching assembly for the species",
    )
    install_parser.add_argument("--provider", default="GENCODE", help="Provider name for genomepy install/search")
    install_parser.add_argument("--threads", type=int, default=8, help="Threads for index creation")
    install_parser.add_argument("--force", action="store_true", help="Overwrite an existing normalized assembly")
    install_parser.add_argument("--keep-alt", action="store_true", help="Keep alternative contigs if supported by provider")
    install_parser.add_argument(
        "--indexers",
        nargs="+",
        default=["hisat2", "star"],
        choices=["hisat2", "star"],
        help="Which aligner indexes to require and normalize",
    )
    install_parser.add_argument(
        "--mask",
        default="soft",
        choices=["soft", "hard", "none"],
        help="Masking mode passed to genomepy",
    )
    install_parser.add_argument(
        "--ucsc-annotation",
        help="UCSC annotation style passed through to genomepy, e.g. refGene or ensGene",
    )
    install_parser.add_argument("--dry-run", action="store_true", help="Resolve targets and print actions without downloading")

    return parser.parse_args(argv)


def import_genomepy():
    try:
        import genomepy
    except ImportError:
        print(
            "The genome helper requires the 'genomepy' package in the active environment. Aborting...",
            file=sys.stderr,
        )
        sys.exit(1)
    return genomepy


def search_rows(genomepy, species, provider):
    rows = list(genomepy.search(species, provider=provider))
    if not rows:
        print(f"No assemblies found for species query '{species}' on provider '{provider}'.", file=sys.stderr)
        sys.exit(1)
    return rows


def provider_rows(genomepy, provider):
    rows = list(genomepy.search("", provider=provider))
    if not rows:
        print(f"No assemblies could be retrieved from provider '{provider}'.", file=sys.stderr)
        sys.exit(1)
    return rows


def resolve_default_assembly(genomepy, species, provider):
    rows = search_rows(genomepy, species, provider)
    species_norm = species.strip().lower()

    def score(row):
        row_species = str(row[5]).strip().lower() if len(row) > 5 else ""
        name = str(row[0]).strip().lower()
        exact_species = row_species == species_norm
        starts_with_species = row_species.startswith(species_norm)
        name_matches_species = name == species_norm.replace(" ", "_")
        return (exact_species, starts_with_species, name_matches_species)

    rows = sorted(rows, key=score, reverse=True)
    return rows[0][0], rows


def print_rows(rows, limit):
    header = ["name", "provider", "accession", "tax_id", "annotation", "species", "other_info"]
    print("\t".join(header))
    for row in rows[:limit]:
        print("\t".join(str(field) for field in row))


def validate_requested_assemblies(genomepy, assembly_names, provider):
    rows = provider_rows(genomepy, provider)
    available = {str(row[0]): row for row in rows}

    for assembly_name in assembly_names:
        if assembly_name in available:
            continue

        assembly_name_lower = assembly_name.lower()
        suggestions = [
            row[0]
            for row in rows
            if assembly_name_lower in str(row[0]).lower() or str(row[0]).lower() in assembly_name_lower
        ]
        suggestion_text = f" Did you mean: {', '.join(suggestions[:5])}?" if suggestions else ""
        print(
            f"Assembly '{assembly_name}' is not available from provider '{provider}'.{suggestion_text}",
            file=sys.stderr,
        )
        sys.exit(1)


def ensure_plugins(genomepy):
    genomepy.manage_plugins("enable", ["star", "hisat2"])


def copy_text_file(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_annotation(annotation_file, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if str(annotation_file).endswith(".gz"):
        with gzip.open(annotation_file, "rt") as src_handle, open(dst, "w") as dst_handle:
            shutil.copyfileobj(src_handle, dst_handle)
    else:
        shutil.copy2(annotation_file, dst)


def parse_gtf_attributes(attribute_text):
    attributes = {}
    for field in attribute_text.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, value = field.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def create_gene_bed_files(gtf_file, assembly_name, aux_dir):
    detailed_rows = []
    legacy_rows = []

    with open(gtf_file, "r") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id") or attributes.get("ID")
            gene_symbol = (
                attributes.get("gene_name")
                or attributes.get("gene_symbol")
                or attributes.get("Name")
                or gene_id
            )

            if not gene_id:
                continue

            detailed_rows.append((chrom, start, end, gene_id, "0", strand, gene_symbol))
            legacy_rows.append((chrom, start, end, gene_id, "0", strand))

    detailed_rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    legacy_rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))

    detailed_path = aux_dir / f"{assembly_name}_genes.sorted.bed"
    legacy_path = aux_dir / f"{assembly_name}_refseq_genes.sorted.bed"

    with open(detailed_path, "w") as detailed_handle:
        for row in detailed_rows:
            detailed_handle.write("\t".join(map(str, row)) + "\n")

    with open(legacy_path, "w") as legacy_handle:
        for row in legacy_rows:
            legacy_handle.write("\t".join(map(str, row)) + "\n")


def genome_install_root(genome):
    if getattr(genome, "genome_dir", None):
        return Path(genome.genome_dir)
    return Path(genome.filename).resolve().parent


def has_complete_star_index(index_dir):
    required = ["Genome", "SA", "SAindex", "genomeParameters.txt", "chrLength.txt", "chrName.txt"]
    index_path = Path(index_dir)
    return index_path.is_dir() and all((index_path / name).exists() for name in required)


def has_complete_hisat2_index(index_dir, assembly_name):
    index_path = Path(index_dir)
    required = [index_path / f"{assembly_name}.{i}.ht2" for i in range(1, 9)]
    no_intermediates = not any(index_path.glob("*.rf"))
    return index_path.is_dir() and all(path.exists() for path in required) and no_intermediates


def find_star_dir(genome_dir):
    markers = {"Genome", "SA", "SAindex", "genomeParameters.txt"}
    candidates = []
    for path in Path(genome_dir).rglob("*"):
        if path.is_dir():
            child_names = {child.name for child in path.iterdir()}
            if markers.issubset(child_names):
                candidates.append(path)
    if not candidates:
        return None
    return sorted(candidates, key=lambda path: (len(path.parts), str(path)))[0]


def find_hisat2_dir(genome_dir):
    candidates = []
    for path in Path(genome_dir).rglob("*"):
        if path.is_dir():
            ht2_files = list(path.glob("*.ht2")) + list(path.glob("*.ht2l"))
            if ht2_files:
                candidates.append(path)
    if not candidates:
        return None
    return sorted(candidates, key=lambda path: (len(path.parts), str(path)))[0]


def normalize_genome_install(genome, assembly_name, assembly_root, force, indexers):
    final_root = Path(assembly_root) / assembly_name
    if final_root.exists():
        if not force:
            print(
                f"Assembly '{assembly_name}' already exists at {final_root}. Use --force to overwrite.",
                file=sys.stderr,
            )
            sys.exit(1)
        shutil.rmtree(final_root)

    fasta_dir = final_root / "fasta"
    annotation_dir = final_root / "annotation"
    aux_dir = final_root / "aux"
    star_dir = final_root / "star"
    hisat2_dir = final_root / "hisat2"
    notes_dir = final_root / "notes"

    fasta_dir.mkdir(parents=True, exist_ok=True)
    annotation_dir.mkdir(parents=True, exist_ok=True)
    aux_dir.mkdir(parents=True, exist_ok=True)
    notes_dir.mkdir(parents=True, exist_ok=True)

    copy_text_file(genome.genome_file, fasta_dir / "genome.fa")
    copy_text_file(genome.index_file, fasta_dir / "genome.fa.fai")
    copy_text_file(genome.sizes_file, aux_dir / f"{assembly_name}_chrom_sizes.2_column")
    copy_text_file(genome.readme_file, notes_dir / "genomepy.README.txt")

    if getattr(genome, "annotation_gtf_file", None):
        copy_annotation(genome.annotation_gtf_file, annotation_dir / "genes.gtf")
        create_gene_bed_files(annotation_dir / "genes.gtf", assembly_name, aux_dir)

    install_root = genome_install_root(genome)

    if "star" in indexers:
        star_source = find_star_dir(install_root)
        if star_source is None or not has_complete_star_index(star_source):
            print(f"Incomplete STAR index for '{assembly_name}' under {install_root}", file=sys.stderr)
            sys.exit(1)
        shutil.copytree(star_source, star_dir, dirs_exist_ok=True)

    if "hisat2" in indexers:
        hisat2_source = find_hisat2_dir(install_root)
        if hisat2_source is None or not has_complete_hisat2_index(hisat2_source, assembly_name):
            print(f"Incomplete HISAT2 index for '{assembly_name}' under {install_root}", file=sys.stderr)
            sys.exit(1)
        shutil.copytree(hisat2_source, hisat2_dir, dirs_exist_ok=True)

    return final_root


def install_assembly(genomepy, assembly_name, provider, assembly_root, threads, force, keep_alt, mask, ucsc_annotation, indexers):
    assembly_root = Path(assembly_root)
    assembly_root.mkdir(parents=True, exist_ok=True)
    tmp_root = assembly_root / ".genomepy_tmp"
    tmp_root.mkdir(parents=True, exist_ok=True)

    install_kwargs = {
        "annotation": True,
        "provider": provider,
        "genomes_dir": str(tmp_root),
        "localname": assembly_name,
        "threads": threads,
        "force": force,
        "keep_alt": keep_alt,
        "mask": mask,
    }
    if ucsc_annotation:
        install_kwargs["ucsc_annotation"] = ucsc_annotation

    ensure_plugins(genomepy)
    wanted_indexers = set(indexers)
    for plugin_name in ["hisat2", "star"]:
        if plugin_name in wanted_indexers:
            genomepy.manage_plugins("enable", [plugin_name])
        else:
            genomepy.manage_plugins("disable", [plugin_name])

    genome = genomepy.install_genome(assembly_name, **install_kwargs)
    normalized_root = normalize_genome_install(genome, assembly_name, assembly_root, force, indexers)

    tmp_assembly_root = tmp_root / assembly_name
    if tmp_assembly_root.exists():
        shutil.rmtree(tmp_assembly_root)

    return normalized_root


def genomes_main(argv, workflow_root, workflow_config_file, default_site_config):
    args = parse_genomes_arguments(argv)
    site_config_file = Path(args.site_config).expanduser().resolve() if args.site_config else default_site_config
    config = load_site_settings(workflow_root, workflow_config_file, site_config_file)
    genomepy = import_genomepy()

    if args.genomes_command == "list":
        rows = search_rows(genomepy, args.species, args.provider)
        print_rows(rows, args.limit)
        return

    targets = list(args.assembly)
    if not targets:
        resolved_assembly, rows = resolve_default_assembly(genomepy, args.species, args.provider)
        print(f"Resolved default assembly for '{args.species}' on provider '{args.provider}' to '{resolved_assembly}'.")
        print_rows(rows, 5)
        targets = [resolved_assembly]
    else:
        validate_requested_assemblies(genomepy, targets, args.provider)

    for assembly_name in targets:
        final_root = Path(config["genome_assembly_dir"]) / assembly_name
        print(f"Preparing assembly '{assembly_name}' in {final_root}")
        if args.dry_run:
            continue
        installed_path = install_assembly(
            genomepy=genomepy,
            assembly_name=assembly_name,
            provider=args.provider,
            assembly_root=config["genome_assembly_dir"],
            threads=args.threads,
            force=args.force,
            keep_alt=args.keep_alt,
            indexers=args.indexers,
            mask=args.mask,
            ucsc_annotation=args.ucsc_annotation,
        )
        print(f"Installed '{assembly_name}' to {installed_path}")
