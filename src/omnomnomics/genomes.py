import argparse
import gzip
import math
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import yaml


def check_config_file_header(config_file_path, expected_header):
    with open(config_file_path, "r") as config_file:
        header_lines = [line.strip() for line in config_file if line.strip()]
        if len(header_lines) < 1 or expected_header not in header_lines[0]:
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


SPECIES_ALIASES = {
    "human": "homo sapiens",
    "mouse": "mus musculus",
}
INSTALLED_SPECIES_HINTS = {
    "homo sapiens": ("grch", "hg", "human", "homo_sapiens", "homo sapiens"),
    "mus musculus": ("grcm", "mm", "mouse", "mus_musculus", "mus musculus"),
}
UCSC_BLACKLIST_ASSEMBLY_ALIASES = {
    "GRCh38": "hg38",
    "GRCh38.p14": "hg38",
    "GRCm38": "mm10",
    "GRCm39": "mm39",
}
DEFAULT_MEME_MOTIF_DATABASE_URL = (
    "https://jaspar.elixir.no/download/data/2024/CORE/"
    "JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
)
DEFAULT_MEME_MOTIF_DATABASE_NAME = "JASPAR2024_CORE_vertebrates_non-redundant.meme"
MEME_MOTIF_DATABASE_ENV_VARS = (
    "OMNOMNOMICS_MEME_MOTIF_DATABASE",
    "MEME_MOTIF_DATABASE",
    "JASPAR_MOTIF_DATABASE",
)


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

    list_parser = subparsers.add_parser("list", aliases=["ls"], help="List available assemblies for a species")
    list_parser.add_argument("--species", default="homo sapiens", help="Species search term")
    list_parser.add_argument("--provider", default="GENCODE", help="Provider name for genomepy search")
    list_parser.add_argument("--limit", type=int, default=25, help="Maximum number of rows to show")

    installed_parser = subparsers.add_parser("installed", aliases=["local"], help="List locally installed assemblies")
    installed_parser.add_argument(
        "--species",
        help="Optional species hint to filter installed assemblies by assembly name, e.g. human or mouse",
    )
    installed_parser.add_argument("--limit", type=int, default=100, help="Maximum number of rows to show")

    blacklist_parser = subparsers.add_parser(
        "blacklist",
        help="Download and cache an ENCODE blacklist BED for an installed assembly",
    )
    blacklist_parser.add_argument("--assembly", required=True, help="Installed omnomnomics assembly name")
    blacklist_parser.add_argument("--provider", default="UCSC", help="Provider name for blacklist lookup")
    blacklist_parser.add_argument("--threads", type=int, default=1, help="Threads passed to genomepy")
    blacklist_parser.add_argument("--force", action="store_true", help="Refresh an existing cached blacklist BED")

    motifs_parser = subparsers.add_parser(
        "motifs",
        aliases=["motif-db"],
        help="Download or show the cached MEME motif database used by post-DE peak analysis",
    )
    motifs_parser.add_argument("--force", action="store_true", help="Refresh the cached motif database")
    motifs_parser.add_argument("--url", default=DEFAULT_MEME_MOTIF_DATABASE_URL, help="MEME motif database URL")
    motifs_parser.add_argument("--name", default=DEFAULT_MEME_MOTIF_DATABASE_NAME, help="Cached database filename")
    motifs_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the cache path without downloading",
    )

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
        "--skip-blacklist",
        action="store_true",
        help="Do not cache a genomepy-provided blacklist BED during assembly install",
    )
    install_parser.add_argument(
        "--skip-motif-db",
        action="store_true",
        help="Do not cache the default MEME motif database during assembly install",
    )
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


def normalize_species_name(species):
    return SPECIES_ALIASES.get(species.strip().lower(), species)


def import_genomepy():
    try:
        return load_genomepy()
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)


def load_genomepy():
    try:
        import genomepy
    except ImportError as exc:
        raise RuntimeError(
            "The genome helper requires the 'genomepy' package in the active environment. Aborting...",
        ) from exc
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


def filter_installed_rows_by_species(rows, species):
    if not species:
        return rows

    species_norm = normalize_species_name(species).strip().lower()
    hint_tokens = INSTALLED_SPECIES_HINTS.get(species_norm)
    if not hint_tokens:
        compact_species = species_norm.replace(" ", "_")
        hint_tokens = (species_norm, compact_species)

    filtered_rows = [
        row
        for row in rows
        if any(token in str(row[0]).lower() for token in hint_tokens)
    ]
    return filtered_rows


def list_installed_rows(assembly_root, species=None):
    assembly_root = Path(assembly_root)
    if not assembly_root.is_dir():
        print(f"Genome assembly root '{assembly_root}' does not exist.", file=sys.stderr)
        sys.exit(1)

    rows = []
    for path in sorted(assembly_root.iterdir()):
        if not path.is_dir() or path.name.startswith("."):
            continue
        fasta_ok = (path / "fasta" / "genome.fa").is_file()
        gtf_ok = (path / "annotation" / "genes.gtf").is_file()
        star_ok = has_complete_star_index(path / "star")
        hisat2_ok = has_complete_hisat2_index(path / "hisat2", path.name)
        rows.append(
            (
                path.name,
                "yes" if fasta_ok else "no",
                "yes" if gtf_ok else "no",
                "yes" if star_ok else "no",
                "yes" if hisat2_ok else "no",
                str(path),
            )
        )

    return filter_installed_rows_by_species(rows, species)


def print_installed_rows(rows, limit):
    header = ["name", "fasta", "gtf", "star", "hisat2", "path"]
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


def blacklist_aux_path(assembly_root, assembly_name):
    return Path(assembly_root) / assembly_name / "aux" / f"{assembly_name}.blacklist.bed"


def find_cached_blacklist_bed(assembly_root, assembly_name):
    assembly_dir = Path(assembly_root) / assembly_name
    aux_dir = assembly_dir / "aux"
    candidates = [
        blacklist_aux_path(assembly_root, assembly_name),
        aux_dir / f"{assembly_name}.blacklist.bed.gz",
    ]
    candidates.extend(find_blacklist_beds(aux_dir))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def find_blacklist_beds(root):
    return sorted(Path(root).rglob("*blacklist*.bed*"))


def copy_blacklist_bed(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    src = Path(src)
    dst = Path(dst)
    if src.resolve() == dst.resolve():
        return dst
    if src.name.endswith(".gz"):
        with gzip.open(src, "rt") as src_handle, open(dst, "w") as dst_handle:
            shutil.copyfileobj(src_handle, dst_handle)
    else:
        shutil.copy2(src, dst)
    return dst


def meme_motif_database_dir(assembly_root):
    return Path(assembly_root) / "motif_databases"


def meme_motif_database_path(assembly_root, name=DEFAULT_MEME_MOTIF_DATABASE_NAME):
    return meme_motif_database_dir(assembly_root) / name


def looks_like_meme_database(path):
    try:
        with open(path, "r", errors="replace") as handle:
            head = "".join(handle.readline() for _ in range(50))
    except OSError:
        return False
    return "MEME version" in head and "ALPHABET" in head


def find_cached_meme_motif_database(assembly_root, name=DEFAULT_MEME_MOTIF_DATABASE_NAME):
    primary = meme_motif_database_path(assembly_root, name)
    candidates = [primary]
    motif_dir = meme_motif_database_dir(assembly_root)
    if motif_dir.is_dir():
        candidates.extend(sorted(motif_dir.glob("*.meme")))
        candidates.extend(sorted(motif_dir.glob("*_meme.txt")))
        candidates.extend(sorted(motif_dir.glob("*.txt")))

    seen = set()
    for candidate in candidates:
        candidate = Path(candidate)
        if candidate in seen:
            continue
        seen.add(candidate)
        if candidate.is_file() and candidate.stat().st_size > 0 and looks_like_meme_database(candidate):
            return candidate
    return None


def find_environment_meme_motif_database():
    for env_var in MEME_MOTIF_DATABASE_ENV_VARS:
        env_value = os.environ.get(env_var, "").strip()
        if not env_value:
            continue
        env_path = Path(env_value).expanduser()
        if env_path.is_file() and looks_like_meme_database(env_path):
            return env_path

    search_roots = []
    try:
        import glob

        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            search_roots.extend(
                [
                    Path(conda_prefix) / "share" / "meme" / "motif_databases",
                    Path(conda_prefix) / "share" / "meme" / "db" / "motif_databases",
                    Path(conda_prefix) / "share" / "meme" / "db",
                    Path(conda_prefix) / "share" / "meme",
                ]
            )
            search_roots.extend(
                Path(path)
                for path in glob.glob(str(Path(conda_prefix) / "share" / "meme-*" / "db" / "motif_databases"))
            )
            search_roots.extend(
                Path(path)
                for path in glob.glob(str(Path(conda_prefix) / "share" / "meme-*" / "motif_databases"))
            )
    except Exception:
        pass
    search_roots.extend([Path("/usr/local/share"), Path("/usr/share")])

    patterns = [
        "JASPAR*CORE*vertebrates*non-redundant*.meme",
        "JASPAR*CORE*vertebrates*.meme",
        "HOCOMOCO*.meme",
        "*.meme",
    ]
    for root in search_roots:
        if not root.is_dir():
            continue
        for pattern in patterns:
            for candidate in sorted(root.rglob(pattern)):
                candidate_parts = set(candidate.parts)
                if "doc" in candidate_parts or "examples" in candidate_parts:
                    continue
                if candidate.is_file() and candidate.stat().st_size > 0 and looks_like_meme_database(candidate):
                    return candidate
    return None


def copy_meme_motif_database(src, dst):
    src = Path(src).expanduser()
    dst = Path(dst)
    if not src.is_file() or not looks_like_meme_database(src):
        raise FileNotFoundError(f"MEME motif database '{src}' was not found or is not a valid MEME-format file.")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    return dst


def download_meme_motif_database(url, dst, force=False):
    dst = Path(dst)
    if dst.is_file() and not force and looks_like_meme_database(dst):
        return dst

    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = dst.with_name(f"{dst.name}.{os.getpid()}.tmp")
    try:
        with urllib.request.urlopen(url, timeout=120) as response, open(tmp_path, "wb") as handle:
            shutil.copyfileobj(response, handle)
        if not looks_like_meme_database(tmp_path):
            raise ValueError(f"Downloaded motif database from '{url}' is not a valid MEME-format file.")
        tmp_path.replace(dst)
    except (OSError, urllib.error.URLError, ValueError):
        if tmp_path.exists():
            tmp_path.unlink()
        raise
    return dst


def resolve_meme_motif_database(
    assembly_root,
    requested="auto",
    force=False,
    download=True,
    url=DEFAULT_MEME_MOTIF_DATABASE_URL,
    name=DEFAULT_MEME_MOTIF_DATABASE_NAME,
):
    requested = str(requested or "auto").strip()
    if requested and requested.lower() != "auto":
        requested_path = Path(requested).expanduser()
        if requested_path.is_file() and looks_like_meme_database(requested_path):
            return requested_path
        return None

    cached = find_cached_meme_motif_database(assembly_root, name=name)
    if cached and not force:
        return cached

    cache_path = meme_motif_database_path(assembly_root, name=name)
    env_db = find_environment_meme_motif_database()
    if env_db and not force:
        return copy_meme_motif_database(env_db, cache_path)

    if download:
        return download_meme_motif_database(url, cache_path, force=force)

    return None


def cache_blacklist_from_install(install_root, assembly_root, assembly_name):
    candidates = find_blacklist_beds(install_root)
    if not candidates:
        return None
    return copy_blacklist_bed(candidates[0], blacklist_aux_path(assembly_root, assembly_name))


def resolve_blacklist_bed(assembly_name, assembly_root, provider="UCSC", threads=1, force=False):
    assembly_root = Path(assembly_root)
    assembly_dir = assembly_root / assembly_name
    if not assembly_dir.is_dir():
        raise FileNotFoundError(f"Assembly '{assembly_name}' is not installed under {assembly_root}.")

    cached = find_cached_blacklist_bed(assembly_root, assembly_name)
    if cached and not force:
        return cached

    genomepy = load_genomepy()
    provider_assembly = UCSC_BLACKLIST_ASSEMBLY_ALIASES.get(assembly_name, assembly_name)
    tmp_root = assembly_root / ".genomepy_blacklist_tmp"
    local_name = f"{assembly_name}.omnomnomics.blacklist"
    tmp_genome_dir = tmp_root / local_name
    tmp_candidates = find_blacklist_beds(tmp_genome_dir)
    if tmp_candidates and not force:
        cached_path = copy_blacklist_bed(tmp_candidates[0], blacklist_aux_path(assembly_root, assembly_name))
        shutil.rmtree(tmp_genome_dir)
        return cached_path

    tmp_root.mkdir(parents=True, exist_ok=True)
    genomepy.manage_plugins("disable", ["hisat2", "star"])
    genomepy.manage_plugins("enable", ["blacklist"])
    genomepy.install_genome(
        provider_assembly,
        provider=provider,
        genomes_dir=str(tmp_root),
        localname=local_name,
        annotation=False,
        threads=threads,
        force=force,
    )

    candidates = sorted(tmp_genome_dir.rglob("*blacklist*.bed*"))
    if not candidates:
        raise FileNotFoundError(
            f"No blacklist BED found after genomepy install for '{assembly_name}' in {tmp_genome_dir}."
        )

    cached_path = copy_blacklist_bed(candidates[0], blacklist_aux_path(assembly_root, assembly_name))
    if tmp_genome_dir.exists():
        shutil.rmtree(tmp_genome_dir)
    return cached_path


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
    required_ht2 = [index_path / f"{assembly_name}.{i}.ht2" for i in range(1, 9)]
    required_ht2l = [index_path / f"{assembly_name}.{i}.ht2l" for i in range(1, 9)]
    no_intermediates = not any(index_path.glob("*.rf"))
    has_small_index = all(path.exists() for path in required_ht2)
    has_large_index = all(path.exists() for path in required_ht2l)
    return index_path.is_dir() and (has_small_index or has_large_index) and no_intermediates


def require_executable(name):
    executable = shutil.which(name)
    if not executable:
        print(
            f"Requested indexer '{name}' is not available in the active environment.",
            file=sys.stderr,
        )
        sys.exit(1)
    return executable


def fasta_total_bases(fai_path):
    total = 0
    with open(fai_path, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                total += int(fields[1])
    return total


def star_sa_index_nbases(fai_path):
    total_bases = fasta_total_bases(fai_path)
    if total_bases <= 0:
        return 1
    return max(1, min(14, int(math.log2(total_bases) / 2 - 1)))


def run_checked_command(command, label):
    print(f"Running {label}: {' '.join(map(str, command))}", file=sys.stderr)
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"{label} failed with exit code {exc.returncode}", file=sys.stderr)
        sys.exit(exc.returncode)


def build_star_index(star_dir, fasta_path, gtf_path, fai_path, threads):
    require_executable("STAR")
    star_dir = Path(star_dir)
    star_dir.mkdir(parents=True, exist_ok=True)
    command = [
        "STAR",
        "--runThreadN",
        str(threads),
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        str(star_dir),
        "--genomeFastaFiles",
        str(fasta_path),
        "--genomeSAindexNbases",
        str(star_sa_index_nbases(fai_path)),
    ]
    if gtf_path and Path(gtf_path).exists():
        command.extend(["--sjdbGTFfile", str(gtf_path), "--sjdbOverhang", "100"])
    run_checked_command(command, "STAR genomeGenerate")


def build_hisat2_index(hisat2_dir, fasta_path, assembly_name, threads):
    require_executable("hisat2-build")
    hisat2_dir = Path(hisat2_dir)
    hisat2_dir.mkdir(parents=True, exist_ok=True)
    run_checked_command(
        [
            "hisat2-build",
            "-p",
            str(threads),
            str(fasta_path),
            str(hisat2_dir / assembly_name),
        ],
        "hisat2-build",
    )


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


def normalize_genome_install(genome, assembly_name, assembly_root, force, indexers, threads, cache_blacklist=True):
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
    if cache_blacklist:
        cache_blacklist_from_install(install_root, assembly_root, assembly_name)

    if "star" in indexers:
        star_source = find_star_dir(install_root)
        if star_source is not None and has_complete_star_index(star_source):
            shutil.copytree(star_source, star_dir, dirs_exist_ok=True)
        else:
            print(f"Building STAR index for '{assembly_name}' because genomepy did not provide one.", file=sys.stderr)
            build_star_index(
                star_dir,
                fasta_dir / "genome.fa",
                annotation_dir / "genes.gtf",
                fasta_dir / "genome.fa.fai",
                threads,
            )
        if not has_complete_star_index(star_dir):
            print(f"Incomplete STAR index for '{assembly_name}' under {install_root}", file=sys.stderr)
            sys.exit(1)

    if "hisat2" in indexers:
        hisat2_source = find_hisat2_dir(install_root)
        if hisat2_source is not None and has_complete_hisat2_index(hisat2_source, assembly_name):
            shutil.copytree(hisat2_source, hisat2_dir, dirs_exist_ok=True)
        else:
            print(f"Building HISAT2 index for '{assembly_name}' because genomepy did not provide one.", file=sys.stderr)
            build_hisat2_index(hisat2_dir, fasta_dir / "genome.fa", assembly_name, threads)
        if not has_complete_hisat2_index(hisat2_dir, assembly_name):
            print(f"Incomplete HISAT2 index for '{assembly_name}' under {install_root}", file=sys.stderr)
            sys.exit(1)

    return final_root


def install_assembly(
    genomepy,
    assembly_name,
    provider,
    assembly_root,
    threads,
    force,
    keep_alt,
    mask,
    ucsc_annotation,
    indexers,
    cache_blacklist=True,
):
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

    wanted_indexers = set(indexers)
    enabled_plugins = [plugin_name for plugin_name in ["hisat2", "star"] if plugin_name in wanted_indexers]
    if cache_blacklist:
        enabled_plugins.append("blacklist")
    disabled_plugins = [plugin_name for plugin_name in ["hisat2", "star"] if plugin_name not in wanted_indexers]

    if disabled_plugins:
        genomepy.manage_plugins("disable", disabled_plugins)
    if enabled_plugins:
        genomepy.manage_plugins("enable", enabled_plugins)

    genome = genomepy.install_genome(assembly_name, **install_kwargs)
    normalized_root = normalize_genome_install(
        genome,
        assembly_name,
        assembly_root,
        force,
        indexers,
        threads,
        cache_blacklist=cache_blacklist,
    )

    tmp_assembly_root = tmp_root / assembly_name
    if tmp_assembly_root.exists():
        shutil.rmtree(tmp_assembly_root)

    return normalized_root


def genomes_main(argv, workflow_root, workflow_config_file, default_site_config):
    args = parse_genomes_arguments(argv)
    site_config_file = Path(args.site_config).expanduser().resolve() if args.site_config else default_site_config
    config = load_site_settings(workflow_root, workflow_config_file, site_config_file)

    genomes_command = {
        "ls": "list",
        "local": "installed",
        "motif-db": "motifs",
    }.get(args.genomes_command, args.genomes_command)

    if genomes_command == "installed":
        rows = list_installed_rows(config["genome_assembly_dir"], getattr(args, "species", None))
        print_installed_rows(rows, args.limit)
        return

    if genomes_command == "blacklist":
        blacklist_bed = resolve_blacklist_bed(
            assembly_name=args.assembly,
            assembly_root=config["genome_assembly_dir"],
            provider=args.provider,
            threads=args.threads,
            force=args.force,
        )
        print(f"Blacklist BED for '{args.assembly}': {blacklist_bed}")
        return

    if genomes_command == "motifs":
        motif_path = meme_motif_database_path(config["genome_assembly_dir"], name=args.name)
        if args.dry_run:
            print(f"MEME motif database cache path: {motif_path}")
            return
        try:
            cached_motif_db = resolve_meme_motif_database(
                config["genome_assembly_dir"],
                force=args.force,
                url=args.url,
                name=args.name,
            )
        except Exception as exc:
            print(f"Could not cache MEME motif database: {exc}", file=sys.stderr)
            sys.exit(1)
        print(f"MEME motif database: {cached_motif_db}")
        return

    genomepy = import_genomepy()
    args.species = normalize_species_name(args.species)

    if genomes_command == "list":
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
            cache_blacklist=not args.skip_blacklist,
        )
        print(f"Installed '{assembly_name}' to {installed_path}")
        cached_blacklist = find_cached_blacklist_bed(config["genome_assembly_dir"], assembly_name)
        if cached_blacklist:
            print(f"Cached blacklist BED: {cached_blacklist}")
        elif not args.skip_blacklist:
            print(f"No blacklist BED was cached for '{assembly_name}'. Use 'omnomnomics genomes blacklist' to backfill one if needed.")
        if not args.skip_motif_db:
            try:
                cached_motif_db = resolve_meme_motif_database(config["genome_assembly_dir"])
                print(f"Cached MEME motif database: {cached_motif_db}")
            except Exception as exc:
                print(f"No MEME motif database was cached: {exc}", file=sys.stderr)
