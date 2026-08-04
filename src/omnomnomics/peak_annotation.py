from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path


GeneKey = tuple[str, str, str]
GeneSpan = tuple[int, int]
PromoterInterval = tuple[str, int, int]


def parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in str(attr_string).strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def _update_span(spans: dict[GeneKey, GeneSpan], key: GeneKey, start: int, end: int) -> None:
    previous = spans.get(key)
    if previous is None:
        spans[key] = (start, end)
    else:
        spans[key] = (min(previous[0], start), max(previous[1], end))


def _read_gene_models(
    gtf_file: str | Path,
    exon_output: str | Path | None = None,
) -> tuple[
    dict[GeneKey, GeneSpan],
    dict[GeneKey, GeneSpan],
    dict[GeneKey, GeneSpan],
    dict[GeneKey, str],
    dict[GeneKey, set[int]],
]:
    native_gene_spans: dict[GeneKey, GeneSpan] = {}
    transcript_spans: dict[GeneKey, GeneSpan] = {}
    exon_spans: dict[GeneKey, GeneSpan] = {}
    gene_labels: dict[GeneKey, str] = {}
    transcript_tss: dict[GeneKey, set[int]] = {}

    exon_handle = open(exon_output, "w", encoding="utf-8") if exon_output else None
    try:
        with open(gtf_file, encoding="utf-8") as gtf_handle:
            for line in gtf_handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _, feature, start_s, end_s, _, strand, _, attrs_raw = fields
                if strand not in {"+", "-"} or feature not in {"gene", "transcript", "mRNA", "exon"}:
                    continue
                try:
                    start = int(start_s)
                    end = int(end_s)
                except ValueError:
                    continue
                if end < start:
                    continue

                attrs = parse_gtf_attributes(attrs_raw)
                gene_id = attrs.get("gene_id") or attrs.get("gene_name")
                if not gene_id:
                    continue
                gene_name = attrs.get("gene_name") or gene_id
                key = (chrom, gene_id, strand)
                gene_label = f"{gene_name}|{gene_id}"
                gene_labels.setdefault(key, gene_label)

                bed_start = max(start - 1, 0)
                bed_end = end
                if feature == "gene":
                    _update_span(native_gene_spans, key, bed_start, bed_end)
                elif feature in {"transcript", "mRNA"}:
                    _update_span(transcript_spans, key, bed_start, bed_end)
                    tss = bed_start if strand == "+" else bed_end - 1
                    transcript_tss.setdefault(key, set()).add(tss)
                elif feature == "exon":
                    _update_span(exon_spans, key, bed_start, bed_end)
                    if exon_handle is not None:
                        exon_handle.write(
                            f"{chrom}\t{bed_start}\t{bed_end}\t{gene_label}\t0\t{strand}\n"
                        )
    finally:
        if exon_handle is not None:
            exon_handle.close()

    return native_gene_spans, transcript_spans, exon_spans, gene_labels, transcript_tss


def _resolved_gene_spans(
    native_gene_spans: dict[GeneKey, GeneSpan],
    transcript_spans: dict[GeneKey, GeneSpan],
    exon_spans: dict[GeneKey, GeneSpan],
) -> dict[GeneKey, GeneSpan]:
    resolved: dict[GeneKey, GeneSpan] = {}
    for key in set(native_gene_spans) | set(transcript_spans) | set(exon_spans):
        resolved[key] = (
            native_gene_spans.get(key)
            or transcript_spans.get(key)
            or exon_spans[key]
        )
    return resolved


def _tss_positions_for_gene(
    key: GeneKey,
    gene_span: GeneSpan,
    transcript_tss: dict[GeneKey, set[int]],
) -> Iterable[int]:
    positions = transcript_tss.get(key)
    if positions:
        return sorted(positions)
    return [gene_span[0] if key[2] == "+" else gene_span[1] - 1]


def _promoter_bounds(
    tss: int,
    strand: str,
    promoter_upstream: int,
    promoter_downstream: int,
) -> tuple[int, int]:
    if strand == "+":
        return max(tss - promoter_upstream, 0), tss + promoter_downstream
    boundary = tss + 1
    return max(boundary - promoter_downstream, 0), boundary + promoter_upstream


def build_gtf_annotation_sources(
    gtf_file: str | Path,
    work_dir: str | Path,
    promoter_upstream: int = 1000,
    promoter_downstream: int = 1000,
) -> dict[str, str]:
    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)
    genes_path = work_path / "genes.unsorted.bed"
    exons_path = work_path / "exons.unsorted.bed"
    promoters_path = work_path / "promoters.unsorted.bed"
    tss_path = work_path / "tss.unsorted.bed"

    native, transcripts, exons, labels, transcript_tss = _read_gene_models(
        gtf_file,
        exon_output=exons_path,
    )
    gene_spans = _resolved_gene_spans(native, transcripts, exons)

    with (
        open(genes_path, "w", encoding="utf-8") as genes_out,
        open(promoters_path, "w", encoding="utf-8") as promoters_out,
        open(tss_path, "w", encoding="utf-8") as tss_out,
    ):
        for key in sorted(gene_spans):
            chrom, _, strand = key
            start, end = gene_spans[key]
            label = labels[key]
            genes_out.write(f"{chrom}\t{start}\t{end}\t{label}\t0\t{strand}\n")
            for tss in _tss_positions_for_gene(key, gene_spans[key], transcript_tss):
                prom_start, prom_end = _promoter_bounds(
                    tss,
                    strand,
                    promoter_upstream,
                    promoter_downstream,
                )
                if prom_end > prom_start:
                    promoters_out.write(
                        f"{chrom}\t{prom_start}\t{prom_end}\t{label}\t0\t{strand}\n"
                    )
                tss_out.write(f"{chrom}\t{tss}\t{tss + 1}\t{label}\t0\t{strand}\n")

    return {
        "genes": str(genes_path),
        "exons": str(exons_path),
        "promoters": str(promoters_path),
        "tss": str(tss_path),
    }


def promoter_intervals_by_gene(
    gtf_file: str | Path,
    promoter_upstream: int = 1000,
    promoter_downstream: int = 1000,
) -> dict[str, list[PromoterInterval]]:
    native, transcripts, exons, labels, transcript_tss = _read_gene_models(gtf_file)
    gene_spans = _resolved_gene_spans(native, transcripts, exons)
    promoter_map: dict[str, list[PromoterInterval]] = {}
    for key in sorted(gene_spans):
        chrom, _, strand = key
        label = labels[key]
        intervals = set(promoter_map.get(label, []))
        for tss in _tss_positions_for_gene(key, gene_spans[key], transcript_tss):
            prom_start, prom_end = _promoter_bounds(
                tss,
                strand,
                promoter_upstream,
                promoter_downstream,
            )
            if prom_end > prom_start:
                intervals.add((chrom, prom_start, prom_end))
        if intervals:
            promoter_map[label] = sorted(intervals)
    return promoter_map
