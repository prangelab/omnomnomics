from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RULE = ROOT / "src" / "omnomnomics" / "workflow" / "rules" / "12.call_DE.smk"
TEMPLATE = ROOT / "src" / "omnomnomics" / "workflow" / "templates" / "de_core.R.tmpl"


def test_rna_de_renderer_passes_genome_version():
    rule_source = RULE.read_text()
    template_source = TEMPLATE.read_text()

    assert '"__GENOME_VERSION__"' in rule_source
    assert "genome_version <- __GENOME_VERSION__" in template_source


def test_rna_de_uses_genome_context_for_species_detection():
    source = TEMPLATE.read_text()

    assert "detect_species_from_context <- function" in source
    assert 'grepl("^(grch|hg|human|homo)"' in source
    assert 'grepl("^(grcm|mm|mouse|mus)"' in source
    assert "detect_species_from_context(rownames(dds), genome_version)" in source


def test_rna_de_preserves_symbol_like_gene_ids():
    source = TEMPLATE.read_text()

    assert "is_symbol_like <- !is_ensembl_gene & !is_other_accession" in source
    assert "symbols[is_symbol_like] <- gene_ids[is_symbol_like]" in source


def test_de_templates_use_version_compatible_msigdbr_arguments():
    chrom_template = ROOT / "src" / "omnomnomics" / "workflow" / "templates" / "de_core_chrom.R.tmpl"

    for source in (TEMPLATE.read_text(), chrom_template.read_text()):
        assert "load_msigdb_collection <- function" in source
        assert '"collection" %in% supported_args' in source
        assert "call_args$subcollection" in source


def test_de_templates_support_apeglm_factor_contrasts():
    chrom_template = ROOT / "src" / "omnomnomics" / "workflow" / "templates" / "de_core_chrom.R.tmpl"

    for source in (TEMPLATE.read_text(), chrom_template.read_text()):
        assert "resolve_factor_coefficient <- function" in source
        assert "shrink_factor_contrast <- function" in source
        assert "direction * coefficient_res$log2FoldChange" in source
