import os
from pathlib import Path
import shutil
import subprocess
import unittest


ROOT = Path(__file__).parents[1]


class ExplorerPlotControlsTests(unittest.TestCase):
    def run_r(self, script):
        rscript = shutil.which("Rscript")
        self.assertIsNotNone(rscript, "Explorer plot tests require the declared R environment")
        result = subprocess.run([rscript, "--vanilla", "-"], input=script, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=ROOT,
                                env=dict(os.environ), timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_pca_uses_large_circles_and_named_group_colors(self):
        self.run_r('''
source("src/omnomnomics/workflow/R/de_shared/qc_plots.R")
scores <- data.frame(PC1=c(-2,-1,1,2), PC2=c(1,-1,1,-1), condition=c("A","A","B","B"))
colors <- c(A="#123456", B="#ABCDEF")
p <- plot_pca_scores(scores, color_by="condition", color_palette=colors)
points <- ggplot2::ggplot_build(p)$data[[1]]
stopifnot(all(points$size == 6), all(points$shape == 19), identical(points$colour, unname(colors[scores$condition])))
scores$condition <- "A"
one_group <- ggplot2::ggplot_build(plot_pca_scores(scores, color_by="condition", color_palette=colors))$data[[1]]
stopifnot(all(one_group$colour == colors[["A"]]))
names(scores)[3] <- "Cell group"
non_syntactic <- ggplot2::ggplot_build(plot_pca_scores(scores, color_by="Cell group", color_palette=colors))$data[[1]]
stopifnot(all(non_syntactic$colour == colors[["A"]]))
''')

    def test_ora_and_gsea_keep_their_own_statistics(self):
        self.run_r('''
source("src/omnomnomics/workflow/R/de_shared/enrichment_plots.R")
ora <- data.frame(Description=c("Term A","Term B"), p.adjust=c(.001,.01), Count=c(8,12), method="ORA", direction="up")
dotplot <- ggplot2::ggplot_build(plot_cp_top_terms(ora))
stopifnot(length(dotplot$data) == 2, nrow(dotplot$data[[2]]) == 2, all(sort(dotplot$data[[2]]$x) == c(2,3)))
gsea <- data.frame(Description=c("Term A","Term B"), p.adjust=c(.001,.01), NES=c(2.4,-1.8), method="GSEA", direction="ranked")
bars <- ggplot2::ggplot_build(plot_cp_top_terms(gsea))$data[[2]]
stopifnot(nrow(bars) == 2, all(sort(bars$x) == c(-1.8,2.4)))
mixed <- rbind(gsea, transform(gsea, method="ORA"))
failure <- tryCatch(plot_cp_top_terms(mixed), error=identity)
stopifnot(inherits(failure, "error"))
''')

    def test_group_color_overrides_and_gradient_presets(self):
        self.run_r('''
source("src/omnomnomics/workflow/R/de_shared/explorer_colors.R")
source("src/omnomnomics/workflow/R/shiny_app/modules/mod_settings.R")
stopifnot(explorer_color_id("condition", "A-B") != explorer_color_id("condition", "A.B"))
gradient <- explorer_palette(5, "gradient", "#175A81", "#DA845C")
stopifnot(length(gradient) == 5, toupper(gradient[[1]]) == "#175A81", toupper(gradient[[5]]) == "#DA845C")
md <- shiny::reactiveVal(data.frame(sample_id=c("a1","a2","b1","b2"), condition=c("A","A","B","B")))
shiny::testServer(mod_settings_server, args=list(metadata=shiny::reactive(md())), {
  session$setInputs(group_by="condition", preset="coastal", pca_size=5)
  override <- stats::setNames(list("#123456"), explorer_color_id("condition", "A"))
  do.call(session$setInputs, override)
  colors <- colors_for(c("B", "A"), "condition")
  stopifnot(colors[["A"]] == "#123456", colors[["B"]] == "#D77A63")
  stopifnot(identical(colors_for("A", "condition"), colors))
})
''')


if __name__ == "__main__":
    unittest.main()
