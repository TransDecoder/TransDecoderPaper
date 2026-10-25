# Shared helpers for the paper notebooks: paths, figure style, labels, number formatting.
# Each notebook starts with  source(file.path("R", "common.R"))  (knitted from notebooks/).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---- paths ------------------------------------------------------------------------------------
REPO <- normalizePath(file.path(getwd(), ".."))
data_path <- function(...) file.path(REPO, "data", ...)
fig_dir <- file.path(REPO, "figures")
read_tab <- function(...) readr::read_tsv(data_path(...), show_col_types = FALSE, comment = "")

# Link to an evidence record, e.g. evidence("E25") -> [E25](../evidence/E25_final_release_reruns.md)
evidence <- function(id) {
  f <- list.files(file.path(REPO, "evidence"), pattern = paste0("^", id, "_.*\\.md$"))
  stopifnot(length(f) == 1)
  sprintf("[%s](../evidence/%s)", id, f)
}

# ---- figure style (Genome Biology: max width 170 mm, 7-8 pt text) ------------------------------
MM <- 1 / 25.4
FULL_W <- 170 * MM
HALF_W <- 82 * MM

theme_td <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          axis.text = element_text(size = base_size - 0.5, colour = "black"),
          legend.key.size = unit(3, "mm"), legend.text = element_text(size = base_size - 0.5),
          legend.background = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = base_size, face = "plain"),
          plot.title = element_text(size = base_size, face = "plain"),
          plot.tag = element_text(size = base_size + 2, face = "bold"))
}
theme_set(theme_td())

# Save a figure as PDF and PNG under figures/<stem>.{pdf,png}; returns the PNG path for knitr.
save_fig <- function(p, stem, width = FULL_W, height = 4) {
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_f <- file.path(fig_dir, paste0(stem, ".pdf"))
  png_f <- file.path(fig_dir, paste0(stem, ".png"))
  ggsave(pdf_f, p, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_f, p, width = width, height = height, units = "in", dpi = 300)
  invisible(png_f)
}

# ---- labels and colours (same as the earlier matplotlib figures) -------------------------------
ORGS <- c("pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii")
ORG_LABEL <- c(pfalciparum = "P. falciparum", tthermophila = "T. thermophila", spombe = "S. pombe",
               athal = "A. thaliana", dmel = "D. melanogaster", mouse = "M. musculus",
               human = "H. sapiens", creinhardtii = "C. reinhardtii")
ORG_GC <- c(pfalciparum = 22.1, tthermophila = 25.9, spombe = 37.9, athal = 42.2, dmel = 50.9,
            mouse = 51.7, human = 52.3, creinhardtii = 68.3)
ORG_COLOR <- setNames(c("#a6761d", "#e7298a", "#666666", "#66a61e", "#e6ab02", "#7570b3",
                        "#1b9e77", "#d95f02"), ORGS)
ASM_COLOR <- c(lraa = "#1b9e77", stringtie = "#d95f02", trinity = "#7570b3")
ASM_LABEL <- c(lraa = "LRAA (long reads)", stringtie = "StringTie", trinity = "Trinity (genome-guided)")
COMPLETE_C <- "#2c7fb8"
FRAG_C <- "#f03b20"
TD_C <- "#08519c"
GREY <- "#999999"

org_factor <- function(x) factor(x, levels = ORGS, labels = ORG_LABEL[ORGS])

# ---- numbers in text ---------------------------------------------------------------------------
f3 <- function(x) formatC(x, format = "f", digits = 3)
f1 <- function(x) formatC(x, format = "f", digits = 1)
pct <- function(x, digits = 1) paste0(formatC(100 * x, format = "f", digits = digits), "%")
num <- function(x) formatC(x, format = "d", big.mark = ",")

# Check a value quoted in the manuscript against the value computed here; mismatches are listed in
# the notebook's "Checks" section (check_report()) instead of stopping the knit.
.checks <- new.env()
.checks$rows <- list()
check_text <- function(label, computed, quoted, tol = 0) {
  ok <- isTRUE(abs(as.numeric(computed) - as.numeric(quoted)) <= tol + 1e-9)
  .checks$rows[[length(.checks$rows) + 1]] <- tibble(claim = label, computed = as.character(computed),
                                                     manuscript = as.character(quoted), ok = ok)
  invisible(ok)
}
check_report <- function() {
  if (length(.checks$rows) == 0) return(knitr::asis_output("No quoted values checked.\n"))
  knitr::kable(bind_rows(.checks$rows))
}
