# Shared settings for the manuscript figures: paths, theme, palette, tree-rooting helper and export.
# Figure scripts write PDF + PNG (FIG_FORMATS, e.g. "pdf,png,eps,tiff") into FIG_OUTPUT_DIR (default figures/out).

# Labels contain non-ASCII characters (\u escapes); Rscript may start in the C locale.
if (!l10n_info()[["UTF-8"]]) invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

FIG_DIR <- get_script_dir()
PROJ <- normalizePath(file.path(FIG_DIR, ".."), mustWork = TRUE)
OUT_DIR <- Sys.getenv("FIG_OUTPUT_DIR", file.path(FIG_DIR, "out"))
FORMATS <- trimws(strsplit(Sys.getenv("FIG_FORMATS", "pdf,png"), ",")[[1]])

# Default locations of the HPC outputs, relative to the project root
G0_OUT <- file.path("results", "g0")
G1_OUT <- file.path("results", "g1")
MISSCHECK_DIR <- file.path("results", "misscheck")

in_path <- function(...) {
  p <- file.path(PROJ, ...)
  if (!file.exists(p)) stop("missing input: ", p, call. = FALSE)
  p
}
read_tsv <- function(p, ...) read.delim(p, check.names = FALSE, stringsAsFactors = FALSE, ...)

FONT <- "Arial"

theme_fig <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "#52514e"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#52514e"),
      axis.text = element_text(colour = "#52514e"),
      axis.title = element_text(colour = "#0b0b0b"),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "#ecebe7"),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size),
      plot.tag = element_text(face = "bold", size = base_size + 2, family = FONT),
      plot.margin = margin(2, 4, 2, 2, "mm")
    )
}

COL <- list(
  ink = "#0b0b0b",
  ink2 = "#52514e",
  muted = "#9b9a96",
  band = "#ecebe7",
  within = "#2a78d6",
  between = "#eb6834",
  sublineage = "#eb6834",
  interval = c("#86b6ef", "#3987e5", "#1c5cab", "#0d366b"),  # same day, 1-7, 8-30, >30 days
  specimen = c("CSF" = "#1baf7a", "Blood" = "#eda100", "Other" = "#6250d6"),
  hiv = c("Positive" = "#d9d8d4", "Negative" = "#0b0b0b"),
  multi = c("Yes" = "#0b0b0b", "No" = "#f4f3f0")
)

# Distance axis: pseudo-log (linear near 0, log above ~1) so that 0 stays on the axis.
SNP_BREAKS <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
snp_trans <- scales::pseudo_log_trans(sigma = 1, base = 10)
scale_x_snp <- function(limits = c(0, 1400), ...) {
  scale_x_continuous(trans = snp_trans, breaks = SNP_BREAKS, labels = scales::label_comma(),
                     limits = limits, expand = expansion(mult = c(0.015, 0.01)), ...)
}

# Vertical offsets for points sharing a row and a value, so that tied distances stay readable.
stack_offset <- function(row, value, step = 0.14, cap = 0.4) {
  key <- paste(row, value)
  off <- ave(seq_along(key), key, FUN = function(i) (seq_along(i) - (length(i) + 1) / 2) * step)
  pmax(-cap, pmin(cap, off))
}

# Sampling dates of the isolates with the clinical errata applied (a correction is applied only if its archived value matches); the archived dates and the files read are kept as attributes.
isolate_dates <- function() {
  p_iso <- in_path("data", "isolates.csv")
  p_err <- in_path("data", "errata_clinical.tsv")
  iso <- read.csv(p_iso, colClasses = "character", encoding = "UTF-8", check.names = FALSE)
  names(iso)[1] <- sub("^﻿", "", names(iso)[1])
  archived <- setNames(iso$sampling_date, iso$sample)
  d <- archived
  er <- read.delim(p_err, colClasses = "character", na.strings = character())
  er <- er[er$table == "isolates" & er$field == "sampling_date", ]
  for (k in seq_len(nrow(er))) {
    stopifnot(identical(unname(d[er$key_value[k]]), er$archived_value[k]))
    d[[er$key_value[k]]] <- er$corrected_value[k]
  }
  attr(d, "archived") <- archived
  attr(d, "files") <- c(p_iso, p_err)
  d
}

# Root an unrooted tree on the edge that separates the tips in `side` from all other tips, with the root edge split evenly for display.
# ape::root() moves to a neighbouring edge when `side` straddles the arbitrary basal node, so the tree is first rooted inside the complement and the result is checked against `side`.
root_on_split <- function(tr0, side) {
  other <- setdiff(tr0$tip.label, side)
  t1 <- ape::root(tr0, outgroup = other[1], resolve.root = TRUE, edgelabel = TRUE)
  tr <- ape::root(t1, outgroup = side, resolve.root = TRUE, edgelabel = TRUE)
  re <- which(tr$edge[, 1] == ape::Ntip(tr) + 1)
  kids <- lapply(tr$edge[re, 2], function(x) if (x <= ape::Ntip(tr)) tr$tip.label[x] else ape::extract.clade(tr, x)$tip.label)
  stopifnot("the root separates exactly the tips of the display split" = length(re) == 2 && any(vapply(kids, function(k) setequal(k, side), logical(1))) &&
              any(vapply(kids, function(k) setequal(k, other), logical(1))),
            "rooting changes neither the unrooted topology nor the total branch length" = isTRUE(all.equal(ape::unroot(tr), ape::unroot(tr0), use.edge.length = FALSE)) &&
              isTRUE(all.equal(sum(tr$edge.length), sum(tr0$edge.length))))
  tr$edge.length[re] <- sum(tr$edge.length[re]) / 2
  tr
}

# PDF device (FIG_PDF_DEVICE = "quartz" or "cairo"): Quartz on macOS lays text out with CoreText and embeds Arial; cairo_pdf elsewhere.
PDF_DEVICE <- Sys.getenv("FIG_PDF_DEVICE", if (capabilities("aqua")) "quartz" else "cairo")
quartz_pdf <- function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height, ...)
pdf_device <- if (PDF_DEVICE == "quartz") quartz_pdf else function(filename, width, height, ...) grDevices::cairo_pdf(filename, width, height, onefile = TRUE, ...)

save_fig <- function(p, name, width_cm, height_cm) {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  for (f in FORMATS) {
    out <- file.path(OUT_DIR, paste0(name, ".", f))
    args <- list(filename = out, plot = p, width = width_cm, height = height_cm, units = "cm")
    switch(f,
      pdf = do.call(ggsave, c(args, device = pdf_device)),
      png = do.call(ggsave, c(args, device = "png", dpi = 300)),
      eps = do.call(ggsave, c(args, device = cairo_ps, fallback_resolution = 600)),
      tiff = do.call(ggsave, c(args, device = "tiff", dpi = 600, compression = "lzw")),
      stop("unknown format: ", f, call. = FALSE))
    message("wrote ", out)
  }
}
