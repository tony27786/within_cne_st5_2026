# Supplementary trees, one PDF with two pages: (1) the L1 maximum-likelihood tree of the 160 representative isolates with every isolate labelled, at the display root of Figure 3A; (2) the display tree of the 160 isolates and the 5 H99 sequencing samples rooted on H99, with the H99-ST5 branch drawn shortened and marked as broken.
# Inputs: G1 output (04_tree/ml/L1.treefile, 06_display/tree/display.treefile) and data/sublineage_crosscov.tsv. Outputs: FigureS_trees.pdf and FigureS_trees_p1.png / _p2.png.
# Run: Rscript figures/figS_trees.R

local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
  source(file.path(dir, "fig_common.R"))
})
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
})

p_tree <- in_path(G1_OUT, "04_tree", "ml", "L1.treefile")
p_disp <- in_path(G1_OUT, "06_display", "tree", "display.treefile")
p_sub <- in_path("data", "sublineage_crosscov.tsv")
tr0 <- read.tree(p_tree)
disp <- read.tree(p_disp)
subtab <- read_tsv(p_sub)
sub_members <- strsplit(subtab$members[subtab$group == "sub-lineage"], ",")[[1]]

support_ok <- function(lab) {
  sp <- do.call(rbind, strsplit(ifelse(lab %in% c("", "Root") | is.na(lab), "NA/NA", lab), "/"))
  ok <- suppressWarnings(as.numeric(sp[, 1]) >= 80 & as.numeric(sp[, 2]) >= 95)
  ok[is.na(ok)] <- FALSE
  ok
}

# Rooting (same rule as fig3.R)
h99 <- grep("^CnH99", disp$tip.label, value = TRUE)
dr <- root(disp, outgroup = h99, resolve.root = TRUE, edgelabel = TRUE)
kids <- dr$edge[dr$edge[, 1] == Ntip(dr) + 1, 2]
is_h99_side <- vapply(kids, function(x) all((if (x <= Ntip(dr)) dr$tip.label[x] else extract.clade(dr, x)$tip.label) %in% h99), logical(1))
ing <- extract.clade(dr, kids[!is_h99_side])
sides <- lapply(ing$edge[ing$edge[, 1] == Ntip(ing) + 1, 2],
                function(x) if (x <= Ntip(ing)) ing$tip.label[x] else extract.clade(ing, x)$tip.label)
root_side <- sides[[which.min(lengths(sides))]]
tr <- root_on_split(tr0, root_side)   # fig_common.R; same rule as fig3.R, checked against the 43 / 117 isolates

# Shared drawing
tips_below <- function(t) lapply(seq_len(Ntip(t) + t$Nnode), function(n) if (n <= Ntip(t)) t$tip.label[n] else extract.clade(t, n)$tip.label)

draw_tree <- function(t, label_fun, scale_y) {
  tb <- tips_below(t)
  p <- ggtree(t, ladderize = TRUE)
  p$data$grp <- ifelse(vapply(tb, function(x) all(x %in% sub_members), logical(1))[p$data$node], "sub", "other")
  p$data$supported <- c(rep(FALSE, Ntip(t)), support_ok(t$node.label))[p$data$node]
  p$data$tiplab <- ifelse(p$data$isTip, label_fun(p$data$label), NA)
  xmax <- max(p$data$x)
  p + aes(colour = grp) +
    geom_nodepoint(aes(subset = supported), size = 0.25, colour = COL$ink2) +
    geom_tiplab(aes(label = tiplab), size = 1.5, family = FONT, offset = 0.004 * xmax, na.rm = TRUE) +
    scale_colour_manual(values = c(sub = COL$sublineage, other = COL$ink2), guide = "none") +
    scale_y_continuous(limits = c(scale_y, Ntip(t) + 1), expand = c(0, 0)) +
    theme(text = element_text(family = FONT), plot.margin = margin(3, 3, 3, 3, "mm"),
          plot.tag = element_text(face = "bold", size = 9, family = FONT))   # same tag style as theme_fig()
}
label_main <- function(x) ifelse(x == "Cn13", "Cn13 (mixed DNA)", x)

# Page 1: L1 tree, 160 representative isolates
x1 <- max(node.depth.edgelength(tr)[seq_len(Ntip(tr))])
p1 <- draw_tree(tr, label_main, scale_y = -4) +
  annotate("segment", x = 0, xend = 1e-5, y = -1.5, yend = -1.5, linewidth = 0.35, colour = COL$ink) +
  annotate("text", x = 0, y = -3.2, label = "0.00001 substitutions per site", hjust = 0, size = 1.9, family = FONT, colour = COL$ink) +
  xlim(0, x1 * 1.18) +
  labs(tag = "A")

# Page 2: display tree with H99, broken H99-ST5 branch
dr$node.label[is.na(dr$node.label)] <- ""
st5_node <- kids[!is_h99_side]
e_st5 <- which(dr$edge[, 2] == st5_node)
e_h99 <- which(dr$edge[, 2] == kids[is_h99_side])
true_len <- dr$edge.length[e_st5] + dr$edge.length[e_h99]
st5_depth <- max(node.depth.edgelength(ing)[seq_len(Ntip(ing))])
dshow <- dr
dshow$edge.length[e_h99] <- 0.1 * st5_depth
dshow$edge.length[e_st5] <- 0.6 * st5_depth
brk_x <- 0.1 * st5_depth + 0.3 * st5_depth              # middle of the shortened ST5 stem (root at x = 0)
label_disp <- function(x) ifelse(x == "Cn13", "Cn13 (mixed DNA)", ifelse(x %in% h99, paste0(x, " (H99)"), x))
p2_base <- draw_tree(dshow, label_disp, scale_y = -6)
st5_y <- p2_base$data$y[p2_base$data$node == st5_node]
x2 <- max(p2_base$data$x)
p2 <- p2_base +
  annotate("segment", x = brk_x - 0.02 * x2 + c(-0.006, 0.006) * x2, xend = brk_x + 0.02 * x2 + c(-0.006, 0.006) * x2,
           y = st5_y - 2.5, yend = st5_y + 2.5, linewidth = 0.35, colour = COL$ink) +
  annotate("text", x = brk_x, y = st5_y - 4, hjust = 0.5, vjust = 1, size = 1.9, family = FONT, colour = COL$ink, lineheight = 0.95,
           label = sprintf("H99\u2013ST5 branch shortened\n(%.5f substitutions per site)", true_len)) +
  annotate("segment", x = 0, xend = 1e-5, y = -2.5, yend = -2.5, linewidth = 0.35, colour = COL$ink) +
  annotate("text", x = 0, y = -4.6, label = "0.00001 substitutions per site", hjust = 0, size = 1.9, family = FONT, colour = COL$ink) +
  xlim(0, x2 * 1.16) +
  labs(tag = "B")

# Write: two-page PDF and one PNG per page
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
W_CM <- 17.5; H_CM <- 22
pdf_device(file.path(OUT_DIR, "FigureS_trees.pdf"), width = W_CM / 2.54, height = H_CM / 2.54)   # fig_common.R (Quartz on macOS)
print(p1); print(p2)
invisible(dev.off())
# 600 dpi: the tip labels are small, and the Word document scales the page to the text width
ggsave(file.path(OUT_DIR, "FigureS_trees_p1.png"), p1, width = W_CM, height = H_CM, units = "cm", dpi = 600)
ggsave(file.path(OUT_DIR, "FigureS_trees_p2.png"), p2, width = W_CM, height = H_CM, units = "cm", dpi = 600)
message("Supplementary trees done")
