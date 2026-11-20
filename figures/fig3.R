# Figure 3: population structure of the 160 representative isolates. A: L1 maximum-likelihood tree rooted for display at the split where H99 attaches (43 / 117 isolates), with columns for HIV status, sampling year, specimen and patients with more than one sequenced isolate. B: PCA (PC1 vs PC2) coloured by the 15-isolate clade defined on the tree.
# Inputs: G1 output (04_tree/ml/L1.treefile, 06_display/tree/display.treefile, 03_l1/pca_scores.tsv, 03_l1/pca_variance.tsv), data/strain_metadata.tsv, data/strain_hiv.tsv, data/sublineage_crosscov.tsv; isolates table and errata for the sampling year (fig_common.R).
# Run: Rscript figures/fig3.R

local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
  source(file.path(dir, "fig_common.R"))
})
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggnewscale)
  library(ggrepel)
})

p_tree <- in_path(G1_OUT, "04_tree", "ml", "L1.treefile")
p_disp <- in_path(G1_OUT, "06_display", "tree", "display.treefile")
p_pca <- in_path(G1_OUT, "03_l1", "pca_scores.tsv")
p_var <- in_path(G1_OUT, "03_l1", "pca_variance.tsv")
p_meta <- in_path("data", "strain_metadata.tsv")
p_hiv <- in_path("data", "strain_hiv.tsv")
p_sub <- in_path("data", "sublineage_crosscov.tsv")

tr0 <- read.tree(p_tree)
disp <- read.tree(p_disp)
pca <- read_tsv(p_pca)
pvar <- read_tsv(p_var)
meta <- read_tsv(p_meta)
hiv <- read_tsv(p_hiv)
subtab <- read_tsv(p_sub)

# Display root: the split where the H99 clade attaches in the display tree
h99 <- grep("^CnH99", disp$tip.label, value = TRUE)
dr <- root(disp, outgroup = h99, resolve.root = TRUE)
kids <- dr$edge[dr$edge[, 1] == Ntip(dr) + 1, 2]
ingroup <- kids[vapply(kids, function(x) x > Ntip(dr) && !all(extract.clade(dr, x)$tip.label %in% h99), logical(1))]
ing <- extract.clade(dr, ingroup)
sides <- lapply(ing$edge[ing$edge[, 1] == Ntip(ing) + 1, 2],
                function(x) if (x <= Ntip(ing)) ing$tip.label[x] else extract.clade(ing, x)$tip.label)
root_side <- sides[[which.min(lengths(sides))]]

tr <- root_on_split(tr0, root_side)   # fig_common.R; checks that the root separates exactly these isolates

# Tip annotation
rep <- hiv[hiv$wgs_role == "representative", ]
multi_pat <- unique(hiv$patient_id[hiv$wgs_role == "within_patient_extra"])
md <- merge(rep[, c("sample", "patient_id", "hiv_status_clean")],
            meta[, c("sample", "sampling_date", "specimen")], by = "sample", all.x = TRUE)
md$HIV <- factor(c(positive = "Positive", negative = "Negative")[md$hiv_status_clean], levels = names(COL$hiv))
# Sampling year from the isolate dates with the clinical errata applied (fig_common.R).
idates <- isolate_dates()
md$Year <- as.integer(substr(unname(idates[md$sample]), 1, 4))
md$Specimen <- factor(ifelse(md$specimen == "Cerebrospinal fluid", "CSF", ifelse(md$specimen == "Whole blood", "Blood", "Other")),
                      levels = names(COL$specimen))
md$Multi <- factor(ifelse(md$patient_id %in% multi_pat, "Yes", "No"), levels = names(COL$multi))
sub_members <- strsplit(subtab$members[subtab$group == "sub-lineage"], ",")[[1]]

sup <- do.call(rbind, strsplit(ifelse(tr$node.label %in% c("", "Root") | is.na(tr$node.label), "NA/NA", tr$node.label), "/"))
supported <- suppressWarnings(as.numeric(sup[, 1]) >= 80 & as.numeric(sup[, 2]) >= 95)
supported[is.na(supported)] <- FALSE

# Panel A: tree with annotation columns
# Rectangular layout (a fan layout compresses the short internal branches); branch lengths are not transformed.
xmax <- max(node.depth.edgelength(tr)[seq_len(Ntip(tr))])
col_w <- 0.045                        # column width as a fraction of the tree depth (gheatmap 'width')
gap <- 0.008 * xmax
col_off <- function(k) 0.02 * xmax + k * (col_w * xmax + gap)
heat <- data.frame(row.names = md$sample, HIV = md$HIV, Year = md$Year, Specimen = md$Specimen, Multi = md$Multi)
col_names <- c("HIV", "Year", "Specimen", ">1 isolate")
sub_node <- getMRCA(tr, sub_members)

# Branch colour: the 15-isolate clade (its stem edge and all descendants) in the PCA colour, other branches dark grey.
# Branch width: edges above clades of >= 10 isolates are drawn thicker so that the backbone can be followed.
tips_below <- lapply(seq_len(Ntip(tr) + tr$Nnode), function(n) if (n <= Ntip(tr)) tr$tip.label[n] else extract.clade(tr, n)$tip.label)
n_desc <- lengths(tips_below)
in_sub <- vapply(tips_below, function(t) all(t %in% sub_members), logical(1))   # sub_node and everything below it
pt <- ggtree(tr, ladderize = TRUE)
pt$data$grp <- ifelse(in_sub[pt$data$node], "sub", "other")
pt$data$w <- ifelse(n_desc[pt$data$node] >= 10, 0.45, 0.22)
pt$data$supported <- c(rep(FALSE, Ntip(tr)), supported)[pt$data$node]
pt$data$marked <- pt$data$label %in% c("Cn13", "Cn14", "Cn101", "Cn103")
tip_y <- setNames(pt$data$y[pt$data$isTip], pt$data$label[pt$data$isTip])
pt <- pt + aes(colour = grp, linewidth = w) +
  scale_colour_manual(values = c(sub = COL$sublineage, other = COL$ink2), guide = "none") +
  scale_linewidth_identity() +
  geom_nodepoint(aes(subset = supported), size = 0.3, colour = COL$ink2, inherit.aes = TRUE) +
  geom_tippoint(aes(subset = marked), size = 0.6, colour = COL$ink)

add_col <- function(p, k, fill_scale) {
  gheatmap(p, heat[, k, drop = FALSE], offset = col_off(k - 1), width = col_w, colnames = FALSE, color = NA) +
    fill_scale + new_scale_fill()
}
pt <- add_col(pt, 1, scale_fill_manual(values = COL$hiv, breaks = names(COL$hiv), name = "HIV status",
                                       guide = guide_legend(order = 1)))
# Colour bar drawn as a vector gradient (the Quartz PDF device flips raster images vertically).
pt <- add_col(pt, 2, scale_fill_gradient(low = "#cde2fb", high = "#0d366b", breaks = c(2015, 2020, 2025), name = "Sampling year",
                                         guide = guide_colourbar(order = 2, barheight = unit(10, "mm"), barwidth = unit(2.5, "mm"), display = "gradient")))
pt <- add_col(pt, 3, scale_fill_manual(values = COL$specimen, breaks = names(COL$specimen), name = "Specimen",
                                       guide = guide_legend(order = 3)))
pt <- add_col(pt, 4, scale_fill_manual(values = COL$multi, breaks = names(COL$multi), name = "Patient with >1\nsequenced isolate",
                                       guide = guide_legend(order = 4)))

col_mid <- xmax + col_off(0:3) + col_w * xmax / 2
lab_x <- xmax + col_off(4) + 0.04 * xmax
pair_lab <- data.frame(
  y = c(mean(tip_y[c("Cn13", "Cn14")]), mean(tip_y[c("Cn101", "Cn103")])),
  label = c("Cn13 (mixed DNA) / Cn14", "Cn101 / Cn103")
)
sub_y <- range(tip_y[sub_members])
n_tip <- Ntip(tr)

pA <- pt +
  annotate("text", x = col_mid, y = n_tip + 1.5, label = col_names, angle = 90, hjust = 0, size = 2, family = FONT, colour = COL$ink) +
  geom_text(data = pair_lab, aes(x = lab_x, y = y, label = label), hjust = 0, size = 2, family = FONT, colour = COL$ink,
            inherit.aes = FALSE) +
  annotate("segment", x = lab_x, xend = lab_x, y = sub_y[1] - 0.4, yend = sub_y[2] + 0.4, linewidth = 0.8, colour = COL$sublineage) +
  annotate("text", x = lab_x + 0.015 * xmax, y = mean(sub_y), label = sprintf("%d-isolate clade", length(sub_members)),
           hjust = 0, size = 2, family = FONT, colour = COL$ink) +
  annotate("segment", x = 0, xend = 1e-5, y = -2, yend = -2, linewidth = 0.35, colour = COL$ink) +
  annotate("text", x = 0, y = -4.5, label = "0.00001 substitutions per site", hjust = 0, size = 2,
           family = FONT, colour = COL$ink) +
  scale_y_continuous(limits = c(-6, n_tip + 22), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, lab_x + 0.62 * xmax), expand = c(0, 0)) +
  labs(tag = "A") +
  theme(text = element_text(family = FONT, size = 7),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        plot.tag = element_text(face = "bold", size = 9), plot.margin = margin(1, 1, 1, 1, "mm"))

# Panel B: PCA
ve <- setNames(pvar$explained_variance, pvar$PC)
lab_sub <- sprintf("%d-isolate clade (defined on the tree)", length(sub_members))
lab_oth <- sprintf("Other isolates (%d)", nrow(pca) - length(sub_members))
pca$group <- factor(ifelse(pca$sample %in% sub_members, lab_sub, lab_oth), levels = c(lab_oth, lab_sub))
pca <- pca[order(pca$group), ]
lab_pca <- pca[pca$sample == "Cn13", ]

pB <- ggplot(pca, aes(PC1, PC2, colour = group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = COL$band) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = COL$band) +
  geom_point(size = 1.2, alpha = 0.85, stroke = 0) +
  geom_text_repel(data = lab_pca, aes(label = "Cn13 (mixed DNA)"), colour = COL$ink, size = 2, family = FONT,
                  min.segment.length = 0, segment.size = 0.25, nudge_x = 15, nudge_y = 12, seed = 1) +
  annotate("text", x = max(pca$PC1), y = max(pca$PC2) + 18, label = lab_sub, hjust = 1, size = 2, family = FONT, colour = COL$ink) +
  annotate("segment", x = 2, xend = median(pca$PC1[pca$group == lab_oth]), y = -78, yend = -4, linewidth = 0.25, colour = COL$ink2) +
  annotate("text", x = 3, y = -85, label = lab_oth, hjust = 0, size = 2, family = FONT, colour = COL$ink) +
  scale_colour_manual(values = setNames(c(COL$muted, COL$sublineage), c(lab_oth, lab_sub)), name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[["PC1"]]), y = sprintf("PC2 (%.1f%%)", 100 * ve[["PC2"]]), tag = "B") +
  theme_fig() +
  theme(panel.grid.major.x = element_blank(), aspect.ratio = 1.1, legend.position = "none")

# Assemble
# Legends of A go to the free area under B; B is labelled directly.
fig <- pA + pB + guide_area() +
  plot_layout(design = "AB\nAC", widths = c(2, 1), heights = c(1, 1), guides = "collect") &
  theme(text = element_text(family = FONT), legend.justification = c(0, 1), legend.spacing.y = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0))
save_fig(fig, "Figure3", width_cm = 17.5, height_cm = 16)
message("Figure 3 done")
