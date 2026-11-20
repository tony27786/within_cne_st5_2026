# Figure 2: sampling timeline of the patients with more than one isolate (A), technical noise, intra-patient and inter-patient SNP distances (B, C), plus the supplementary plot of all 40 intra-patient pairs (FigureS_within_pairs).
# Inputs: G1 02_l0/main/{ruler_A_pairs, ruler_B_pairs, ruler_B_by_patient, pairwise_L0, anchors}.tsv; isolates table and errata (sampling dates with the errata applied, via fig_common.R).
# Run: Rscript figures/fig2.R

local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
  source(file.path(dir, "fig_common.R"))
})

L0 <- file.path(G1_OUT, "02_l0", "main")
p_A <- in_path(L0, "ruler_A_pairs.tsv")
p_B <- in_path(L0, "ruler_B_pairs.tsv")
p_Bp <- in_path(L0, "ruler_B_by_patient.tsv")
p_L0 <- in_path(L0, "pairwise_L0.tsv")
p_anc <- in_path(L0, "anchors.tsv")

ra_all <- read_tsv(p_A)
ra <- ra_all[grepl("^CnH99", ra_all$sample_a) & grepl("^CnH99", ra_all$sample_b), ]  # file also holds Cn169 and Cn169b
rb <- read_tsv(p_B)
# Intervals from the isolate dates with the clinical errata applied.
idates <- isolate_dates()
dd <- as.integer(abs(as.Date(unname(idates[rb$isolate_2])) - as.Date(unname(idates[rb$isolate_1]))))
rb$delta_days <- dd
rbp <- read_tsv(p_Bp)
l0 <- read_tsv(p_L0)
anc <- read_tsv(p_anc)
anchor <- setNames(anc$value, anc$anchor)

rc <- l0[l0$pair_type == "rep_rep", ]
max_b <- max(rb$D_raw)
zero_c <- rc[rc$D_raw == 0, ]
n_le <- sum(rc$D_raw <= max_b)

# Dotted: 95th percentile of the H99 replicate pairs; dashed: the 5-SNP reference value (the same-day different-specimen maximum of 2, rounded up).
ref_lines <- list(
  geom_vline(xintercept = as.numeric(anchor[["t_tech"]]), linetype = "dotted", linewidth = 0.35, colour = COL$ink2),
  geom_vline(xintercept = as.numeric(anchor[["t_event"]]), linetype = "dashed", linewidth = 0.35, colour = COL$ink2)
)
spec_short <- function(x) {
  out <- c("Cerebrospinal fluid" = "CSF", "Whole blood" = "blood", "Drainage fluid" = "drainage fluid", "Sputum" = "sputum")[x]
  ifelse(is.na(out), tolower(x), out)
}

# Figure 2A: sampling timeline of the patients with more than one isolate
exc <- rb[rb$D_raw > as.numeric(anchor[["t_event"]]), ]   # the three within-patient pairs > 5 SNPs (panels A and B)
exc <- exc[order(exc$D_raw), ]
# Patient order of panel A and the supplementary pair plot: by largest within-patient distance, then number of pairs.
pat_order <- rbp[order(rbp$D_raw_max, rbp$n_pairs, rbp$patient_id), ]

# One row per isolate: days since the patient's first isolate, specimen group, and whether it comes from a later admission than the patient's first isolate.
timeline_data <- function(rb, dates) {
  ti <- unique(rbind(
    data.frame(patient = rb$patient_id, isolate = rb$isolate_1, admission = rb$admission_1, specimen = rb$specimen_1),
    data.frame(patient = rb$patient_id, isolate = rb$isolate_2, admission = rb$admission_2, specimen = rb$specimen_2)
  ))
  ti$date <- as.Date(unname(dates[ti$isolate]))
  ti <- ti[order(ti$patient, ti$date, ti$isolate), ]
  first <- match(ti$patient, ti$patient)   # row of each patient's first isolate
  ti$day <- as.integer(ti$date - ti$date[first])
  ti$later <- ti$admission != ti$admission[first]
  ti$group <- ifelse(ti$specimen == "Cerebrospinal fluid", "CSF", ifelse(ti$specimen == "Whole blood", "Blood", "Other"))
  ti
}

# Rows run top to bottom in patient order; same-day isolates are stacked in a taller row; square-root day scale; the later isolate of each pair > 5 SNPs carries the distance.
timeline_plot <- function(ti, pat_order, exc) {
  same_day <- tapply(duplicated(paste(ti$patient, ti$day)), ti$patient, any)[pat_order]
  h <- ifelse(same_day, 2, 1)
  mid <- setNames(sum(h) - (cumsum(h) - h) - h / 2, pat_order)
  ti$y <- mid[ti$patient] - stack_offset(ti$patient, ti$day, step = 0.8, cap = 0.4)
  span <- aggregate(day ~ patient, ti, max)
  span$y <- mid[span$patient]
  span$exc <- span$patient %in% exc$patient_id
  last <- ti$patient %in% exc$patient_id & ti$day == span$day[match(ti$patient, span$patient)]
  ti$label <- ifelse(last, sprintf("%s (%d SNPs)", ti$isolate, exc$D_raw[match(ti$patient, exc$patient_id)]), ti$isolate)
  ti$adm <- factor(ifelse(ti$later, "later", "same"), levels = c("same", "later"))
  ggplot(ti, aes(day, y)) +
    geom_segment(data = span[!span$exc, ], aes(x = 0, xend = day, y = y, yend = y), linewidth = 0.3, colour = COL$muted) +
    geom_segment(data = span[span$exc, ], aes(x = 0, xend = day, y = y, yend = y), linewidth = 0.45, colour = COL$ink2) +
    # White label boxes hide the line; points are drawn last so that close neighbours are not covered.
    geom_label(aes(label = label), hjust = 0, nudge_x = 0.1, size = 1.9, family = FONT, colour = COL$ink, fill = "white",
               border.colour = NA, label.padding = unit(0.15, "mm"), label.r = unit(0, "mm")) +
    geom_point(aes(colour = group, shape = adm), size = 1.5, stroke = 0.6, fill = "white") +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 2, 5, 10, 20, 30, 40, 50), limits = c(0, 70),
                       expand = expansion(mult = c(0.015, 0.01))) +
    scale_y_continuous(breaks = unname(mid), labels = names(mid), limits = c(0, sum(h)), expand = c(0, 0)) +
    scale_colour_manual(values = COL$specimen, breaks = names(COL$specimen), name = "Specimen",
                        labels = c("CSF", "Blood", "Other (sputum, pleural effusion, drainage fluid)")) +
    scale_shape_manual(values = c(same = 16, later = 21), name = "Admission",
                       labels = c("Same as the first isolate", "Later admission")) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 1.8)),
           shape = guide_legend(order = 2, override.aes = list(size = 1.6, colour = COL$ink2))) +
    labs(x = "Days since the patient's first isolate", y = NULL, tag = "A") +
    theme_fig() +
    theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "#f4f3f0"),
          axis.text.y = element_text(size = 6, colour = COL$ink),
          legend.position = "inside", legend.position.inside = c(0.995, 0.995), legend.justification = c(1, 1),
          legend.box = "horizontal", legend.box.just = "top",
          legend.background = element_rect(fill = "white", colour = NA))
}

iso_tab <- read.csv(attr(idates, "files")[1], colClasses = "character", encoding = "UTF-8", check.names = FALSE)
names(iso_tab)[1] <- sub("^﻿", "", names(iso_tab)[1])
ti <- timeline_data(rb, idates)
ti$role <- iso_tab$wgs_role[match(ti$isolate, iso_tab$sample)]
pA <- timeline_plot(ti, pat_order$patient_id, exc)

# Figure 2B: noise floor, the patients with all pairs <= 5 SNPs, and the three exceptions
low <- rb[!rb$patient_id %in% exc$patient_id, ]

row_h <- c(1, 2.5, rep(1.05, nrow(exc)))          # H99, the 23 patients (stacked dots), one row per exception
row_top <- sum(row_h) - (cumsum(row_h) - row_h)   # rows run top to bottom
row_mid <- row_top - row_h / 2
lab_rows <- c(sprintf("H99 pairs\n(5 sequencing samples, %d pairs)", nrow(ra)),
              sprintf("%d patients\n(%d pairs, all \u2264 5)", length(unique(low$patient_id)), nrow(low)),
              exc$patient_id)
da <- rbind(
  data.frame(row = 1, d = ra$D_raw, grp = "h99"),
  data.frame(row = 2, d = low$D_raw, grp = "within"),
  data.frame(row = 2 + seq_len(nrow(exc)), d = exc$D_raw, grp = "within")
)
da$y <- row_mid[da$row] + ifelse(da$row == 1, stack_offset(da$row, da$d, step = 0.13, cap = 0.42),
                                 ifelse(da$row == 2, stack_offset(da$row, da$d, step = 0.17, cap = 1.15), 0))
exc_lab <- data.frame(
  y = row_mid[2 + seq_len(nrow(exc))],
  d = exc$D_raw,
  label = sprintf("%s \u2192 %s, %d SNPs\n%d days apart, %s, %s",
                  exc$isolate_1, exc$isolate_2, exc$D_raw, exc$delta_days,
                  ifelse(exc$admission_axis == "same", "same admission", "two admissions"),
                  ifelse(exc$specimen_1 == exc$specimen_2, paste("both", spec_short(exc$specimen_1)),
                         paste(spec_short(exc$specimen_1), "\u2192", spec_short(exc$specimen_2))))
)

pB <- ggplot(da, aes(d, y)) +
  ref_lines +
  geom_point(aes(colour = grp), size = 1.25, alpha = 0.9, stroke = 0) +
  geom_label(data = exc_lab, aes(label = label), hjust = 1, nudge_x = -0.06, size = 2.1, family = FONT, colour = COL$ink,
             lineheight = 0.95, fill = "white", border.colour = NA, label.padding = unit(0.5, "mm")) +
  annotate("text", x = 1350, y = row_mid[1], hjust = 1, size = 2.1, family = FONT, colour = COL$ink2, lineheight = 0.95,
           label = "Dotted line: H99 95th percentile (2 SNPs)\nDashed line: reference value (5 SNPs)") +
  scale_x_snp() +
  scale_y_continuous(breaks = row_mid, labels = lab_rows, limits = c(0, sum(row_h)), expand = c(0, 0)) +
  scale_colour_manual(values = c(h99 = COL$muted, within = COL$within), guide = "none") +
  labs(x = NULL, y = NULL, tag = "B") +
  theme_fig() +
  theme(axis.text.y = element_text(colour = COL$ink, lineheight = 0.9),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "#f4f3f0"))

# Figure 2C: intra-patient vs inter-patient distances (representative isolates)
lab_w <- sprintf("Intra-patient pairs (%d)", nrow(rb))
lab_c <- sprintf("Inter-patient pairs, representative isolates (%s)", format(nrow(rc), big.mark = ","))
dc <- rbind(data.frame(group = lab_w, d = rb$D_raw), data.frame(group = lab_c, d = rc$D_raw))
dc$group <- factor(dc$group, levels = c(lab_w, lab_c))
frac_le <- n_le / nrow(rc)
# Tags for the three exception steps of the within-patient curve, matching the rows of panel B.
exc_tag <- data.frame(d = exc$D_raw, y = vapply(exc$D_raw, function(v) mean(rb$D_raw < v), numeric(1)), label = exc$patient_id)

pC <- ggplot(dc, aes(d, colour = group)) +
  annotate("rect", xmin = 0, xmax = max_b, ymin = -Inf, ymax = Inf, fill = COL$band, alpha = 0.7) +
  ref_lines +
  stat_ecdf(geom = "step", linewidth = 0.6, pad = FALSE) +
  geom_text(data = exc_tag, aes(d, y - 0.035, label = label), inherit.aes = FALSE, hjust = 0.5, size = 1.9,
            family = FONT, colour = COL$ink) +
  annotate("segment", x = max_b, xend = max_b, y = 0, yend = frac_le, linewidth = 0.3, colour = COL$between) +
  annotate("point", x = max_b, y = frac_le, size = 1.4, colour = COL$between) +
  annotate("segment", x = 175, xend = max_b * 0.97, y = 0.45, yend = frac_le + 0.02, linewidth = 0.25, colour = COL$ink2) +
  annotate("text", x = 170, y = 0.45, hjust = 1, vjust = 0.5, size = 2.1, family = FONT, colour = COL$ink, lineheight = 0.95,
           label = sprintf("%.1f%% of inter-patient pairs (%s/%s)\nare \u2264 %d SNPs, the intra-patient maximum",
                           100 * frac_le, format(n_le, big.mark = ","), format(nrow(rc), big.mark = ","), max_b)) +
  annotate("segment", x = 2.6, xend = 0.08, y = 0.1, yend = 0.012, linewidth = 0.25, colour = COL$ink2) +
  annotate("label", x = 2.8, y = 0.1, hjust = 0, vjust = 0.5, size = 2.1, family = FONT, colour = COL$ink2, lineheight = 0.95,
           fill = COL$band, border.colour = NA, label.padding = unit(0.6, "mm"),
           label = sprintf("%d/%s inter-patient pairs at 0 SNPs:\nCn13\u2013Cn14 (Cn13 mixed DNA) and Cn101\u2013Cn103",
                           nrow(zero_c), format(nrow(rc), big.mark = ","))) +
  annotate("text", x = 14, y = 1.035, size = 2.1, family = FONT, colour = COL$ink2,
           label = sprintf("overlap of ranges, 0\u2013%d", max_b)) +
  scale_x_snp() +
  scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_colour_manual(values = setNames(c(COL$within, COL$between), c(lab_w, lab_c)), name = NULL) +
  labs(x = "SNP differences between two isolates (pseudo-log scale)", y = "Cumulative proportion of pairs", tag = "C") +
  theme_fig() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, -2, 0, "mm"))

fig <- pA / pB / pC + plot_layout(heights = c(1.45, 1, 1.15))
save_fig(fig, "Figure2", width_cm = 17.5, height_cm = 22)

# Supplementary figure: the 40 intra-patient pairs
interval_lab <- c("0" = "Same day", "1-7" = "1\u20137 days", "8-30" = "8\u201330 days", ">30" = ">30 days")
kind_lab <- c(
  "same.same" = "Same admission, same specimen type",
  "same.different" = "Same admission, different specimen type",
  "different.same" = "Different admission, same specimen type",
  "different.different" = "Different admission, different specimen type"
)
kind_shape <- setNames(c(16, 17, 1, 2), kind_lab)

pat_lab <- setNames(sprintf("%s (%d)", pat_order$patient_id, pat_order$n_strains), pat_order$patient_id)
ds <- data.frame(
  patient = factor(pat_lab[rb$patient_id], levels = pat_lab),
  d = rb$D_raw,
  interval = factor(interval_lab[rb$interval_class], levels = interval_lab),
  kind = factor(kind_lab[paste(rb$admission_axis, rb$site_axis, sep = ".")], levels = kind_lab),
  pair = paste(rb$isolate_1, rb$isolate_2, sep = "\u2013")
)
# Rows with tied distances (P097, P100) get extra height so tied points stay apart.
ties <- tapply(paste(ds$patient, ds$d), ds$patient, function(k) max(table(k)))
srow_h <- pmax(1, 0.45 * ties[levels(ds$patient)] + 0.3)
srow_mid <- setNames(cumsum(srow_h) - srow_h / 2, levels(ds$patient))
ds$y <- srow_mid[as.character(ds$patient)] + stack_offset(ds$patient, ds$d, step = 0.45, cap = 0.9)
ds_exc <- ds[ds$d > 5, ]

pS <- ggplot(ds, aes(d, y)) +
  ref_lines +
  geom_point(aes(colour = interval, shape = kind), size = 1.5, stroke = 0.55) +
  geom_text(data = ds_exc, aes(label = pair), hjust = -0.18, size = 2.1, family = FONT, colour = COL$ink2) +
  scale_x_snp() +
  scale_y_continuous(breaks = unname(srow_mid), labels = pat_lab, limits = c(-0.1, sum(srow_h) + 0.1), expand = c(0, 0)) +
  scale_colour_manual(values = setNames(COL$interval, interval_lab), drop = FALSE, name = "Interval between samplings") +
  scale_shape_manual(values = kind_shape, drop = FALSE, name = "Admission and specimen type") +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 1.8)),
         shape = guide_legend(order = 2, override.aes = list(size = 1.6))) +
  labs(x = "SNP differences between two isolates (pseudo-log scale)", y = "Patient (number of isolates)") +
  theme_fig() +
  theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "#f4f3f0"),
        axis.text.y = element_text(size = 6),
        legend.position = "inside", legend.position.inside = c(0.995, 0.02), legend.justification = c(1, 0),
        legend.box = "horizontal", legend.box.just = "bottom",
        legend.background = element_rect(fill = "white", colour = NA))

save_fig(pS, "FigureS_within_pairs", width_cm = 17.5, height_cm = 12)
message("Figure 2 and the supplementary pair plot done")
