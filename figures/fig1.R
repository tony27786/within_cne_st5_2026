# Figure 1: flow of patients and isolates, with the sequenced samples and the analysis sets.
# Inputs: the cleaned clinical tables (isolates, admissions), data/samples.tsv, the G1 rulers (ruler_B_pairs, pairwise_L0) and the A1 design table; sampling dates with the errata applied (fig_common.R).
# Run: Rscript figures/fig1.R

local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
  source(file.path(dir, "fig_common.R"))
})

p_iso <- in_path("data", "isolates.csv")
p_adm <- in_path("data", "admissions.csv")
p_smp <- in_path("data", "samples.tsv")
p_rb <- in_path(G1_OUT, "02_l0", "main", "ruler_B_pairs.tsv")
p_l0 <- in_path(G1_OUT, "02_l0", "main", "pairwise_L0.tsv")
p_des <- in_path("association", "design.tsv")

read_clean <- function(p, cols) {
  x <- read.csv(p, colClasses = "character", encoding = "UTF-8", check.names = FALSE)
  names(x)[1] <- sub("^\ufeff", "", names(x)[1])
  x[, cols]
}
iso <- read_clean(p_iso, c("sample", "patient_id", "sampling_date", "st", "wgs_role"))
adm <- read_clean(p_adm, c("admission_id", "patient_id"))
smp <- read_tsv(p_smp, colClasses = "character")
rb <- read_tsv(p_rb)
l0 <- read_tsv(p_l0)
des <- read_tsv(p_des, colClasses = "character")

# The released tables hold only the patients with an ST5 isolate (and P001). Totals of the patients without a confirmed ST5 isolate, whose records are not released:
# 27 patients (28 admissions, 31 isolates), 21 with other STs only and 6 whose isolates have no assigned ST.
n_excl_pat <- 27; n_excl_adm <- 28; n_excl_iso <- 31; n_excl_other <- 21; n_excl_untyped <- 6
n_arch_iso <- nrow(iso) + n_excl_iso + 1        # + isolates of the excluded patients + Cn83
n_pat_clean <- length(unique(adm$patient_id)) + n_excl_pat
n_arch_pat <- n_pat_clean + 1                   # + the patient without a medical record
idates <- isolate_dates()                       # fig_common.R: sampling dates with the clinical errata applied
months <- range(substr(idates[iso$sample], 1, 7))

has_st5 <- tapply(iso$st == "5", iso$patient_id, any)
st5_pat <- names(has_st5)[has_st5]
reps <- iso[iso$wgs_role == "representative", ]
extras <- iso[iso$wgs_role == "within_patient_extra", ]
grp <- table(ifelse(des$cm == "1", ifelse(des$hiv_positive == "1", "pos", "neg"), "noncm"))
n_multi <- length(unique(extras$patient_id))
cross <- sort(unique(rb$patient_id[rb$admission_axis == "different"]))
n_genome <- nrow(reps) + nrow(extras)
# ST5 culture tubes by laboratory MLST: the ST5 isolates of the cleaned table + Cn83 (removed with the patient without a record) + Cn169b (second tube of Cn169).
n_st5_tubes <- sum(iso$st == "5") + 2
n_c <- sum(l0$pair_type == "rep_rep")
n_a1 <- sum(des$a1_main == "1")

# Boxes: x, y = centre in a 100 x 100 canvas; w, h = size
fmt <- function(n) format(n, big.mark = ",")
boxes <- data.frame(
  id = c("arch", "ex1", "clean", "ex2", "st5", "ex3", "main", "multi", "seq", "g1", "g2", "a1"),
  x = c(32, 81, 32, 81, 32, 81, 32, 49.5, 82, 16, 49, 82),
  y = c(95, 87, 79, 71, 63, 55, 45, 31.5, 37, 11, 11, 11),
  w = c(46, 34, 46, 34, 46, 34, 52, 27, 32, 30, 30, 32),
  h = c(5.5, 5.5, 5, 8.5, 4, 6, 8, 8.5, 19, 10, 10, 10),
  label = c(
    "",   # drawn below with the species name in italics
    "Excluded: 1 patient without a medical record\n(isolate Cn83, ST5)",
    sprintf("%d patients, %d admissions\nwith cleaned clinical data", n_pat_clean, nrow(adm) + n_excl_adm),
    sprintf("Excluded: %d patients without a\nconfirmed ST5 isolate (%d isolates):\n%d with other STs only, %d whose\nisolates had no assigned ST",
            n_excl_pat, n_excl_iso, n_excl_other, n_excl_untyped),
    sprintf("%d patients with at least one ST5 isolate", length(st5_pat)),
    "Described separately: 1 patient whose first\nisolate was ST2 (Cn76, ST2; later Cn84, ST5)",
    sprintf("%d patients in the ST5 analysis\nCM %d (HIV-positive %d, HIV-negative %d); non-CM %d\n%d representative isolates (earliest archived isolate per patient)",
            nrow(reps), sum(grp[c("pos", "neg")]), grp[["pos"]], grp[["neg"]], grp[["noncm"]], nrow(reps)),
    sprintf("Of these, %d patients with \u2265 2 ST5\nisolates: %d additional isolates,\n%d intra-patient pairs (%d patients\nwith isolates from different admissions)",
            n_multi, nrow(extras), nrow(rb), length(cross)),
    sprintf("Whole-genome sequencing: %d samples\n\n%d ST5 isolates (%d representative + %d additional)\nCn169b, second culture tube of Cn169\nCn76 and Cn84 (case report)\n%d H99 sequencing samples (reference control)\n\nST5 isolates by MLST: %d tubes \u2212 Cn169b\n\u2212 Cn83 \u2212 Cn84 = %d",
            nrow(smp), n_genome, nrow(reps), nrow(extras), sum(smp$role == "technical_replicate"), n_st5_tubes, n_genome),
    sprintf("SNP distances\n%d ST5 isolates\n%d intra-patient pairs\n%s inter-patient pairs\n(representative isolates)\nH99, 10 pairs", n_genome, nrow(rb), fmt(n_c)),
    sprintf("Population structure\n%d representative isolates\nML tree, PCA, genomic\nrelationship matrix", nrow(reps)),
    sprintf("Genome\u2013phenotype analysis\n%d patients (Cn13 excluded:\nmixed DNA); complete cases\nper outcome (Table 3)", n_a1)
  ),
  stringsAsFactors = FALSE
)
b <- function(id) boxes[boxes$id == id, ]
bottom <- function(id) b(id)$y - b(id)$h / 2
top <- function(id) b(id)$y + b(id)$h / 2
left <- function(id) b(id)$x - b(id)$w / 2

arr <- arrow(length = unit(1.4, "mm"), type = "closed")
main_x <- b("arch")$x
v <- data.frame(x = main_x, xend = main_x,
                y = c(bottom("arch"), bottom("clean"), bottom("st5")),
                yend = c(top("clean"), top("st5"), top("main")))
# Exclusions go to the right; the patients with several isolates are a side branch of the 160, not a step towards the analysis sets.
h <- data.frame(x = main_x, xend = c(left("ex1"), left("ex2"), left("ex3"), left("multi")),
                y = c(b("ex1")$y, b("ex2")$y, b("ex3")$y, b("multi")$y))
h$yend <- h$y
# from the analysis population down to the three analysis sets
split_y <- 22
fan <- data.frame(x = c(b("g1")$x, b("g2")$x, b("a1")$x), xend = c(b("g1")$x, b("g2")$x, b("a1")$x),
                  y = split_y, yend = top("g1"))

pS1 <- ggplot() +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
            fill = ifelse(boxes$id %in% c("ex1", "ex2", "ex3"), "#f4f3f0", "white"), colour = COL$ink2, linewidth = 0.3) +
  geom_text(data = boxes, aes(x, y, label = label), size = 2.1, family = FONT, colour = COL$ink, lineheight = 0.95) +
  annotate("text", x = b("arch")$x, y = b("arch")$y + 1, parse = TRUE, size = 2.1, family = FONT, colour = COL$ink,
           label = sprintf("\"%d patients with archived \"*italic(\"C. neoformans\")*\" isolates\"", n_arch_pat)) +
  annotate("text", x = b("arch")$x, y = b("arch")$y - 1.1, size = 2.1, family = FONT, colour = COL$ink,
           label = sprintf("%d isolates, sampled %s to %s", n_arch_iso, months[1], months[2])) +
  geom_segment(data = v, aes(x = x, xend = xend, y = y, yend = yend), arrow = arr, linewidth = 0.3, colour = COL$ink2) +
  geom_segment(data = h, aes(x = x, xend = xend, y = y, yend = yend), arrow = arr, linewidth = 0.3, colour = COL$ink2) +
  annotate("segment", x = main_x, xend = main_x, y = bottom("main"), yend = split_y, linewidth = 0.3, colour = COL$ink2) +
  annotate("segment", x = b("g1")$x, xend = b("a1")$x, y = split_y, yend = split_y, linewidth = 0.3, colour = COL$ink2) +
  geom_segment(data = fan, aes(x = x, xend = xend, y = y, yend = yend), arrow = arr, linewidth = 0.3, colour = COL$ink2) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5.5, 98.5), expand = c(0, 0)) +
  theme_void(base_family = FONT)

save_fig(pS1, "Figure1", width_cm = 17.5, height_cm = 15)
message("Figure 1 done")
