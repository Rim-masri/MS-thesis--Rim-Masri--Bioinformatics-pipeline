nano make_pheno_width_majdool.R
base <- "/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
book <- file.path(base, "Book1.csv")
groups <- file.path(base, "sample_groups_majdool.tsv")
out_plink2 <- file.path(base, "pheno_width_majdool.plink2.txt")

df <- read.csv(book, check.names = FALSE, stringsAsFactors = FALSE)
names(df) <- trimws(names(df))

grp <- read.table(groups, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(grp) <- trimws(names(grp))

serial   <- trimws(as.character(df[[1]]))
cultivar <- trimws(as.character(df[[4]]))
quality  <- trimws(as.character(df[[5]]))
width    <- suppressWarnings(as.numeric(df$width))

grp$sample_id <- trimws(as.character(grp$sample_id))

df$id <- paste(cultivar, quality, serial, sep = "_")
df$width_num <- width

mj <- df[tolower(cultivar) == "majdool" & !is.na(df$width_num), c("id", "width_num")]
names(mj) <- c("IID", "width")

mj <- merge(grp[, "sample_id", drop = FALSE], mj,
            by.x = "sample_id", by.y = "IID", all.x = TRUE)

mj <- mj[!is.na(mj$width), ]

out <- mj[, c("sample_id", "width")]
names(out) <- c("IID", "width")

write.table(out, out_plink2, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Created:", out_plink2, "\n")
cat("Rows:", nrow(out), "\n")
print(head(out, 10))

nano make_pheno_width_majdool.sh
#!/bin/bash
set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
SRC="${BASE}/pheno_width_majdool.plink2.txt"
OUT="${BASE}/pheno_width_majdool.iid.txt"

if [[ ! -f "${SRC}" ]]; then
    echo "ERROR: source phenotype file not found: ${SRC}"
    exit 1
fi

awk 'BEGIN{OFS="\t"} NR==1{print "IID","width"; next} NR>1{print $1,$2}' "${SRC}" > "${OUT}"

echo "Created phenotype file:"
ls -lh "${OUT}"
echo
echo "Preview:"
head "${OUT}"

Rscript make_pheno_width_majdool.R
bash make_pheno_width_majdool.sh

nano merge_and_plot_width_majdool_genomewide.R
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outdir <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

libdir <- Sys.getenv("R_LIBS_USER")
if (libdir == "") {
  libdir <- "/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/Rlib"
}
dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(libdir, .libPaths()))

if (!requireNamespace("qqman", quietly = TRUE)) {
  stop("Package qqman is not available in R_LIBS_USER. Install it first into Rlib.")
}
library(qqman)

gwas <- read.table(
  infile,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

cat("Columns read:\n")
print(colnames(gwas))
cat("Rows before filtering:", nrow(gwas), "\n")

if (!("#CHROM" %in% colnames(gwas))) stop("Column '#CHROM' not found.")
if (!("POS" %in% colnames(gwas))) stop("Column 'POS' not found.")
if (!("ID" %in% colnames(gwas))) stop("Column 'ID' not found.")
if (!("P" %in% colnames(gwas))) stop("Column 'P' not found.")

gwas$CHR_LABEL <- trimws(as.character(gwas$`#CHROM`))
gwas$P <- suppressWarnings(as.numeric(gwas$P))
gwas$BP <- suppressWarnings(as.integer(gwas$POS))
gwas$SNP <- ifelse(
  is.na(gwas$ID) | gwas$ID == ".",
  paste0(gwas$CHR_LABEL, ":", gwas$BP),
  gwas$ID
)

chr_levels <- c(
  as.character(1:14),
  "Chr14_male.SDR-oriented",
  as.character(15:18)
)

gwas$CHR_FACTOR <- factor(gwas$CHR_LABEL, levels = chr_levels, ordered = TRUE)

gwas <- gwas[
  is.finite(gwas$P) &
    gwas$P > 0 &
    gwas$P <= 1 &
    !is.na(gwas$BP) &
    !is.na(gwas$CHR_FACTOR),
]

cat("Rows after filtering:", nrow(gwas), "\n")
if (nrow(gwas) == 0) stop("No valid rows after filtering.")

gwas$CHR <- as.integer(gwas$CHR_FACTOR)

n_tests <- nrow(gwas)
bonf_p <- 0.05 / n_tests
bonf_line <- -log10(bonf_p)
max_y <- max(-log10(gwas$P), bonf_line, na.rm = TRUE)

write.table(
  gwas[, c("CHR_LABEL", "CHR", "BP", "SNP", "P")],
  file = file.path(outdir, "majdool_width_manhattan_input.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

png(
  filename = file.path(outdir, "majdool_width_genomewide_manhattan.png"),
  width = 3200,
  height = 1600,
  res = 240
)

par(mar = c(5, 5, 4, 2) + 0.1)

manhattan(
  gwas,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  col = c("#4C78A8", "#F58518"),
  suggestiveline = FALSE,
  genomewideline = FALSE,
  ylim = c(0, max_y + 0.5),
  chrlabs = chr_levels,
  main = paste0(
    "Majdool seed width GWAS (genome-wide)\nBonferroni threshold = ",
    signif(bonf_p, 3)
  )
)

abline(h = bonf_line, col = "red", lwd = 2, lty = 2)

dev.off()

writeLines(
  c(
    paste("Total tests:", n_tests),
    paste("Bonferroni p-value:", bonf_p),
    paste("Bonferroni -log10 threshold:", bonf_line),
    paste("Max plotted y:", max_y),
    paste("Chromosomes in plot:", paste(chr_levels, collapse = ", "))
  ),
  con = file.path(outdir, "majdool_width_plot_summary.txt")
)

nano run_all_gwas_width_majdool.slurm
#!/bin/bash
#SBATCH --job-name=majdool_width_genomewide
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/majdool_width_genomewide_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/majdool_width_genomewide_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
STITCH_BASE="/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/stitch"
PHENO="${BASE}/pheno_width_majdool.iid.txt"
OUTBASE="${BASE}/gwas_majdool_width_genomewide"
MERGED="${OUTBASE}/majdool_width_all_chunks.merged.glm.linear"
PLOT_R="${BASE}/merge_and_plot_width_majdool_genomewide.R"
R_LIB_DIR="${BASE}/Rlib"

mkdir -p "${BASE}/logs"
mkdir -p "${OUTBASE}"
mkdir -p "${R_LIB_DIR}"

export R_LIBS_USER="${R_LIB_DIR}"

module load plink/2.0
module load R/4.3.2/gnu-12.2.0

if [[ ! -f "${PHENO}" ]]; then
    echo "ERROR: phenotype file not found: ${PHENO}"
    exit 1
fi

if [[ ! -f "${PLOT_R}" ]]; then
    echo "ERROR: plotting script not found: ${PLOT_R}"
    exit 1
fi

find "${STITCH_BASE}" -type f -name "*.vcf.gz" | sort > "${OUTBASE}/vcf_list.txt"
NVCF=$(wc -l < "${OUTBASE}/vcf_list.txt")

if [[ "${NVCF}" -eq 0 ]]; then
    echo "ERROR: no VCF files found under ${STITCH_BASE}"
    exit 1
fi

while read -r VCF; do
    REL=$(echo "${VCF}" | sed "s#${STITCH_BASE}/##" | sed 's#/#_#g' | sed 's/\.vcf\.gz$//')
    PREFIX="${OUTBASE}/${REL}"

    plink2 \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --vcf "${VCF}" \
      --pheno iid-only "${PHENO}" \
      --pheno-name width \
      --glm hide-covar allow-no-covars \
      --out "${PREFIX}"

done < "${OUTBASE}/vcf_list.txt"

find "${OUTBASE}" -type f -name "*.width.glm.linear" | sort > "${OUTBASE}/glm_list.txt"
FIRST=$(head -n 1 "${OUTBASE}/glm_list.txt")
head -n 1 "${FIRST}" > "${MERGED}"

while read -r F; do
    tail -n +2 "${F}" >> "${MERGED}"
done < "${OUTBASE}/glm_list.txt"

Rscript "${PLOT_R}" "${MERGED}" "${OUTBASE}"

nano plot_only_width_majdool.slurm
#!/bin/bash
#SBATCH --job-name=plot_width_majdool_only
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/plot_width_majdool_only_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/plot_width_majdool_only_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
MERGED="${BASE}/gwas_majdool_width_genomewide/majdool_width_all_chunks.merged.glm.linear"
OUTBASE="${BASE}/gwas_majdool_width_genomewide"
PLOT_R="${BASE}/merge_and_plot_width_majdool_genomewide.R"
R_LIB_DIR="${BASE}/Rlib"

mkdir -p "${BASE}/logs"
mkdir -p "${R_LIB_DIR}"

export R_LIBS_USER="${R_LIB_DIR}"

module load R/4.3.2/gnu-12.2.0

if [[ ! -f "${MERGED}" ]]; then
    echo "ERROR: merged GWAS file not found: ${MERGED}"
    exit 1
fi

Rscript "${PLOT_R}" "${MERGED}" "${OUTBASE}"

sbatch run_all_gwas_width_majdool.slurm
