cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia

nano make_pheno_length_saqi.R
base <- "/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
book <- file.path(base, "Book1.csv")
groups <- file.path(base, "sample_groups_saqi.tsv")
out_plink2 <- file.path(base, "pheno_length_saqi.plink2.txt")

df <- read.csv(book, check.names = FALSE, stringsAsFactors = FALSE)
names(df) <- trimws(names(df))

grp <- read.table(groups, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(grp) <- trimws(names(grp))

df$serial   <- trimws(as.character(df[[1]]))
df$cultivar <- trimws(as.character(df[[4]]))
df$quality  <- trimws(as.character(df[[5]]))
df$length   <- suppressWarnings(as.numeric(df$length))
grp$sample_id <- trimws(as.character(grp$sample_id))

df$id <- paste(df$cultivar, df$quality, df$serial, sep = "_")

sq <- df[tolower(df$cultivar) == "saqi" & !is.na(df$length), c("id", "length")]
names(sq) <- c("IID", "length")

sq <- merge(grp[, "sample_id", drop = FALSE], sq, by.x = "sample_id", by.y = "IID", all.x = TRUE)
sq <- sq[!is.na(sq$length), ]

out <- sq[, c("sample_id", "length")]
names(out) <- c("IID", "length")

write.table(out, out_plink2, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Created:", out_plink2, "\n")
cat("Rows:", nrow(out), "\n")
print(head(out, 10))

module load R/4.3.2/gnu-12.2.0
Rscript make_pheno_length_saqi.R

nano make_pheno_length_saqi.sh

#!/bin/bash
set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
SRC="${BASE}/pheno_length_saqi.plink2.txt"
OUT="${BASE}/pheno_length_saqi.iid.txt"

if [[ ! -f "${SRC}" ]]; then
    echo "ERROR: source phenotype file not found: ${SRC}"
    exit 1
fi

awk 'BEGIN{OFS="\t"} NR==1{print "IID","length"; next} NR>1{print $1,$2}' "${SRC}" > "${OUT}"

echo "Created phenotype file:"
ls -lh "${OUT}"
echo
echo "Preview:"
head "${OUT}"

bash make_pheno_length_saqi.sh

nano merge_and_plot_length_saqi_genomewide.R

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

chr_counts <- table(gwas$CHR_FACTOR)
cat("Chromosome counts kept:\n")
print(chr_counts)

n_tests <- nrow(gwas)
bonf_p <- 0.05 / n_tests
bonf_line <- -log10(bonf_p)
max_y <- max(-log10(gwas$P), bonf_line, na.rm = TRUE)

write.table(
  gwas[, c("CHR_LABEL", "CHR", "BP", "SNP", "P")],
  file = file.path(outdir, "saqi_length_manhattan_input.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

png(
  filename = file.path(outdir, "saqi_length_genomewide_manhattan.png"),
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
    "Saqi seed length GWAS (genome-wide)\nBonferroni threshold = ",
    signif(bonf_p, 3)
  )
)

abline(h = bonf_line, col = "red", lwd = 2, lty = 2)
mtext(
  paste0("Bonferroni line: -log10(", signif(bonf_p, 3), ") = ", round(bonf_line, 2)),
  side = 3,
  line = 0.2,
  col = "red",
  cex = 0.9
)

dev.off()

writeLines(
  c(
    paste("Total tests:", n_tests),
    paste("Bonferroni p-value:", bonf_p),
    paste("Bonferroni -log10 threshold:", bonf_line),
    paste("Max plotted y:", max_y),
    paste("Chromosomes in plot:", paste(chr_levels, collapse = ", "))
  ),
  con = file.path(outdir, "saqi_length_plot_summary.txt")
)

cat("Done.\n")
cat("Bonferroni p-value:", bonf_p, "\n")
cat("Bonferroni line:", bonf_line, "\n")
cat("Plot written to:", file.path(outdir, "saqi_length_genomewide_manhattan.png"), "\n")

nano run_all_gwas_length_saqi.slurm
#!/bin/bash
#SBATCH --job-name=saqi_length_genomewide
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/saqi_length_genomewide_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/saqi_length_genomewide_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
STITCH_BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/stitch"
PHENO="${BASE}/pheno_length_saqi.iid.txt"
OUTBASE="${BASE}/gwas_saqi_length_genomewide"
MERGED="${OUTBASE}/saqi_length_all_chunks.merged.glm.linear"
PLOT_R="${BASE}/merge_and_plot_length_saqi_genomewide.R"
R_LIB_DIR="${BASE}/Rlib"

mkdir -p "${BASE}/logs"
mkdir -p "${OUTBASE}"
mkdir -p "${R_LIB_DIR}"

export R_LIBS_USER="${R_LIB_DIR}"

module load plink/2.0
module load R/4.3.2/gnu-12.2.0

echo "BASE=${BASE}"
echo "STITCH_BASE=${STITCH_BASE}"
echo "PHENO=${PHENO}"
echo "OUTBASE=${OUTBASE}"
echo "MERGED=${MERGED}"
echo "PLOT_R=${PLOT_R}"
echo "R_LIBS_USER=${R_LIBS_USER}"

if [[ ! -f "${PHENO}" ]]; then
    echo "ERROR: phenotype file not found: ${PHENO}"
    exit 1
fi

if [[ ! -f "${PLOT_R}" ]]; then
    echo "ERROR: plotting script not found: ${PLOT_R}"
    exit 1
fi

echo "Phenotype preview:"
head "${PHENO}"

echo "Finding all VCF chunks..."
find "${STITCH_BASE}" -type f -name "*.vcf.gz" | sort > "${OUTBASE}/vcf_list.txt"

NVCF=$(wc -l < "${OUTBASE}/vcf_list.txt")
echo "Number of VCF chunks found: ${NVCF}"

if [[ "${NVCF}" -eq 0 ]]; then
    echo "ERROR: no VCF files found under ${STITCH_BASE}"
    exit 1
fi

while read -r VCF; do
    REL=$(echo "${VCF}" | sed "s#${STITCH_BASE}/##" | sed 's#/#_#g' | sed 's/\.vcf\.gz$//')
    PREFIX="${OUTBASE}/${REL}"

    echo "----------------------------------------"
    echo "Running GWAS on: ${VCF}"
    echo "Output prefix: ${PREFIX}"

    plink2 \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --vcf "${VCF}" \
      --pheno iid-only "${PHENO}" \
      --pheno-name length \
      --glm hide-covar allow-no-covars \
      --out "${PREFIX}"

done < "${OUTBASE}/vcf_list.txt"

echo "Collecting GWAS result files..."
find "${OUTBASE}" -type f -name "*.length.glm.linear" | sort > "${OUTBASE}/glm_list.txt"

NGLM=$(wc -l < "${OUTBASE}/glm_list.txt")
echo "Number of .length.glm.linear files found: ${NGLM}"

if [[ "${NGLM}" -eq 0 ]]; then
    echo "ERROR: no .length.glm.linear files were produced."
    exit 1
fi

FIRST=$(head -n 1 "${OUTBASE}/glm_list.txt")
head -n 1 "${FIRST}" > "${MERGED}"

while read -r F; do
    tail -n +2 "${F}" >> "${MERGED}"
done < "${OUTBASE}/glm_list.txt"

echo "Merged GWAS file created:"
ls -lh "${MERGED}"

echo "Checking chromosome counts in merged file..."
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="#CHROM") c=i} NR>1{n[$c]++} END{for(k in n) print k, n[k]}' "${MERGED}" | sort -V

echo "Generating Manhattan plot..."
Rscript "${PLOT_R}" "${MERGED}" "${OUTBASE}"

echo "Final outputs:"
ls -lh "${OUTBASE}"

echo "Finished at $(date)"

nano plot_only_length_saqi.slurm
#!/bin/bash
#SBATCH --job-name=plot_length_saqi_only
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/plot_length_saqi_only_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/plot_length_saqi_only_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
MERGED="${BASE}/gwas_saqi_length_genomewide/saqi_length_all_chunks.merged.glm.linear"
OUTBASE="${BASE}/gwas_saqi_length_genomewide"
PLOT_R="${BASE}/merge_and_plot_length_saqi_genomewide.R"
R_LIB_DIR="${BASE}/Rlib"

mkdir -p "${BASE}/logs"
mkdir -p "${R_LIB_DIR}"

export R_LIBS_USER="${R_LIB_DIR}"

module load R/4.3.2/gnu-12.2.0

if [[ ! -f "${MERGED}" ]]; then
    echo "ERROR: merged GWAS file not found: ${MERGED}"
    exit 1
fi

if [[ ! -f "${PLOT_R}" ]]; then
    echo "ERROR: plotting script not found: ${PLOT_R}"
    exit 1
fi

Rscript "${PLOT_R}" "${MERGED}" "${OUTBASE}"

echo "Plot-only job finished at $(date)"

