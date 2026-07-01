nano prep_saqi_chr01_permutation.sh
#!/usr/bin/env bash
set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
WORK="/ibex/project/c2042/user_masrir/datepalm_metaxenia"

mkdir -p "${BASE}/merged" "${BASE}/pruned_chr01" "${BASE}/permutation_chr01" "${BASE}/logs" "${BASE}/plots"

echo "BASE=${BASE}"
echo "WORK=${WORK}"

echo "[check 1] Confirm stitched Chr01 VCFs exist"
find "${BASE}/stitch/Chr01" -type f -name "*.vcf.gz" | head
test -n "$(find "${BASE}/stitch/Chr01" -type f -name "*.vcf.gz" | head -n 1)"

echo "[check 2] Confirm sample groups file exists"
test -s "${BASE}/sample_groups_saqi.tsv"
head -n 3 "${BASE}/sample_groups_saqi.tsv"

echo "[1/5] Concatenate chr01 VCF chunks"
bcftools concat -Oz -o "${BASE}/merged/Chr01.concat.vcf.gz" ${BASE}/stitch/Chr01/*/*.vcf.gz
bcftools index -f "${BASE}/merged/Chr01.concat.vcf.gz"

echo "[2/5] Filter INFO_SCORE > 0.3 and SNPs only"
bcftools view -v snps -i 'INFO/INFO_SCORE>0.3' "${BASE}/merged/Chr01.concat.vcf.gz" \
  -Oz -o "${BASE}/merged/Chr01.info30.snps.vcf.gz"
bcftools index -f "${BASE}/merged/Chr01.info30.snps.vcf.gz"

echo "[3/5] Convert VCF to PLINK2"
module load plink/2.0
plink2 \
  --vcf "${BASE}/merged/Chr01.info30.snps.vcf.gz" \
  --allow-extra-chr \
  --set-all-var-ids '@:#' \
  --max-alleles 2 \
  --make-pgen \
  --out "${BASE}/pruned_chr01/chr01_info30"

echo "[4/5] LD pruning"
plink2 \
  --pfile "${BASE}/pruned_chr01/chr01_info30" \
  --allow-extra-chr \
  --indep-pairwise 50 5 0.2 \
  --out "${BASE}/pruned_chr01/chr01_info30_pruned"

echo "[5/5] Export pruned hard genotype matrix"
plink2 \
  --pfile "${BASE}/pruned_chr01/chr01_info30" \
  --allow-extra-chr \
  --extract "${BASE}/pruned_chr01/chr01_info30_pruned.prune.in" \
  --export A \
  --out "${BASE}/pruned_chr01/chr01_info30_pruned_A"

echo "Done."
echo "Check these files:"
ls -lh "${BASE}/pruned_chr01"/chr01_info30_pruned.prune.in
ls -lh "${BASE}/pruned_chr01"/chr01_info30_pruned_A.raw

nano permute_saqi_chr01_fisher.R
args <- commandArgs(trailingOnly = TRUE)

raw_file   <- args[1]
group_file <- args[2]
out_dir    <- args[3]
nperm      <- as.integer(args[4])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grp <- read.table(group_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("sample_id", "group") %in% names(grp)))

x <- read.table(raw_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

meta_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
stopifnot(all(meta_cols %in% names(x)))

sample_ids_raw <- x$IID
geno <- x[, !(names(x) %in% meta_cols), drop = FALSE]

clean_varname <- function(v) sub("_.*$", "", v)
var_ids <- clean_varname(names(geno))

map <- data.frame(
  colname = names(geno),
  varid = var_ids,
  stringsAsFactors = FALSE
)

parts <- strsplit(map$varid, ":", fixed = TRUE)
map$CHROM <- vapply(parts, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1))
map$POS   <- vapply(parts, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))

keep_samples <- sample_ids_raw %in% grp$sample_id
x <- x[keep_samples, , drop = FALSE]
sample_ids <- x$IID
geno <- x[, !(names(x) %in% meta_cols), drop = FALSE]

grp <- grp[match(sample_ids, grp$sample_id), , drop = FALSE]
stopifnot(all(sample_ids == grp$sample_id))

geno_mat <- as.matrix(geno)
storage.mode(geno_mat) <- "numeric"

run_one_scan <- function(group_vec, geno_mat) {
  hi <- which(group_vec == "high")
  lo <- which(group_vec == "low")

  pvals <- rep(NA_real_, ncol(geno_mat))
  odds  <- rep(NA_real_, ncol(geno_mat))

  for (j in seq_len(ncol(geno_mat))) {
    g <- geno_mat[, j]

    hi_nonmiss <- !is.na(g[hi])
    lo_nonmiss <- !is.na(g[lo])

    high_alt <- sum(g[hi][hi_nonmiss], na.rm = TRUE)
    low_alt  <- sum(g[lo][lo_nonmiss], na.rm = TRUE)

    high_called <- sum(hi_nonmiss)
    low_called  <- sum(lo_nonmiss)

    high_ref <- 2 * high_called - high_alt
    low_ref  <- 2 * low_called  - low_alt

    tab <- matrix(c(high_alt, high_ref, low_alt, low_ref), nrow = 2, byrow = TRUE)

    if (any(is.na(tab)) || any(tab < 0) || sum(tab) == 0) next

    ft <- fisher.test(tab)
    pvals[j] <- ft$p.value
    odds[j]  <- unname(ft$estimate)
  }

  list(pvals = pvals, odds = odds)
}

obs <- run_one_scan(grp$group, geno_mat)

obs_df <- data.frame(
  CHROM = map$CHROM,
  POS = map$POS,
  SNP = map$varid,
  odds_ratio = obs$odds,
  p_value = obs$pvals,
  stringsAsFactors = FALSE
)

write.table(
  obs_df,
  file = file.path(out_dir, "chr01_pruned_observed_fisher.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

set.seed(1)
perm_min_p <- numeric(nperm)

n_high <- sum(grp$group == "high")
n_low  <- sum(grp$group == "low")

for (b in seq_len(nperm)) {
  perm_group <- sample(c(rep("high", n_high), rep("low", n_low)))
  ans <- run_one_scan(perm_group, geno_mat)
  perm_min_p[b] <- min(ans$pvals, na.rm = TRUE)
}

perm_df <- data.frame(perm = seq_len(nperm), min_p = perm_min_p)
write.table(
  perm_df,
  file = file.path(out_dir, "chr01_perm_min_pvalues.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

thr_p <- as.numeric(quantile(perm_min_p, 0.05, na.rm = TRUE))
thr_log10 <- -log10(thr_p)

writeLines(
  c(
    paste0("nperm: ", nperm),
    paste0("threshold_p: ", format(thr_p, scientific = TRUE)),
    paste0("threshold_log10: ", format(thr_log10, scientific = FALSE))
  ),
  con = file.path(out_dir, "chr01_perm_threshold.txt")
)


nano run_chr01_permutation_saqi.slurm
#!/usr/bin/env bash
#SBATCH --job-name=chr01_perm_saqi
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/chr01_perm_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/chr01_perm_%j.err

set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
WORK="/ibex/project/c2042/user_masrir/datepalm_metaxenia"

mkdir -p "${BASE}/logs" "${BASE}/permutation_chr01"

module load R

RAW_FILE="${BASE}/pruned_chr01/chr01_info30_pruned_A.raw"
GROUP_FILE="${BASE}/sample_groups_saqi.tsv"
OUT_DIR="${BASE}/permutation_chr01"
NPERM=100

test -s "${RAW_FILE}"
test -s "${GROUP_FILE}"
test -f "${WORK}/permute_saqi_chr01_fisher.R"

Rscript --vanilla "${WORK}/permute_saqi_chr01_fisher.R" \
  "${RAW_FILE}" \
  "${GROUP_FILE}" \
  "${OUT_DIR}" \
  "${NPERM}"
  
nano plot_saqi_whole_genome_manhattan.R
args <- commandArgs(trailingOnly = TRUE)

infile  <- args[1]
thrfile <- args[2]
outfile <- args[3]

libdir <- "/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/Rlib"
dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(libdir, .libPaths()))

if (!requireNamespace("qqman", quietly = TRUE)) {
  install.packages("qqman", repos = "https://cloud.r-project.org", lib = libdir)
}
library(qqman)

x <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("CHROM", "POS", "p_value") %in% names(x))) {
  stop("Input file must contain CHROM, POS, and p_value columns.")
}

thr_lines <- readLines(thrfile)
thr_p <- as.numeric(sub("threshold_p:\\s*", "", thr_lines[grep("^threshold_p:", thr_lines)]))
if (length(thr_p) != 1 || is.na(thr_p) || thr_p <= 0 || thr_p >= 1) {
  stop("Could not read threshold_p from threshold file.")
}

if (!("SNP" %in% names(x))) {
  x$SNP <- paste0(x$CHROM, ":", x$POS)
}

chr_order <- c(
  "Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10",
  "Chr11","Chr12","Chr13","Chr14","Chr14_male.SDR-oriented","Chr15","Chr16","Chr17","Chr18"
)

x <- x[!is.na(x$p_value) & is.finite(x$p_value) & x$p_value > 0 & x$p_value <= 1, , drop = FALSE]
x <- x[x$CHROM %in% chr_order, , drop = FALSE]

x$CHROM <- factor(x$CHROM, levels = chr_order, ordered = TRUE)
x$BP <- as.integer(x$POS)
x$P <- as.numeric(x$p_value)
x <- x[order(x$CHROM, x$BP), , drop = FALSE]
x$CHR <- as.integer(x$CHROM)

if (!("SNP" %in% names(x))) {
  x$SNP <- paste0(x$CHROM, ":", x$POS)
}

emp_line <- -log10(thr_p)
max_y <- max(-log10(x$P), emp_line, na.rm = TRUE) + 0.5

png(outfile, width = 3600, height = 1600, res = 250)
par(mar = c(5, 5, 4, 2) + 0.1)

manhattan(
  x,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  col = c("gray30", "gray70"),
  chrlabs = chr_order,
  suggestiveline = FALSE,
  genomewideline = FALSE,
  ylim = c(0, max_y),
  main = "Saqi whole-genome Manhattan plot"
)

abline(h = emp_line, col = "red", lwd = 2, lty = 2)

dev.off()

nano plot_saqi_whole_genome_manhattan.slurm

#!/usr/bin/env bash
#SBATCH --job-name=saqi_manhattan
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/saqi_manhattan_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/logs/saqi_manhattan_%j.err

set -euo pipefail

BASE="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia"
WORK="/ibex/project/c2042/user_masrir/datepalm_metaxenia"

mkdir -p "${BASE}/logs" "${BASE}/plots"

module load R

INPUT_TSV="${BASE}/Saqi_whole_genome.fisher_results.tsv"
THR_FILE="${BASE}/permutation_chr01/chr01_perm_threshold.txt"
OUT_PNG="${BASE}/plots/Saqi_whole_genome_empirical_manhattan.png"

test -s "${INPUT_TSV}"
test -s "${THR_FILE}"
test -f "${WORK}/plot_saqi_whole_genome_manhattan.R"

Rscript --vanilla "${WORK}/plot_saqi_whole_genome_manhattan.R" \
  "${INPUT_TSV}" \
  "${THR_FILE}" \
  "${OUT_PNG}"
  
bash prep_saqi_chr01_permutation.sh
sbatch run_chr01_permutation_saqi.slurm
sbatch plot_saqi_whole_genome_manhattan.slurm
