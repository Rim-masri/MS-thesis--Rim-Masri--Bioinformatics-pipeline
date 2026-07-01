cd /ibex/scratch/projects/c2042/user_masrir/majdool_stitch/

nano build_majdool_bam_list.sh
#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch
BAMDIR=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results/mapped_dup_rm
OUT=${WD}/lists/list_of_bams_for_stitch_majdool_no_header.tsv

mkdir -p "${WD}/lists"

find "${BAMDIR}" -maxdepth 1 -type f \
  \( -name '*majdool*high*.RG.bam' -o -name '*majdool*low*.RG.bam' \) \
  | sed 's/\.bam\.RG\.bam$/.RG.bam/' \
  | sort -u > "${OUT}"

echo "Wrote: ${OUT}"
wc -l "${OUT}"
head "${OUT}"

nano prepare_stitch_inputs_majdool.sh

#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch
VCF=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results/ALL_samples_V3.Genotype.scatter300.vcf.gz
REF=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/data/reference_genome/Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta
FAI=${REF}.fai

CHRS=(
  Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
  Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 ChrC ChrM Chr14_male.SDR-oriented
)

module load samtools
module load bcftools

mkdir -p "${WD}/lists" "${WD}/logs"
samtools faidx "${REF}"

for CHR in "${CHRS[@]}"; do
  echo "Preparing ${CHR}"

  bcftools view -r "${CHR}" -v snps -m2 -M2 "${VCF}" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    > "${WD}/lists/${CHR}.pos.tsv"

  L=$(awk -v c="${CHR}" '$1==c{print $2}' "${FAI}")
  test -n "${L}"

  awk -v C="${CHR}" -v L="${L}" 'BEGIN{
    chunk=5000000; id=1;
    for (s=1; s<=L; s+=chunk) {
      e=s+chunk-1; if (e>L) e=L;
      printf "%s\t%d\t%d\tchunk%03d\n", C, s, e, id++;
    }
  }' > "${WD}/lists/${CHR}.chunks.5Mb.tsv"

  echo "${CHR} POS:"
  wc -l "${WD}/lists/${CHR}.pos.tsv"
  echo "${CHR} chunks:"
  wc -l "${WD}/lists/${CHR}.chunks.5Mb.tsv"
done

nano stitch_majdool.slurm
#!/usr/bin/env bash
#SBATCH --job-name=stitch_majdool
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/logs/stitch_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/logs/stitch_%A_%a.err
#SBATCH --export=ALL

set -euo pipefail

module load bcftools
module load samtools
module load tabix
module load R
module load stitch/1.8.4

WD=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch
CHR=${CHR:?ERROR: CHR not set}

BAMLIST=${WD}/lists/list_of_bams_for_stitch_majdool_no_header.tsv
POSFILE=${WD}/lists/${CHR}.pos.tsv
CHUNKFILE=${WD}/lists/${CHR}.chunks.5Mb.tsv

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHUNKFILE}" || true)
test -n "${LINE}"

read -r CHR_FROM_FILE START END CHUNK_ID <<< "${LINE}"
test "${CHR_FROM_FILE}" = "${CHR}"

OUTDIR=${WD}/stitch/${CHR}/${CHUNK_ID}_${START}_${END}
TMPDIR=${WD}/stitch/tmp_${CHR}_${CHUNK_ID}_${START}_${END}
mkdir -p "${OUTDIR}" "${TMPDIR}"

Rscript /ibex/sw/rl9c/stitch/1.8.4/rl9.4_conda/STITCH/STITCH.R \
  --chr="${CHR}" \
  --bamlist="${BAMLIST}" \
  --posfile="${POSFILE}" \
  --outputdir="${OUTDIR}" \
  --K=12 \
  --nGen=100 \
  --nCores="${SLURM_CPUS_PER_TASK}" \
  --method=diploid \
  --regionStart="${START}" \
  --regionEnd="${END}" \
  --buffer=250000 \
  --inputBundleBlockSize=150 \
  --output_format=bgvcf \
  --iSizeUpperLimit=1000 \
  --tempdir="${TMPDIR}"

VCF_OUT=$(find "${OUTDIR}" -type f -name "*.vcf.gz" | sort | head -n 1 || true)
test -n "${VCF_OUT}"
test -s "${VCF_OUT}"

bcftools index -f -t "${VCF_OUT}"
bcftools stats "${VCF_OUT}" > "${OUTDIR}/stitch.${CHR}.${START}_${END}.stats.txt"
rm -rf "${TMPDIR}"

nano submit_stitch_all_majdool.sh
#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch
cd "${WD}"

CHRS=(
  Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
  Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 ChrC ChrM Chr14_male.SDR-oriented
)

for CHR in "${CHRS[@]}"; do
  N=$(wc -l < "lists/${CHR}.chunks.5Mb.tsv")
  echo "Submitting ${CHR} with ${N} chunks"
  sbatch --array=1-"${N}" --export=ALL,CHR="${CHR}" stitch_majdool.slurm
done

nano make_sample_groups_majdool.sh
#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch
module load bcftools

VCF=$(find "${WD}/stitch/Chr01" -type f -name "*.vcf.gz" | sort | head -n 1)
test -n "${VCF}"

bcftools query -l "${VCF}" > "${WD}/sample_order_majdool.txt"

awk 'BEGIN { print "sample_id\tgroup" }
{
  if ($0 ~ /high/i) print $0 "\thigh";
  else if ($0 ~ /low/i) print $0 "\tlow";
  else print $0 "\tNA";
}' "${WD}/sample_order_majdool.txt" > "${WD}/sample_groups.tsv"

wc -l "${WD}/sample_order_majdool.txt"
head "${WD}/sample_groups.tsv"

nano extract_dosage_majdool_all.sh
#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch
module load bcftools

CHRS=(
  Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
  Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 ChrC ChrM Chr14_male.SDR-oriented
)

mkdir -p "${WD}/dosage_chunks"

REF_VCF=$(find "${WD}/stitch/Chr01" -type f -name "*.vcf.gz" | sort | head -n 1)
test -n "${REF_VCF}"

bcftools query -l "${REF_VCF}" > "${WD}/sample_order_majdool.txt"
paste <(printf "CHROM\nPOS\n") /dev/null >/dev/null 2>&1 || true

{
  printf "CHROM\tPOS"
  awk '{printf "\t%s", $0} END{printf "\n"}' "${WD}/sample_order_majdool.txt"
} > "${WD}/Majdool_whole_genome.dedup.dosage.tsv"

for CHR in "${CHRS[@]}"; do
  mkdir -p "${WD}/dosage_chunks/${CHR}"
  rm -f "${WD}/dosage_chunks/${CHR}"/*.dosage.tsv

  for vcf in "${WD}"/stitch/${CHR}/*/*.vcf.gz; do
    [ -e "$vcf" ] || continue
    base=$(basename "$vcf" .vcf.gz)

    bcftools query \
      -S "${WD}/sample_order_majdool.txt" \
      -f '%CHROM\t%POS[\t%DS]\n' \
      "$vcf" > "${WD}/dosage_chunks/${CHR}/${base}.dosage.tsv"
  done

  cat "${WD}/dosage_chunks/${CHR}"/*.dosage.tsv > "${WD}/${CHR}.dosage.tsv"

  sort -k1,1 -k2,2n "${WD}/${CHR}.dosage.tsv" \
    | awk 'BEGIN{FS=OFS="\t"} {k=$1 FS $2; if (!seen[k]++) print}' \
    > "${WD}/${CHR}.dedup.dosage.tsv"

  cat "${WD}/${CHR}.dedup.dosage.tsv" >> "${WD}/Majdool_whole_genome.dedup.dosage.tsv"
done

nano fisher_whole_genome_majdool.R
#!/usr/bin/env Rscript

dos <- read.table(
  "Majdool_whole_genome.dedup.dosage.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

grp <- read.table(
  "sample_groups.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(all(c("sample_id", "group") %in% names(grp)))
stopifnot(all(c("CHROM", "POS") %in% names(dos)))

grp$sample_id <- trimws(grp$sample_id)
grp$group <- trimws(grp$group)

grp <- grp[grp$group %in% c("high", "low"), , drop = FALSE]

sample_cols <- setdiff(names(dos), c("CHROM", "POS"))

missing_in_dos <- setdiff(grp$sample_id, sample_cols)
if (length(missing_in_dos) > 0) {
  stop("These samples are in sample_groups.tsv but missing from dosage table: ",
       paste(missing_in_dos, collapse = ", "))
}

missing_in_grp <- setdiff(sample_cols, grp$sample_id)
if (length(missing_in_grp) > 0) {
  message("Ignoring dosage columns not present in sample_groups.tsv: ",
          paste(missing_in_grp, collapse = ", "))
}

grp <- grp[match(intersect(grp$sample_id, sample_cols), grp$sample_id), , drop = FALSE]
grp <- grp[grp$sample_id %in% sample_cols, , drop = FALSE]

ordered_samples <- c(grp$sample_id[grp$group == "high"],
                     grp$sample_id[grp$group == "low"])

dos2 <- dos[, c("CHROM", "POS", ordered_samples), drop = FALSE]

high_n <- sum(grp$group == "high")
low_n  <- sum(grp$group == "low")

mat <- as.matrix(dos2[, ordered_samples, drop = FALSE])
storage.mode(mat) <- "numeric"

high_mat <- mat[, seq_len(high_n), drop = FALSE]
low_mat  <- mat[, (high_n + 1):(high_n + low_n), drop = FALSE]

high_alt <- rowSums(high_mat, na.rm = TRUE)
low_alt  <- rowSums(low_mat, na.rm = TRUE)

high_called <- rowSums(!is.na(high_mat))
low_called  <- rowSums(!is.na(low_mat))

high_ref <- 2 * high_called - high_alt
low_ref  <- 2 * low_called  - low_alt

pvals <- rep(NA_real_, nrow(mat))
odds  <- rep(NA_real_, nrow(mat))

for (i in seq_len(nrow(mat))) {
  tab <- matrix(
    c(high_alt[i], high_ref[i],
      low_alt[i],  low_ref[i]),
    nrow = 2, byrow = TRUE
  )

  if (any(is.na(tab)) || any(tab < 0) || sum(tab) == 0) next

  ft <- fisher.test(tab)
  pvals[i] <- ft$p.value
  odds[i] <- unname(ft$estimate)
}

out <- data.frame(
  CHROM = dos2$CHROM,
  POS = dos2$POS,
  SNP = paste0(dos2$CHROM, ":", dos2$POS),
  high_ALT_sum = high_alt,
  high_REF_sum = high_ref,
  low_ALT_sum = low_alt,
  low_REF_sum = low_ref,
  odds_ratio = odds,
  p_value = pvals,
  FDR = p.adjust(pvals, method = "fdr"),
  stringsAsFactors = FALSE
)

write.table(
  out,
  file = "Majdool_whole_genome.fisher_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Wrote Majdool_whole_genome.fisher_results.tsv\n")
cat("high_n =", high_n, "\n")
cat("low_n =", low_n, "\n")
cat("variants =", nrow(out), "\n")

nano fisher_whole_genome_majdool.slurm
#!/usr/bin/env bash
#SBATCH --job-name=fisher_whole_genome_majdool
#SBATCH --cpus-per-task=1
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/logs/fisher_whole_genome_majdool_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/logs/fisher_whole_genome_majdool_%j.err

set -euo pipefail
module load R

cd /ibex/scratch/projects/c2042/user_masrir/majdool_stitch
Rscript --vanilla fisher_whole_genome_majdool.R

nano Manhattan_Majdool_imputed.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

library(ggplot2)
library(dplyr)
library(readr)

infile <- "/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/Majdool_whole_genome.fisher_results.tsv"
outfile <- "/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/plots/Majdool_whole_genome_imputed_manhattan.png"

df <- read_tsv(infile, show_col_types = FALSE)

keep_order <- c(
  paste0("Chr", sprintf("%02d", 1:14)),
  "Chr14_male.SDR-oriented",
  paste0("Chr", sprintf("%02d", 15:18))
)

df <- df %>%
  filter(CHROM %in% keep_order, !is.na(p_value), p_value > 0, p_value <= 1) %>%
  mutate(
    CHROM = factor(CHROM, levels = keep_order),
    neglogp = -log10(p_value)
  ) %>%
  arrange(CHROM, POS)

n_tests <- nrow(df)
bonf_p <- 0.05 / n_tests
bonf_line <- -log10(bonf_p)

chr_info <- df %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(POS, na.rm = TRUE), .groups = "drop") %>%
  arrange(factor(CHROM, levels = keep_order)) %>%
  mutate(
    start = lag(cumsum(chr_len), default = 0),
    center = start + chr_len / 2,
    color_id = row_number() %% 2
  )

df <- df %>%
  left_join(chr_info %>% select(CHROM, start, color_id), by = "CHROM") %>%
  mutate(bp_cum = POS + start)

p <- ggplot(df, aes(bp_cum, neglogp, color = factor(color_id))) +
  geom_point(size = 0.35, alpha = 0.8) +
  geom_hline(yintercept = bonf_line, linetype = "dashed", color = "red3", linewidth = 0.7) +
  scale_color_manual(values = c("0" = "gray20", "1" = "gray70")) +
  scale_x_continuous(breaks = chr_info$center, labels = chr_info$CHROM, expand = c(0.01, 0)) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "Majdool whole-genome Manhattan plot"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

ggsave(outfile, p, width = 14, height = 6, dpi = 300)

cat("Markers tested:", n_tests, "\n")
cat("Bonferroni p-value:", bonf_p, "\n")
cat("Bonferroni -log10 line:", bonf_line, "\n")
cat("Output written to:", outfile, "\n")

nano Manhattan_Majdool_imputed.slurm
#!/usr/bin/env bash
#SBATCH --job-name=manhattan_majdool
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/logs/manhattan_majdool_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/majdool_stitch/logs/manhattan_majdool_%j.err

set -euo pipefail
module load R

export R_LIBS_USER=/ibex/scratch/projects/c2042/user_masrir/Rlibs/4.5.0
mkdir -p "$R_LIBS_USER"

cd /ibex/scratch/projects/c2042/user_masrir/majdool_stitch
mkdir -p logs plots

Rscript --vanilla -e 'pkgs <- c("ggplot2","dplyr","readr"); for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org", lib = Sys.getenv("R_LIBS_USER"))'
Rscript --vanilla Manhattan_Majdool_imputed.R


bash build_majdool_bam_list.sh
bash prepare_stitch_inputs_majdool.sh
bash submit_stitch_all_majdool.sh
bash make_sample_groups_majdool.sh
bash extract_dosage_majdool_all.sh
sbatch fisher_whole_genome_majdool.slurm
sbatch Manhattan_Majdool_imputed.slurm