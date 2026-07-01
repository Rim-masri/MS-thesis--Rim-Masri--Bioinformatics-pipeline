#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch
BAMDIR=/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results/mapped_dup_rm
OUT=${WD}/lists/list_of_bams_for_stitch_saqi_no_header.tsv

mkdir -p "${WD}/lists"

find "${BAMDIR}" -maxdepth 1 -type f \( -name '*saqi*high*.RG.bam' -o -name '*saqi*low*.RG.bam' \) \
  | sort -u > "${OUT}"

echo "Saqi BAM list written to: ${OUT}"
wc -l "${OUT}"
head "${OUT}"

#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch
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

  if [ -z "${L}" ]; then
    echo "ERROR: ${CHR} not found in ${FAI}"
    exit 1
  fi

  awk -v C="${CHR}" -v L="${L}" 'BEGIN{
    chunk=5000000; id=1;
    for (s=1; s<=L; s+=chunk) {
      e=s+chunk-1;
      if (e>L) e=L;
      printf "%s\t%d\t%d\tchunk%03d\n", C, s, e, id++;
    }
  }' > "${WD}/lists/${CHR}.chunks.5Mb.tsv"

  echo "POS rows for ${CHR}:"
  wc -l "${WD}/lists/${CHR}.pos.tsv"
  echo "Chunks for ${CHR}:"
  wc -l "${WD}/lists/${CHR}.chunks.5Mb.tsv"
done

#!/usr/bin/env bash
#SBATCH --job-name=stitch_saqi
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/logs/stitch_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/logs/stitch_%A_%a.err
#SBATCH --export=ALL

set -euo pipefail

module load bcftools
module load samtools
module load tabix
module load R
module load stitch/1.8.4

WD=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch
CHR=${CHR:?ERROR: CHR not set}
BAMLIST=${WD}/lists/list_of_bams_for_stitch_saqi_no_header.tsv
POSFILE=${WD}/lists/${CHR}.pos.tsv
CHUNKFILE=${WD}/lists/${CHR}.chunks.5Mb.tsv

read CHR_FROM_FILE START END CHUNK_ID < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHUNKFILE}")

OUTDIR=${WD}/stitch/${CHR}/${CHUNK_ID}_${START}_${END}
TMPDIR=${WD}/stitch/tmp_${CHR}_${CHUNK_ID}_${START}_${END}

mkdir -p "${OUTDIR}" "${TMPDIR}" "${WD}/logs"

echo "Starting STITCH"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "CHR env: ${CHR}"
echo "CHR file: ${CHR_FROM_FILE}"
echo "START: ${START}"
echo "END: ${END}"
echo "CHUNK_ID: ${CHUNK_ID}"
echo "BAMLIST: ${BAMLIST}"
echo "POSFILE: ${POSFILE}"
echo "OUTDIR: ${OUTDIR}"

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

VCF_OUT=$(find "${OUTDIR}" -type f \( -name "*.vcf.gz" -o -name "*.vcf" \) | sort | head -n 1)

if [ -z "${VCF_OUT}" ]; then
  echo "ERROR: no VCF output found in ${OUTDIR}"
  ls -lhR "${OUTDIR}"
  exit 1
fi

echo "Found VCF output: ${VCF_OUT}"

if [[ "${VCF_OUT}" == *.vcf ]]; then
  bgzip -f "${VCF_OUT}"
  VCF_OUT="${VCF_OUT}.gz"
fi

bcftools index -f -t "${VCF_OUT}"
bcftools stats "${VCF_OUT}" > "${OUTDIR}/stitch.${CHR}.${START}_${END}.stats.txt"

rm -rf "${TMPDIR}"

echo "Finished STITCH"
echo "Output VCF: ${VCF_OUT}"

#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch
cd "${WD}"

CHRS=(
  Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
  Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 ChrC ChrM Chr14_male.SDR-oriented
)

for CHR in "${CHRS[@]}"; do
  N=$(wc -l < "lists/${CHR}.chunks.5Mb.tsv")
  echo "Submitting ${CHR} with ${N} chunks"
  sbatch --array=1-"${N}" --export=ALL,CHR="${CHR}" stitch_saqi.slurm
done

module load bcftools
bcftools query -l /ibex/scratch/projects/c2042/user_masrir/saqi_stitch/stitch/Chr01/chunk001_1_5000000/stitch.Chr01.1.5000000.vcf.gz > sample_order.txt

awk 'BEGIN{print "sample_id\tgroup"} {
  if ($0 ~ /high/) print $0 "\thigh";
  else if ($0 ~ /low/) print $0 "\tlow";
}' sample_order.txt > sample_groups.tsv

#!/usr/bin/env bash
set -euo pipefail

WD=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch
module load bcftools

CHRS=(
  Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
  Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 ChrC ChrM Chr14_male.SDR-oriented
)

mkdir -p "${WD}/dosage_chunks"

for CHR in "${CHRS[@]}"; do
  mkdir -p "${WD}/dosage_chunks/${CHR}"

  for vcf in "${WD}"/stitch/${CHR}/*/*.vcf.gz; do
    [ -e "$vcf" ] || continue
    base=$(basename "$vcf" .vcf.gz)
    bcftools query -f '%CHROM\t%POS[\t%DS]\n' "$vcf" > "${WD}/dosage_chunks/${CHR}/${base}.dosage.tsv"
  done

  cat "${WD}/dosage_chunks/${CHR}"/*.dosage.tsv > "${WD}/${CHR}.dosage.tsv"

  sort -k1,1 -k2,2n "${WD}/${CHR}.dosage.tsv" \
    | awk 'BEGIN{FS=OFS="\t"} {key=$1 FS $2; if (!seen[key]++) print}' \
    > "${WD}/${CHR}.dedup.dosage.tsv"
done

cat \
  "${WD}/Chr01.dedup.dosage.tsv" "${WD}/Chr02.dedup.dosage.tsv" "${WD}/Chr03.dedup.dosage.tsv" "${WD}/Chr04.dedup.dosage.tsv" \
  "${WD}/Chr05.dedup.dosage.tsv" "${WD}/Chr06.dedup.dosage.tsv" "${WD}/Chr07.dedup.dosage.tsv" "${WD}/Chr08.dedup.dosage.tsv" \
  "${WD}/Chr09.dedup.dosage.tsv" "${WD}/Chr10.dedup.dosage.tsv" "${WD}/Chr11.dedup.dosage.tsv" "${WD}/Chr12.dedup.dosage.tsv" \
  "${WD}/Chr13.dedup.dosage.tsv" "${WD}/Chr14.dedup.dosage.tsv" "${WD}/Chr15.dedup.dosage.tsv" "${WD}/Chr16.dedup.dosage.tsv" \
  "${WD}/Chr17.dedup.dosage.tsv" "${WD}/Chr18.dedup.dosage.tsv" "${WD}/ChrC.dedup.dosage.tsv" "${WD}/ChrM.dedup.dosage.tsv" \
  "${WD}/Chr14_male.SDR-oriented.dedup.dosage.tsv" \
  > "${WD}/Saqi_whole_genome.dedup.dosage.tsv"
  
dos <- read.table("Saqi_whole_genome.dedup.dosage.tsv",
                  header = FALSE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

grp <- read.table("sample_groups.txt",
                  header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)

high_n <- sum(grp$group == "high")
low_n  <- sum(grp$group == "low")

cat("high_n =", high_n, "\n")
cat("low_n =", low_n, "\n")
cat("ncol(dos) =", ncol(dos), "\n")

stopifnot(ncol(dos) == 2 + high_n + low_n)

chrom <- dos[[1]]
pos   <- dos[[2]]

mat <- as.matrix(dos[, -c(1,2)])
storage.mode(mat) <- "numeric"

high_mat <- mat[, 1:high_n, drop = FALSE]
low_mat  <- mat[, (high_n + 1):(high_n + low_n), drop = FALSE]

high_alt <- rowSums(high_mat, na.rm = TRUE)
low_alt  <- rowSums(low_mat, na.rm = TRUE)

high_ref <- 2 * high_n - high_alt
low_ref  <- 2 * low_n  - low_alt

pvals <- numeric(nrow(mat))
odds  <- numeric(nrow(mat))

for (i in seq_len(nrow(mat))) {
  tab <- matrix(c(high_alt[i], high_ref[i],
                  low_alt[i],  low_ref[i]),
                nrow = 2, byrow = TRUE)
  ft <- fisher.test(tab)
  pvals[i] <- ft$p.value
  odds[i]  <- ifelse(is.null(ft$estimate), NA, unname(ft$estimate))
}

out <- data.frame(
  CHROM = chrom,
  POS = pos,
  high_ALT_sum = high_alt,
  high_REF_sum = high_ref,
  low_ALT_sum = low_alt,
  low_REF_sum = low_ref,
  odds_ratio = odds,
  p_value = pvals,
  FDR = p.adjust(pvals, method = "fdr")
)

write.table(out,
            file = "Saqi_whole_genome.fisher_results.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
            
#!/usr/bin/env bash
#SBATCH --job-name=fisher_whole_genome_saqi
#SBATCH --cpus-per-task=1
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/logs/fisher_whole_genome_saqi_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/logs/fisher_whole_genome_saqi_%j.err

set -euo pipefail

module load R

WD=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch
cd "${WD}"

Rscript --vanilla fisher_whole_genome_saqi.R

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

library(ggplot2)
library(dplyr)
library(readr)

infile <- "/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/Saqi_whole_genome.fisher_results.tsv"
outfile <- "/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/plots/Saqi_whole_genome_manhattan.png"

dir.create("/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/plots",
           showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(infile))

df <- read_tsv(infile, show_col_types = FALSE)

keep_order <- c(
  paste0("Chr", sprintf("%02d", 1:14)),
  "Chr14_male.SDR-oriented",
  paste0("Chr", sprintf("%02d", 15:18))
)

df <- df %>%
  filter(CHROM %in% keep_order) %>%
  mutate(
    CHROM = factor(CHROM, levels = keep_order),
    neglogp = -log10(p_value)
  ) %>%
  arrange(CHROM, POS)

stopifnot(nrow(df) > 0)

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

n_markers <- nrow(df)
bonf_threshold <- 0.05 / n_markers

p <- ggplot(df, aes(x = bp_cum, y = neglogp, color = factor(color_id))) +
  geom_point(size = 0.35, alpha = 0.8) +
  geom_hline(
    yintercept = -log10(bonf_threshold),
    linetype = "dashed",
    color = "red3",
    linewidth = 0.7
  ) +
  scale_color_manual(values = c("0" = "gray20", "1" = "gray70")) +
  scale_x_continuous(
    breaks = chr_info$center,
    labels = chr_info$CHROM,
    expand = c(0.01, 0)
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "Saqi whole-genome Manhattan plot"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

ggsave(
  filename = outfile,
  plot = p,
  width = 14,
  height = 6,
  dpi = 300
)

#!/usr/bin/env bash
#SBATCH --job-name=manhattan_saqi
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/logs/manhattan_saqi_%j.out
#SBATCH --error=/ibex/scratch/projects/c2042/user_masrir/saqi_stitch/logs/manhattan_saqi_%j.err

set -euo pipefail
module load R

export R_LIBS_USER=/ibex/scratch/projects/c2042/user_masrir/Rlibs/4.5.0
mkdir -p "$R_LIBS_USER"

cd /ibex/scratch/projects/c2042/user_masrir/saqi_stitch
mkdir -p logs plots

Rscript --vanilla -e 'pkgs <- c("ggplot2","dplyr","readr"); for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org", lib = Sys.getenv("R_LIBS_USER"))'

Rscript --vanilla Manhattan_Saqi_imputed.R

