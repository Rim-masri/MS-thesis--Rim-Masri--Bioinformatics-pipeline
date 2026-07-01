# Sample lists
cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results

bcftools query -l ALL_samples_V3.Genotype.scatter300.vcf.gz > all_samples.txt

grep '^khudri_high_'   all_samples.txt > Khudri_high_samples.txt
grep '^khudri_low_'    all_samples.txt > Khudri_low_samples.txt
grep '^majdool_high_'  all_samples.txt > Majdool_high_samples.txt
grep '^majdool_low_'   all_samples.txt > Majdool_low_samples.txt
grep '^saqi_high_'     all_samples.txt > Saqi_high_samples.txt
grep '^saqi_low_'      all_samples.txt > Saqi_low_samples.txt

nano run_bcftools1_all.sh
#!/bin/bash
set -euo pipefail

module load bcftools

IN_ALL=ALL_samples_V3.Genotype.scatter300.vcf.gz

CHRS=(
  Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
  Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 ChrC ChrM Chr14_male.SDR-oriented
)

KHIGH_SAMP=Khudri_high_samples.txt
KLOW_SAMP=Khudri_low_samples.txt
MHIGH_SAMP=Majdool_high_samples.txt
MLOW_SAMP=Majdool_low_samples.txt
SHIGH_SAMP=Saqi_high_samples.txt
SLOW_SAMP=Saqi_low_samples.txt

DP_MIN=176
DP_MAX=17640

for CHR in "${CHRS[@]}"; do
  echo "=== Processing ${CHR} ==="

  OUT_CHR_SNP="ALL_samples_V3.Genotype.scatter300.${CHR}.biallelicSnps.vcf.gz"
  bcftools view \
    -r "${CHR}" \
    -v snps \
    -m2 -M2 \
    -Oz -o "${OUT_CHR_SNP}" "${IN_ALL}"
  bcftools index -f "${OUT_CHR_SNP}"

  OUT_DP="ALL_samples_V3.Genotyped.${CHR}.biallelicSnps.DP${DP_MIN}_${DP_MAX}.vcf.gz"
  bcftools view \
    -i "INFO/DP>=${DP_MIN} && INFO/DP<=${DP_MAX}" \
    -Oz -o "${OUT_DP}" "${OUT_CHR_SNP}"
  bcftools index -f "${OUT_DP}"

  KHIGH_VCF="KhudriHigh.${CHR}.DP${DP_MIN}_${DP_MAX}.vcf.gz"
  KLOW_VCF="KhudriLow.${CHR}.DP${DP_MIN}_${DP_MAX}.vcf.gz"
  MHIGH_VCF="MajdoolHigh.${CHR}.DP${DP_MIN}_${DP_MAX}.vcf.gz"
  MLOW_VCF="MajdoolLow.${CHR}.DP${DP_MIN}_${DP_MAX}.vcf.gz"
  SHIGH_VCF="SaqiHigh.${CHR}.DP${DP_MIN}_${DP_MAX}.vcf.gz"
  SLOW_VCF="SaqiLow.${CHR}.DP${DP_MIN}_${DP_MAX}.vcf.gz"

  bcftools view -S "${KHIGH_SAMP}" -Oz -o "${KHIGH_VCF}" "${OUT_DP}"
  bcftools index -f "${KHIGH_VCF}"

  bcftools view -S "${KLOW_SAMP}" -Oz -o "${KLOW_VCF}" "${OUT_DP}"
  bcftools index -f "${KLOW_VCF}"

  bcftools view -S "${MHIGH_SAMP}" -Oz -o "${MHIGH_VCF}" "${OUT_DP}"
  bcftools index -f "${MHIGH_VCF}"

  bcftools view -S "${MLOW_SAMP}" -Oz -o "${MLOW_VCF}" "${OUT_DP}"
  bcftools index -f "${MLOW_VCF}"

  bcftools view -S "${SHIGH_SAMP}" -Oz -o "${SHIGH_VCF}" "${OUT_DP}"
  bcftools index -f "${SHIGH_VCF}"

  bcftools view -S "${SLOW_SAMP}" -Oz -o "${SLOW_VCF}" "${OUT_DP}"
  bcftools index -f "${SLOW_VCF}"

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "${KHIGH_VCF}" > "KhudriHigh.${CHR}.AD.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "${KLOW_VCF}"  > "KhudriLow.${CHR}.AD.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "${MHIGH_VCF}" > "MajdoolHigh.${CHR}.AD.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "${MLOW_VCF}"  > "MajdoolLow.${CHR}.AD.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "${SHIGH_VCF}" > "SaqiHigh.${CHR}.AD.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "${SLOW_VCF}"  > "SaqiLow.${CHR}.AD.tsv"
done

echo "All chromosomes and cultivars processed."

nano bcftools1_all.slurm
#!/bin/bash
#SBATCH --job-name=bcf_all_chr
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --output=bcf_all_chr.%j.out
#SBATCH --error=bcf_all_chr.%j.err

set -euo pipefail

echo "Job started on $(date)"
echo "Running on host $(hostname)"

cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results
bash run_bcftools1_all.sh

echo "Job finished on $(date)"

nano make_wholepop_AD.sbatch
#!/bin/bash
#SBATCH --job-name=WG_allPop_AD
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --output=WG_allPop_AD.%j.out
#SBATCH --error=WG_allPop_AD.%j.err

set -euo pipefail

module load bcftools

cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results

IN_ALL=ALL_samples_V3.Genotype.scatter300.vcf.gz
DP_MIN=176
DP_MAX=17640

OUT_VCF="ALL_samples_V3.Genotyped.biallelicSnps.DP${DP_MIN}_${DP_MAX}.vcf.gz"
OUT_TSV="WholeGenome_allPop_AD.tsv"

bcftools view \
  -v snps \
  -m2 -M2 \
  -i "INFO/DP>=${DP_MIN} && INFO/DP<=${DP_MAX}" \
  -Oz -o "${OUT_VCF}" \
  "${IN_ALL}"

bcftools index -f "${OUT_VCF}"

{
  printf "CHROM\tPOS\tREF\tALT"
  bcftools query -l "${OUT_VCF}" | awk '{printf "\t"$1}'
  printf "\n"
} > "${OUT_TSV}"

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
  "${OUT_VCF}" >> "${OUT_TSV}"
  
cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results

chmod +x run_bcftools1_all.sh

sbatch bcftools1_all.slurm
sbatch make_wholepop_AD.sbatch

