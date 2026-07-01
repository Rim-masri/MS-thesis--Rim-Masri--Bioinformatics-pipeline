# Create project structure
cd /ibex/scratch/projects/c2042
mkdir -p user_masrir/datepalm_metaxenia/{results,data}
cd user_masrir/datepalm_metaxenia

mkdir -p results/{mapping,vcf_call}
mkdir -p data/{WGS_data,WGS_filtered_data,reference_genome}

# Copy raw WGS data from sequencing delivery to WGS_data
cd /ibex/scratch/projects/c2042/celiim/KAUST_BCL/NovaSeqX/20251123_LL00134_0007_A232HNKLT4/Lane0/version_01
cp * /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/data/WGS_data/

# FastQC on raw reads (example sample)
cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/data/WGS_data
module load fastqc
sbatch -t 00:30:00 --job-name=fastqc --output=slurm_temp.txt --mem=10GB -N 1 -c 20 \
  --wrap "module load fastqc; fastqc -o fastqc_output/ -t 20 \
  metaxenia_..._R1_001.fastq.gz metaxenia_..._R2_001.fastq.gz"

# Generate Trim Galore adapter-trimming jobs for all non-BLANK libraries
ls -lhtr WGS_data/*_R1_001.fastq.gz | grep -v "BLANK" \
  | awk '{print "sbatch -t 10:30:00 --job-name=trim --output=slurm_temp.txt --mem=10GB -N 1 -c 30 --wrap \"module load trimgalore; trim_galore --fastqc -j 30 --output_dir WGS_filtered_data --paired "$9" "$9"2\""}' \
  | sed 's/R1_001.fastq.gz2/R2_001.fastq.gz"/g' \
  > adapter_trim.sh

bash adapter_trim.sh

# Move reference genome into project and inspect chromosomes
mv Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta \
   /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/data/reference_genome/

cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/data/reference_genome
grep ">" Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta

# Index reference for alignment
module load hisat2
hisat2-build Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta \
             Ajwa_hap1_with_organelles.with_Male_SDR-oriented