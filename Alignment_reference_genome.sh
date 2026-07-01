cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia
ref="data/reference_genome/hisat2_index/Ajwa_hap1_with_organelles.with_Male_SDR-oriented"

# Pairing R1/R2 and submitting HISAT2 jobs
for read1 in data/WGS_filtered_data/*R1_001_val_1.fq.gz; do
  read2=$(echo "$read1" | sed 's/R1_001_val_1/R2_001_val_2/g')
  sort_bam=$(echo "$read1" | sed 's/R1_001_val_1.fq.gz/sorted.bam/g; s|data/WGS_filtered_data|results/mapping2|g')
  bam=$(echo "$read1" | sed 's/_R1_001_val_1.fq.gz/.bam/g; s|data/WGS_filtered_data|results/mapping2|g')
  base_name=$(echo "$read1" | cut -d"/" -f3)

  sbatch -t 1-10:30:00 --job-name="$base_name" \
    --output=hisat2_out/$base_name.out -e hisat2_out/$base_name.err \
    --mem=100GB -N 1 -c 30 \
    --export=ref="$ref",read1="$read1",read2="$read2",bam="$bam",base_name="$base_name" \
    --wrap "module load hisat2 samtools; hisat2 -x \$ref -1 \$read1 -2 \$read2 -p 30 | samtools view -@ 20 -b -h > \$bam"
done

#Duplicate removal
cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results/mapping2
mkdir /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results/mapped_dup_rm

module load samtools

for file in *L000.bam; do
  sbatch -t 1-10:30:00 --job-name="$file" \
    --output=align_stat/$file.out -e align_stat/$file.err \
    --mem=50GB -c 16 --export=file="$file" \
    --wrap "module load samtools; \
      samtools sort -n \$file \
        | samtools fixmate -m - - \
        | samtools sort -@ 20 -o \$file.fix.sort.bam; \
      samtools index \$file.fix.sort.bam; \
      samtools markdup -r \$file.fix.sort.bam \
        ../mapped_dup_rm/\${file}_bwa.fix.sort.dup_rm.bam \
        -f ../mapped_dup_rm/\${file}.dup_rm.txt"
done

#Sample renaming
cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia
module load R/4.0.0/gnu-12.2.0
R

csv <- read.csv("Book1.csv")
csv$new_name <- paste(csv$cultivar, csv$quality, csv$serial, sep = "_")
csv$new_name <- paste0(csv$new_name, ".bam")

all.files <- dir("results/mapped_dup_rm")
bam.files <- grep("sorted.bam_bwa.fix.sort.dup_rm.bam", all.files, value = TRUE)

bam_table <- data.frame(bam.file = bam.files, stringsAsFactors = FALSE)
bam_splittable <- strsplit(bam_table$bam.file, "_")
bam_table$UUID <- sapply(bam_splittable, function(x) x[2])

new_csv <- merge(csv, bam_table, by = "UUID")

# Copy to new names
for (i in seq_len(nrow(new_csv))) {
  old <- file.path("results/mapped_dup_rm", new_csv$bam.file[i])
  new <- file.path("results/mapped_dup_rm", new_csv$new_name[i])
  system(paste("cp", old, new))
}

#Flagstat and coverage statistics
cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results/mapped_dup_rm
module load samtools

# Flagstat per renamed BAM
mkdir renamed_flagstat
for file in *{high,low}*bam; do
  sbatch -t 1-10:30:00 --job-name="$file" \
    --output=renamed_flagstat/$file.out -e renamed_flagstat/$file.err \
    --mem=50GB -c 16 --export=file="$file" \
    --wrap "module load samtools; samtools flagstat \$file > renamed_flagstat/\$file.flagstat.txt"
done

# Coverage per renamed BAM
mkdir coverage
for file in *{high,low}*bam; do
  sbatch -t 1-10:30:00 --job-name="$file" \
    --output=coverage/$file.out -e coverage/$file.err \
    --mem=50GB -c 16 --export=file="$file" \
    --wrap 'module load samtools; samtools depth "$file" | awk "{sum+=\$3} END {print sum/NR}" > coverage/"$file".depth.txt'
done

#GATK per-sample gVCFS
cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia
module load samtools picard gatk

# Reference indices for GATK
samtools faidx data/reference_genome/Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta
picard CreateSequenceDictionary \
  R=data/reference_genome/Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta \
  O=data/reference_genome/Ajwa_hap1_with_organelles.with_Male_SDR-oriented.dict

mkdir results/GATK_gvcf

# Add read groups + HaplotypeCaller
k=0
for file in results/mapped_dup_rm/*{high,low}*.bam; do
  k=$((k+1))
  sample=$(basename "$file" | sed 's/.bam//g')
  file_out=$(echo "$file" | sed 's/.bam/.RG.bam/g')

  sbatch -t 12:00:00 --job-name="$sample" \
    -o results/GATK_gvcf/${sample}.out \
    -e results/GATK_gvcf/${sample}.err \
    --mem=50GB -c 16 \
    --export=file="$file",sample="$sample",file_out="$file_out" \
    --wrap "module load gatk/4.6.2.0 picard samtools; \
      picard AddOrReplaceReadGroups I=\$file O=\$file_out RGID=\$sample RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=\$sample; \
      samtools index \$file_out; \
      gatk --java-options '-Xmx4g' HaplotypeCaller \
        -R data/reference_genome/Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta \
        -I \$file_out \
        -O results/GATK_gvcf/\$sample.g.vcf.gz \
        -ERC GVCF"
done

#Combine gVCFs + Genotype gVCFS

cd /ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/results
ref="/ibex/scratch/projects/c2042/user_masrir/datepalm_metaxenia/data/reference_genome/Ajwa_hap1_with_organelles.with_Male_SDR-oriented.fasta"

module load gatk/4.6.2.0

# SplitIntervals for 300 genome chunks
sbatch -t 12:00:00 --job-name="SplitIntervals" \
  -o SplitIntervals.v2.out -e SplitIntervals.v2.err \
  --mem=50GB -c 32 --export=ref="$ref" \
  --wrap 'gatk SplitIntervals -R "$ref" --scatter-count 300 -O intervals_300'

# GenotypeGVCFs per interval (using ALL_samples_V3.CombineGVCFs.g.vcf.gz)
mkdir -p gg_300 logs3
N=$(ls intervals_300/*-scattered.interval_list | wc -l)

sbatch -t 2-12:00:00 --job-name="V3.GG.scatter300" \
  -o logs3/GG_%A_%a.out -e logs3/GG_%A_%a.err \
  --mem=80G -c 8 --array=0-$(($N-1)) --export=ALL,ref="$ref" \
  --wrap 'INT=$(ls intervals_300/*-scattered.interval_list | sed -n "$((SLURM_ARRAY_TASK_ID+1))p"); \
          OUT=gg_300/$(basename "$INT" .interval_list).vcf.gz; \
          gatk --java-options "-Xmx75g" GenotypeGVCFs \
            -R "$ref" \
            -V ALL_samples_V3.CombineGVCFs.g.vcf.gz \
            -L "$INT" \
            -O "$OUT" \
            --max-alternate-alleles 20'

# GatherVcfs to produce final genotype VCF
ls gg_300/*.vcf.gz | sort > gg_parts_300.list

sbatch -t 1-12:00:00 --job-name="v2.GatherVcfs" \
  -o GatherVcfs.v2.out -e GatherVcfs.v2.err \
  --mem=200GB -c 64 --export=ref="$ref" \
  --wrap 'module load gatk/4.6.2.0 tabix; \
          gatk GatherVcfs -I gg_parts_300.list -O ALL_samples_V3.Genotype.scatter300.vcf.gz; \
          tabix -p vcf ALL_samples_V3.Genotype.scatter300.vcf.gz'
          

