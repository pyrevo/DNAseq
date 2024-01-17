#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
set -x
# set current directory
#data_dir="$( cd -P "$( dirname "$0" )" && pwd )"
base_dir=/media/Data2/mouna_testing_data
data_dir=$base_dir/200515_NS500340_0486_AH5LLLAFX2/Alignment_1/20200516_034826/Fastq
fastq_1=$data_dir/FR1_S16_L001_R1_001.fastq.gz
fastq_2=$data_dir/FR1_S16_L001_R2_001.fastq.gz #If using Illumina paired data
fastq_3=$data_dir/FR1_S16_L002_R1_001.fastq.gz
fastq_4=$data_dir/FR1_S16_L002_R2_001.fastq.gz #If using Illumina paired data
fastq_5=$data_dir/FR1_S16_L003_R1_001.fastq.gz
fastq_6=$data_dir/FR1_S16_L003_R2_001.fastq.gz #If using Illumina paired data
fastq_7=$data_dir/FR1_S16_L004_R1_001.fastq.gz
fastq_8=$data_dir/FR1_S16_L004_R2_001.fastq.gz #If using Illumina paired data
INTERVAL=$base_dir/targets_GeneList_Mouna.bed

# Update with the location of the reference data files
fasta=/media/Data4/ref/ucsc_hg38/hg38.fa
#fasta=/home/massi/Documents/ucsc_hg38/hg38.fa
#dbsnp=$data_dir/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
dbsnp=/media/Data4/ref/GATK_ref/hg38/v0/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
known_1000G_indels=/media/Data4/ref/GATK_ref/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
known_Mills_indels=/media/Data4/ref/GATK_ref/hg38/v0/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
targets=/media/Data2/mouna_testing_data/targets_GeneList_Mouna.bed

# Set SENTIEON_LICENSE if it is not set in the environment
#export SENTIEON_LICENSE=/home/jyoda68/Desktop/PACK1/Sentieon/Linkoping_University_eval.lic

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/opt/sw/bioinfo-tools/sources/sentieon-genomics-202010

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$base_dir/tmp
mkdir -p $SENTIEON_TMPDIR

# It is important to assign meaningful names in actual cases.
# It is particularly important to assign different read group names.
sample="FR1"
group="L001"
platform="ILLUMINA"

# Other settings
nt=24 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$base_dir/results/gene_panel/FR1
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
# Controlling memory usage in BWA
export bwt_max_mem=16G

# speed up memory allocation malloc in bwa
export LD_PRELOAD=$SENTIEON_INSTALL_DIR/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1

#( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort $bam_option -r $fasta -o sorted.bam -t $nt --sam2bam -i -

#Run alignment for 1st input file set
(sentieon bwa mem -M -R '@RG\tID:L001 \tSM:FR1 \tPL:ILLUMINA' \
  -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error' ) \
  | sentieon util sort -o SORTED_BAM_1.bam -t $nt $bam_option -r $fasta --sam2bam -i -
#Run alignment for 2nd input file set
(sentieon bwa mem -M -R '@RG\tID:L002 \tSM:FR1 \tPL:ILLUMINA' \
  -t $nt -K 10000000 $fasta $fastq_3 $fastq_4 || echo -n 'error' ) \
  | sentieon util sort -o SORTED_BAM_2.bam -t $nt $bam_option -r $fasta --sam2bam -i -
#Run alignment for 1st input file set
(sentieon bwa mem -M -R '@RG\tID:L003 \tSM:FR1 \tPL:ILLUMINA' \
  -t $nt -K 10000000 $fasta $fastq_5 $fastq_6 || echo -n 'error' ) \
  | sentieon util sort -o SORTED_BAM_3.bam -t $nt $bam_option -r $fasta --sam2bam -i -
#Run alignment for 2nd input file set
(sentieon bwa mem -M -R '@RG\tID:L004 \tSM:FR1 \tPL:ILLUMINA' \
  -t $nt -K 10000000 $fasta $fastq_7 $fastq_8 || echo -n 'error' ) \
  | sentieon util sort -o SORTED_BAM_4.bam -t $nt $bam_option -r $fasta --sam2bam -i -


#Run dedup on both BAM files
sentieon driver -t $nt -i SORTED_BAM_1.bam -i SORTED_BAM_2.bam -i SORTED_BAM_3.bam -i SORTED_BAM_4.bam \
  --algo LocusCollector --fun score_info score.txt
sentieon driver -t $nt -i SORTED_BAM_1.bam -i SORTED_BAM_2.bam -i SORTED_BAM_3.bam -i SORTED_BAM_4.bam \
  --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt $bam_option deduped.bam







