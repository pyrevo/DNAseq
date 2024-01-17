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
group_1="L001"
group_2="L002"
group_3="L003"
group_4="L004"
platform="ILLUMINA"

# Other settings
nt=24 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$base_dir/results/gene_panel/$sample
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

#Run alignment for both input file sets sequentially
(sentieon bwa mem -M -R "@RG\tID:$group_1\tSM:$sample\tPL:$platform" \
  -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 && sentieon bwa mem -M -R \
  "@RG\tID:$group_2\tSM:$sample\tPL:$platform" \
  -t $nt -K 10000000 $fasta $fastq_3 $fastq_4 && sentieon bwa mem -M -R \
  "@RG\tID:$group_3\tSM:$sample\tPL:$platform" \
  -t $nt -K 10000000 $fasta $fastq_5 $fastq_6 && sentieon bwa mem -M -R \
  "@RG\tID:$group_4\tSM:$sample\tPL:$platform" \
  -t $nt -K 10000000 $fasta $fastq_7 $fastq_8 || echo -n 'error' ) \
  | sentieon util sort $bam_option -r $fasta -o sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2. Metrics
# ******************************************
sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads
# To mark duplicate reads only without removing them, remove "--rmdup" in the second command
# ******************************************
sentieon driver -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
sentieon driver -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt $bam_option deduped.bam 

# ******************************************
# 4. Indel realigner
# This step is optional for haplotyper-based
# caller like HC, but necessary for any
# pile-up based caller. If you want to use
# this step, you need to update the rest of
# the commands to use realigned.bam instead
# of deduped.bam
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels $bam_option realigned.bam

# ******************************************
# 5. Base recalibration
# ******************************************

# Perform recalibration
sentieon driver -r $fasta -t $nt -i deduped.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table

# Perform post-calibration check (optional)
sentieon driver --interval $INTERVAL -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
sentieon driver --interval $INTERVAL -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv   
sentieon plot QualCal -o recal_plots.pdf recal.csv

# ******************************************
# 5b. ReadWriter to output recalibrated bam
# This stage is optional as variant callers
# can perform the recalibration on the fly
# using the before recalibration bam plus
# the recalibration table
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo ReadWriter recaled.bam


# ******************************************
# 6. HC Variant caller
# Note: Sentieon default setting matches versions before GATK 3.7.
# Starting GATK v3.7, the default settings have been updated multiple times. 
# Below shows commands to match GATK v3.7 - 4.1
# Please change according to your desired behavior.
# ******************************************

# Matching GATK 3.7, 3.8, 4.0
#$SENTIEON_INSTALL_DIR/bin/sentieon driver --interval $INTERVAL -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=10 output-hc.vcf.gz
sentieon driver --interval $INTERVAL -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_mode gvcf ${sample}_GVCF.vcf.gz
