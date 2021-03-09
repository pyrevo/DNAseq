################################################################################
## DNAseq Pipeline
## Somatic variant calling with Sentieon DNAseq (from fastq files to vcf)
##
## Authors: Massimiliano Volpe, Jyotirmoy Das
## Email: massimiliano.volpe@liu.se, jyotirmoy.das@liu.se
## Date: 11/02/2021
## Developed on behalf of the Bioinformatics Core Facility, Linköping University
##
## Rules:
## - Raw data should be stored in the same folder
## - Paths to raw data and reference must be set in the config.json file
## - Fastq filename suffix must be set in the config.json file, e.g.:
##      /path_to_fastq/FR1_S16_L001_R1_001.fastq.gz --> "_L001_R1_001.fastq.gz"
##      /path_to_fastq/FR1_S16_L001_R2_001.fastq.gz --> "_L001_R2_001.fastq.gz"
################################################################################


# Functions -------------------------------------------------------------------
def id_maker(path, sep, sample, suffix, d):
    f = "".join([path, sep, sample, suffix])
    # run:flowcell:lane
    # l = subprocess.check_output("zcat " + f + " | head -n 1 | cut -d ':' -f 2,3,4 | sed  s/:/./g | sed 's/@//'", shell=True).strip().decode()
    # flowcell:lane
    l = subprocess.check_output("zcat " + f + " | head -n 1 | cut -d ':' -f 3,4 | sed  s/:/./g | sed 's/@//'", shell=True).strip().decode()
    #d[sample] = '_'.join([sample, group, l])
    d[sample] = l
    return(d)


# Globals ---------------------------------------------------------------------
import subprocess
from collections import defaultdict

configfile:
    "config.json"

#workdir:
#    config['workdir']


R1SUFFIX = config['R1_suffix']
R2SUFFIX = config['R2_suffix']
R3SUFFIX = config['R3_suffix']
R4SUFFIX = config['R4_suffix']
R5SUFFIX = config['R5_suffix']
R6SUFFIX = config['R6_suffix']
R7SUFFIX = config['R7_suffix']
R8SUFFIX = config['R8_suffix']

DATA = config['dataset']
SAMPLES, = glob_wildcards(DATA + "/{sample}" + R1SUFFIX)
RESULTS = config['workdir'] + '/{sample}/'
BAMS = RESULTS + 'bams/'
LOGS = RESULTS + 'logs/'
METRICS = RESULTS + 'metrics/'
PLOTS = RESULTS + 'plots/'
MARKDUP = RESULTS + 'markdup/'
RECAL = RESULTS + 'baserecal/'

fasta = config['reference']
dbsnp = config['dbsnp']
known_Mills_indels = config['known_Mills_indels']
known_1000G_indels = config['known_1000G_indels']


#print(SAMPLES)

#d1 = defaultdict(list)

#for sample in SAMPLES:
#    for file in os.listdir(config['dataset']):
#        if file.startswith(sample):
#            d1[sample].append(file)

#d2 = {x:sorted(d1[x]) for x in d1.keys()}
#print(d2)


# Rules -----------------------------------------------------------------------
rule all:
    input:
        bam = expand(BAMS + "{sample}.deduped.bam", sample=SAMPLES),
        rdt = expand(RECAL + '{sample}.recal_data.table', sample=SAMPLES),
        vcf = expand(RESULTS + '{sample}.dnaseq.vcf.gz', sample=SAMPLES)
        #bam = expand(BAMS + "{sample}.normal_sorted.bam", sample=SAMPLES),
        #normal_ml = expand(LOGS + '{sample}.normal_metrics.log', sample=SAMPLES),


rule bwa:
    input:
        #fq = lambda wildcards: d2[wildcards.sample],
        R1 = DATA + "/{sample}" + R1SUFFIX,
        R2 = DATA + "/{sample}" + R2SUFFIX,
        R3 = DATA + "/{sample}" + R3SUFFIX,
        R4 = DATA + "/{sample}" + R4SUFFIX,
        R5 = DATA + "/{sample}" + R5SUFFIX,
        R6 = DATA + "/{sample}" + R6SUFFIX,
        R7 = DATA + "/{sample}" + R7SUFFIX,
        R8 = DATA + "/{sample}" + R8SUFFIX,
        ref = config['reference']
    output:
        sam1 = temp(BAMS + '{sample}_1.sam'),
        bam1 = temp(BAMS + '{sample}_1.sorted.bam'),
        sam2 = temp(BAMS + '{sample}_2.sam'),
        bam2 = temp(BAMS + '{sample}_2.sorted.bam'),
        sam3 = temp(BAMS + '{sample}_3.sam'),
        bam3 = temp(BAMS + '{sample}_3.sorted.bam'),
        sam4 = temp(BAMS + '{sample}_4.sam'),
        bam4 = temp(BAMS + '{sample}_4.sorted.bam')
    log:
        bwa = LOGS + '{sample}.bwa.log',
        sort = LOGS + '{sample}.sort.log'
    params:
        K = 10000000,
        ID = ["L001", "L002", "L003", "L004"],
        SM = "{sample}",
        PL = config["platform"]
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon bwa mem -M -R '@RG\\tID:{params.ID[0]}\\tSM:{params.SM}\\tPL:{params.PL}' -t {threads} -K {params.K} -o {output.sam1} {input.ref} {input.R1} {input.R2} >> {log.bwa} 2>&1
        sentieon util sort -r {input.ref} -i {output.sam1} -o {output.bam1} -t {threads} --sam2bam >> {log.sort} 2>&1

        sentieon bwa mem -M -R '@RG\\tID:{params.ID[1]}\\tSM:{params.SM}\\tPL:{params.PL}' -t {threads} -K {params.K} -o {output.sam2} {input.ref} {input.R3} {input.R4} >> {log.bwa} 2>&1
        sentieon util sort -r {input.ref} -i {output.sam2} -o {output.bam2} -t {threads} --sam2bam >> {log.sort} 2>&1

        sentieon bwa mem -M -R '@RG\\tID:{params.ID[2]}\\tSM:{params.SM}\\tPL:{params.PL}' -t {threads} -K {params.K} -o {output.sam3} {input.ref} {input.R5} {input.R6} >> {log.bwa} 2>&1
        sentieon util sort -r {input.ref} -i {output.sam3} -o {output.bam3} -t {threads} --sam2bam >> {log.sort} 2>&1

        sentieon bwa mem -M -R '@RG\\tID:{params.ID[3]}\\tSM:{params.SM}\\tPL:{params.PL}' -t {threads} -K {params.K} -o {output.sam4} {input.ref} {input.R7} {input.R8} >> {log.bwa} 2>&1
        sentieon util sort -r {input.ref} -i {output.sam4} -o {output.bam4} -t {threads} --sam2bam >> {log.sort} 2>&1
        """


rule markdup:
    input:
        bam1 = rules.bwa.output.bam1,
        bam2 = rules.bwa.output.bam2,
        bam3 = rules.bwa.output.bam3,
        bam4 = rules.bwa.output.bam4,
        ref = config['reference']
    output:
        ns = MARKDUP + '{sample}.score.txt',
        dm = MARKDUP + '{sample}.dedup_metrics.txt',
        bam = BAMS + '{sample}.deduped.bam',
        # cm = MARKDUP + '{sample}.coverage_metrics'
    log:
        LOGS + '{sample}.dedup.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        # sentieon driver -r {input.ref} -t {threads} -i {output.bam} --algo CoverageMetrics {output.cm} >> {log} 2>&1
        """
        sentieon driver -t {threads} -i {input.bam1} -i {input.bam2} -i {input.bam3} -i {input.bam4} --algo LocusCollector --fun score_info {output.ns} >> {log} 2>&1
        sentieon driver -t {threads} -i {input.bam1} -i {input.bam2} -i {input.bam3} -i {input.bam4} --algo Dedup --rmdup --score_info {output.ns} --metrics {output.dm} {output.bam} >> {log} 2>&1
        """


rule baserecal:
    input:
        bam = rules.markdup.output.bam,
        ref = config['reference'],
        bed = config['interval']
    output:
        rdt = RECAL + '{sample}.recal_data.table',
        post = RECAL + '{sample}.recal_data.table.post',
        recal = RECAL + '{sample}.recal.csv',
        rp = PLOTS + '{sample}.recal_plots.pdf',

    log:
        LOGS + '{sample}.recal.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver -r {input.ref} -t {threads} -i {input.bam} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.rdt} >> {log} 2>&1
        sentieon driver --interval {input.bed} -r {input.ref} -t {threads} -i {input.bam} -q {output.rdt} --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {output.post} >> {log} 2>&1
        sentieon driver --interval {input.bed} -t {threads} --algo QualCal --plot --before {output.rdt} --after {output.post} {output.recal} >> {log} 2>&1
        sentieon plot QualCal -o {output.rp} {output.recal}
        """


rule variant_calling:
    input:
        bam = rules.markdup.output.bam,
        ref = config['reference'],
        bed = config['interval'],
        rdt = rules.baserecal.output.rdt
    output:
        vcf = RESULTS + '{sample}.dnaseq.vcf.gz',
    log:
        LOGS + '{sample}.dnaseq.log'
    params:
        emit = 10,
        call = 10
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        sentieon driver --interval {input.bed} -r {input.ref} -t {threads} -i {input.bam} -q {input.rdt} --algo Haplotyper --dbsnp {dbsnp} --emit_conf={params.emit} --call_conf={params.call} {output.vcf} >> {log} 2>&1
        """


# markdup cm output is commented because output files are too big