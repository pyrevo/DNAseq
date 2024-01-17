#!/usr/bin/env python

l0 = []
l = []

dir = "/media/Data2/mouna_testing_data/results/gene_panel"

with open("files.txt", "r") as f:
    for line in f:
        id = line.split('_')[0]
        if id not in l0:
            l0.append(id)
            l.append('-v {}/{}/{}{}'.format(dir, id, id, "_GVCF.vcf.gz"))
        else:
            continue

#print(l)
com = "sentieon driver -r /media/Data4/ref/ucsc_hg38/hg38.fa --algo GVCFtyper"
lm = " ".join(l)
out = "gene_panel.vcf.gz"

command = '{} {} {}/{}'.format(com, lm, dir, out)
print(command)
