#!/usr/bin/env python

from collections import defaultdict
d = defaultdict(list)

with open("files.txt", "r") as f:
    for line in f:
        id = line.split('_')[0]
        d[id].append(line.strip())

#print(d)

for id in d:
    l = " ".join(d[id])
    print('{} {} {}'.format("./dnaseq_gvcf.sh", l, '\nwait'))

#python scripter.py > runscript.sh