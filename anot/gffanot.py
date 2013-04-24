#!/usr/bin/python
import sys
import re

gff_file = sys.argv[1]
annot_file = sys.argv[2]

gff = open(gff_file,'r').readlines()
ann = open(annot_file,'r').readlines()

for gff_line in gff:
    seq_name = gff_line.rsplit()[8][3:]
    product=""
    pat = re.compile("^.*\\b" + seq_name + "\\b.*$")
    for ann_line in ann:
        if re.findall(pat,ann_line):
            product = ann_line.split("\t")[2][:-1]
            go = ann_line.split("\t")[1]
            break
    print gff_line[:-1] + ";product="+product + ";" + go
