import sys
from Bio import SeqIO

#Corta seq nlen bp de largo upstream (500 por default) a partir de un GFF
#Params: lista.genes, file.gff, [len.margen]

lista = open(sys.argv[1],'r').readlines()
index = SeqIO.index(sys.argv[2], 'fasta')
try:
    nlen = int(sys.argv[3])
except:
    nlen = 500

for i in lista:
    chr = i.split()[0]
    name = i.split()[8].split(";")[0].replace("ID=","")
    strand = i.split()[6]
    start = int(i.split()[3])
    end = int(i.split()[4])

    if strand == "+":
        seq = index[chr].seq[start-nlen:end]
    else:
        seqaux = index[chr].seq[start:end+nlen]
        seq = seqaux.reverse_complement()

    if len(seq) < nlen:
        seq = index[chr].seq[:end] if strand == "+" else index[chr].seq[start:]

    print ">" + name
    n=60
    for i in xrange(0, len(seq), n):
        print seq[i:i+n]
