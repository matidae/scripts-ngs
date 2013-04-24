import sys, re, urllib
from os import popen
from xml.etree import cElementTree

name_file = open(sys.argv[1],'r').readlines()
blast_file = open(sys.argv[2],'r').readlines()
go_file = open(sys.argv[3],'r').readlines()
ipro_file = open(sys.argv[4],'r').readlines()
ec_file = open(sys.argv[5],'r').readlines()
nr_db = "/mnt2/1T/nr/nr"

def getBlast(gene_name):
    pat = re.compile("^.*\\b" + gene_name + "\\b.*$")
    blast_result_mod = blast_result = blast_desc = aln = idn = eva = ""
    n = 0
    for blast_line in blast_file:        
        if re.findall(pat, blast_line):
            n += 1
            gi = blast_line.split("|")[1]
            aln = blast_line.split("\t")[3]
            idn = blast_line.split("\t")[2]
            eva = blast_line.split("\t")[10]
            blast_desc = popen("fastacmd -s " + gi + " -d " + nr_db + "| head -n1" ).read().rstrip()
            blast_result += "blast=aln:" + aln + "|idn:" + idn + "|eva:" +eva + blast_desc.rstrip()  + ";"
            if n>3: break
    blast_result_mod = blast_result.replace(" ","_").replace("_>","|").replace(">","|")
    return blast_result_mod

def getGo(gene_name):
    pat = re.compile("^.*\\b" + gene_name + "\\b.*$")
    go_desc = goc = gop = gof = ""
    for go_line in go_file:
        go_line_mod = go_line.replace(" ","_")
        if re.findall(pat, go_line_mod):
            godb = go_line_mod.split("\t")[4].rstrip()
            if godb == "C":
                goCterm=["cell" ,"cell_part"]
                if go_line_mod.split("\t")[3] not in goCterm:
                    goc += "goC=" + go_line_mod.split("\t")[2] +  "|" + go_line_mod.split("\t")[3] + ";"
            if godb == "P":
                goc += "goP=" + go_line_mod.split("\t")[2] +  "|" + go_line_mod.split("\t")[3] + ";"
            if godb == "F":
                goc += "goF=" + go_line_mod.split("\t")[2] +  "|" + go_line_mod.split("\t")[3] + ";"
    go_desc = goc + gop + gof
    return go_desc

def getIpro(gene_name):
    pat = re.compile("^.*\\b" + gene_name + "\\b.*$")
    ipro_desc= ""
    for ipro_line in ipro_file:
        ipro_line_mod = ipro_line.replace(" ","_")
        if re.findall(pat, ipro_line_mod):
            if len(ipro_line_mod.split("\t"))==6:
                ipro_desc += "ipro=" + ipro_line_mod.split("\t")[2] + "|" + ipro_line_mod.split("\t")[3] + "|" + ipro_line_mod.split("\t")[4] + "|" +  ipro_line_mod.split("\t")[5].rstrip() + ";"
    return ipro_desc

def getEc(gene_name):
    pat = re.compile("^.*\\b" + gene_name + "\\b.*$")
    ec_desc = ""
    for ec_line in ec_file:
        ec_line_mod = ec_line.replace(" ","_")
        if re.findall(pat, ec_line):
            ec_desc += "ec=" + ec_line_mod.split("\t")[1].rstrip() + ";"
    return ec_desc

def buildGffLine(gene_name, blast_desc, go_desc, ipro_desc, ec_desc):
    gff_line="ID=" + gene_name + ";"
    if blast_desc: gff_line +=  blast_desc 
    if go_desc: gff_line +=  go_desc 
    if ipro_desc: gff_line += ipro_desc 
    if ec_desc: gff_line += ec_desc 
    return gff_line

if __name__ == "__main__":
    for gene in name_file:
        gene_name = gene.rstrip()
        blast_desc = getBlast(gene_name)
        go_desc = getGo(gene_name)
        ipro_desc = getIpro(gene_name)
        ec_desc = getEc(gene_name)
    	print buildGffLine(gene_name, blast_desc, go_desc, ipro_desc, ec_desc)
	
