import urllib
import sys
from xml.etree import cElementTree

ec_number = sys.argv[1]

def get_ec_name(ec_number):
    ec = ec_number.split(".")
    if ec[3] == "0": ec[3] = "1"
    if ec[2] == "0": ec[2] = "1"
    if ec[1] == "0": ec[1] = "1"
    xml = urllib.urlopen("ftp://ftp.ebi.ac.uk/pub/databases/intenz/xml/ASCII/EC_"+ec[0]+"/EC_"+ec[0]+"."+ec[1]+"/EC_"+ec[0]+"."+ec[1]+"."+ec[2]+"/EC_"+ec[0]+"."+ec[1]+"."+ec[2]+"."+ec[3]+".xml")
    for event, element in cElementTree.iterparse(xml):
        if ec_number.split(".")[3] != "0":
            if element.tag == '{http://www.ebi.ac.uk/intenz}enzyme':
                for child in element.getchildren():
                    if child.tag == '{http://www.ebi.ac.uk/intenz}accepted_name':
                        if child.text != None:
                            return child.text
                    if child.tag == '{http://www.ebi.ac.uk/intenz}transferred':
                            for subchild in child:
                                if subchild.tag == '{http://www.ebi.ac.uk/intenz}note':
                                    return subchild.text.split(".")[0]

        else:
            if element.tag == '{http://www.ebi.ac.uk/intenz}ec_sub-subclass':
                for child in element.getchildren():
                    if child.tag == '{http://www.ebi.ac.uk/intenz}name':
                        return child.text

ec_file = open(sys.argv[1], 'r').readlines()
for line in ec_file:
    try:
        ec_desc = get_ec_name(line.split(":")[1].rstrip())
    except:
        ec_desc = "-"
    if ec_desc == None :ec_desc="-"
    print line.rstrip() + "|" + ec_desc
