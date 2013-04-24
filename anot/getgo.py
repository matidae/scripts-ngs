import urllib
import sys
from xml.etree import cElementTree as ElementTree

go_number = sys.argv[1]

def get_go_name(go_id):
    xml = urllib.urlopen("http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:"+go_id+"&format=oboxml")
    for event, element in ElementTree.iterparse(xml):
        if element.tag == 'term':
            for child in element.getchildren():
                if child.tag == 'name':
                    name = child.text 
                if child.tag == 'namespace':
                    namespace = child.text
                    return namespace + "\t" + name

try:
    print get_go_name(sys.argv[1])
except:
    print "-"
