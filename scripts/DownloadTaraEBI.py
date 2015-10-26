"""
This is a script to download the submitted fastq files from EBI, for the sequences of the TARA Oceans Study
"""


def verify_link():
    pass


def get_xml(ena_url):
    import urllib2
    import xml.etree.ElementTree as ET

    input_file = urllib2.urlopen(ena_url + "&display=xml")
    xml_content = input_file.read()
    input_file.close()

    root = ET.fromstring(xml_content)

    for entry in root.iter("FILE"):
        checksum = entry.attrib["checksum"] # MD5 checksum to verify
        fastq_name = entry.attrib["filename"] # Name of the file

        print checksum, fastq_name



#http://www.ebi.ac.uk/ena/data/view/ERR598950,ERR599095&display=xml