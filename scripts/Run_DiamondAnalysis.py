"""
This scripts run the read mapping analysis for the Tara oceans samples.
The approach is to download each sample file, then use Bowtie2 to map the reads to the
Resfinder database, store the results in a new file, and delete the downloaded sequence.

The output is a bam file for the mapped reads
"""

from collections import defaultdict
import argparse
import os
import subprocess
from Bio import SeqIO


####USE DIAMOND###


def read_xml(ena_url):
    """
    Recieves a plain URL for a nucleotide present in ENA. For example,
    http://www.ebi.ac.uk/ena/data/view/ERR562552
    It will query the ENA database, obtain the result in xml, and parse
    the sequence file, the type and the checksum value

    :param ena_url:
    :return:
    """
    import urllib2
    import xml.etree.ElementTree as ET

    input_file = urllib2.urlopen(ena_url + "&display=xml")
    xml_content = input_file.read()
    input_file.close()

    root = ET.fromstring(xml_content)

    xml_files = list()

    for entry in root.iter("FILE"):
        checksum = entry.attrib["checksum"]  # MD5 checksum to verify
        fastq_name = entry.attrib["filename"]  # Name of the file
        file_type = entry.attrib["filetype"]

        xml_files.append((file_type, checksum, fastq_name))

    return xml_files

def parse_diamond_results():
    pass

# ---------------- #
program_description = "Script that runs the read mapping analysis, using BWA"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-i", "--input_table", type=str,
                    help="Input tabular file with the Tara oceans sample", required=True)
parser.add_argument("-o", "--output_folder", type=str, help="Output folder", required=True)
parser.add_argument("-d", "--database", type=str, help="Diamond database", required=True)
parser.add_argument("-c", "--cpus", type=int, help="Number of CPUs to use", required=True)

args = parser.parse_args()

# Create the output directory
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

# Read the input table
tara_data = defaultdict()
entry_count = 0

for line in open(args.input_table, 'r'):
    if not line.strip():
        continue

    if entry_count == 0:
        entry_count += 1
        continue

    line = line.rstrip()
    sample_label = line.split("\t")[4]  # I'm using the PANGAEA sample ID, which seems to be unique
    ena_url = line.split("\t")[3]

    tara_data[sample_label] = ena_url
    entry_count += 1

print "Total of %d samples" % (entry_count-1)

# Process each file

summary_table = open(args.output_folder + "/summary_table.txt", 'w')

#Run Diamond, and parse results

for sample in tara_data:
    print "\n# Processing sample: %s" % sample

    #Create output folder
    sample_folder = args.output_folder + "/" + sample

    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)

    #Get the xml file for each entry, get the sequences and run diamond
    for entry in read_xml(tara_data[sample]):
        seq_type, checksum, seq_file = entry

        output_seq_file = sample_folder + "/" + seq_file.split("/")[-1]

        file_url = "ftp://ftp.sra.ebi.ac.uk/vol1/" + seq_file

        #Download seqs file
        print "## Downloading file %s \n" % file_url

        subprocess.call(["wget", file_url, "-O", output_seq_file])
        #print seq_type

        # Files could be either Fastq or SFF
        if seq_type == "fastq":
            base = os.path.basename(output_seq_file)
            file_prefix = os.path.splitext(base)[0]

            fastq_file = sample_folder + "/" + os.path.basename(output_seq_file)
            output_diamond = sample_folder + "/" + file_prefix + ".daa"
            output_tab = sample_folder + "/" + file_prefix + ".m8"

            #Get the number of reads on the file
            print "### Counting records\n"

            proc = subprocess.Popen(["wc", "-l", fastq_file], stdout=subprocess.PIPE)
            line_count = proc.stdout.read()
            read_count = int(line_count.split()[0])/4

            summary_table.write(sample + "\t" + file_prefix + "\t" + str(read_count) + "\t")

            #Run diamond

            print "###### Running diamond\n"

            subprocess.call(["diamond", "blastx", "-p", "32", "-q", fastq_file, "-e", "0.00001", "--sensitive",
                             "-a", output_diamond, "-d", args.database])

            subprocess.call(["diamond", "view", "-a", output_diamond, "-o", output_tab])

            # Delete temporal files
            os.remove(fastq_file)
            os.remove(output_diamond)

        elif seq_type == "sff":
            base = os.path.basename(output_seq_file)
            file_prefix = os.path.splitext(base)[0].split(".")[0]
            sff_fastq_name = sample_folder + "/" + file_prefix + ".fastq"
            output_diamond = sample_folder + "/" + file_prefix + ".daa"
            output_tab = sample_folder + "/" + file_prefix + ".m8"

            # Get number of read on the file
            print "### Converting sff file\n"
            read_count = SeqIO.convert(output_seq_file, "sff", sff_fastq_name, "fastq")
            summary_table.write(sample + "\t" + file_prefix + "\t" + str(read_count) + "\t")

        #     # Run Diamond
            print "###### Running diamond\n"
            subprocess.call(["diamond", "blastx", "-p", "32", "-q", sff_fastq_name, "-e", "0.00001", "--sensitive",
                             "-a", output_diamond, "-d", args.database])

            subprocess.call(["diamond", "view", "-a", output_diamond, "-o", output_tab])

        # Delete the files
            os.remove(output_seq_file)
            os.remove(sff_fastq_name)
            os.remove(output_diamond)


summary_table.close()


#Old bowtie2 code
    # read_check = 0
    # for entry in read_xml(tara_data[sample]):
    #     seq_type, checksum, seq_file = entry
    #
    #     if seq_type == "sff":
    #         print seq_file
    #
    #     elif seq_type == "fastq":
    #         if read_check == 1:
    #             continue
    #
    #         read_check += 1 # To avoid reading the same read pair two times
    #         read1, read2 = ["ftp://ftp.sra.ebi.ac.uk/vol1/" + fastq[2] for fastq in read_xml(tara_data[sample])]
    #
    #         output_read1_file = sample_folder + "/" + read1.split("/")[-1]
    #         output_read2_file = sample_folder + "/" + read2.split("/")[-1]
    #
    #         print " ## Downloading file %s \n" % output_read1_file
    #         subprocess.call(["wget", read1, "-O", output_read1_file])
    #
    #         print " ## Downloading file %s \n" % output_read2_file
    #         subprocess.call(["wget", read1, "-O", output_read2_file])
    #
    #         #Run Bowtie2
    #
    #         #Bajar test de ejemplo
    #         #Probar parametros bowtie2



