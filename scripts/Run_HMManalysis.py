"""
This scripts run the HMM analysis for the Tara oceans samples.
The approach is to download each sample file, then use hmmsearch against the
ResFam database, store the results in a new file, and delete the downloaded sequence.

The output is the result for the hmmsearch, with the sample label as the file name (with prefix hmmsearch)
"""

from collections import defaultdict
import argparse
import os
import subprocess
from Bio import SeqIO


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


def count_fastq_reads(input_file):
    count = 0
    for record in SeqIO.parse(open(input_file, 'r'), "fastq"):
        count += 1

    return count

# ---------------- #
program_description = "Script that runs the HMM analysis"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-i", "--input_table", type=str,
                    help="Input tabular file with the Tara oceans sample", required=True)
parser.add_argument("-o", "--output_folder", type=str, help="Output folder", required=True)
parser.add_argument("-d", "--hmm_files", type=str, help="HMM files to use", required=True)
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

for sample in tara_data:
    print "\n# Processing sample: %s" % sample
    sample_folder = args.output_folder + "/" + sample

    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)

    for entry in read_xml(tara_data[sample]):
        seq_type, checksum, seq_file = entry

        output_seq_file = sample_folder + "/" + seq_file.split("/")[-1]

        file_url = "ftp://ftp.sra.ebi.ac.uk/vol1/" + seq_file

        print "## Downloading file %s \n" % file_url

        # getfile = urllib.URLopener()
        # getfile.retrieve(file_url, output_seq_file)

        subprocess.call(["wget", file_url, "-O", output_seq_file])

        # Unzipping file if fastq
        if seq_type == "fastq":

            print "### Uncompressing file\n"
            subprocess.call(["gunzip", output_seq_file])

            base = os.path.basename(output_seq_file)
            file_prefix = os.path.splitext(base)[0]

            fastq_filename = sample_folder + "/" + file_prefix
            output_faa = sample_folder + "/" + file_prefix + ".faa"
            output_hmm = sample_folder + "/" + file_prefix + ".hmmsearch"
            logfile_hmm = sample_folder + "/" + file_prefix + ".logfile"

            # Get number of reads on the file
            print "#### Counting records\n"
            #read_count = count_fastq_reads(fastq_filename)


            proc = subprocess.Popen(["wc", "-l", fastq_filename], stdout=subprocess.PIPE)
            output = proc.stdout.read()

            read_count = int(output.split(" ")[0])/4
            summary_table.write(sample + "\t" + file_prefix + "\t" + str(read_count) + "\t")

            # Translate sequences

            print "#### Running Transeq\n"
            subprocess.call(["transeq", "-sequence", fastq_filename, "-outseq",
                             output_faa, "-frame", "6", "-clean"])

            # Count proteins
            print "### Counting proteins\n"
            peptide_count = subprocess.check_output('grep -c ">" ' + output_faa)
            peptide_count = peptide_count.rstrip()
            summary_table.write(str(peptide_count) + "\n")

            # Run hmmsearch
            print "##### Running HMM searches\n"
            subprocess.call(["hmmsearch", "--cpu", str(args.cpus), "--cut_ga", "--tblout", output_hmm, "-o", logfile_hmm,
                             args.hmm_files, output_faa])

            # Delete the files
            os.remove(fastq_filename)
            os.remove(output_faa)

        elif seq_type == "sff":
            base = os.path.basename(output_seq_file)
            file_prefix = os.path.splitext(base)[0].split(".")[0]

            sff_fastq_name = sample_folder + "/" + file_prefix + ".fastq"
            output_faa = sample_folder + "/" + file_prefix + ".faa"
            output_hmm = sample_folder + "/" + file_prefix + ".hmmsearch"
            logfile_hmm = sample_folder + "/" + file_prefix + ".logfile"

            # Get number of read on the file
            print "### Converting sff file\n"
            read_count = SeqIO.convert(output_seq_file, "sff", sff_fastq_name, "fastq")
            summary_table.write(sample + "\t" + file_prefix + "\t" + str(read_count) + "\n")

            # Translate the sequences
            print "#### Running Transeq\n"
            subprocess.call(["transeq", "-sequence", sff_fastq_name, "-outseq",
                             output_faa, "-frame", "6", "-clean"])

            # Count proteins
            print "### Counting proteins\n"
            peptide_count = subprocess.check_output('grep -c ">" ' + output_faa)
            peptide_count = peptide_count.rstrip()
            summary_table.write(str(peptide_count) + "\n")

            # Run hmmsearch
            print "##### Running HMM searches\n"
            subprocess.call(["hmmsearch", "--cpu", str(args.cpus), "--cut_ga", "--tblout", output_hmm, "-o", logfile_hmm,
                             args.hmm_files, output_faa])

            # Delete the files
            os.remove(sff_fastq_name)
            os.remove(output_faa)


summary_table.close()
