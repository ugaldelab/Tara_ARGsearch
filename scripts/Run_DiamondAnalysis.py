"""
This scripts run the read mapping analysis for the Tara oceans samples.
The approach is to download each sample file, then use Diamond against the CARD database.
Also it will compare against a set of recA genes from the FuncGene database.

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


def parse_diamond_results(m8_table):
    """
    Used to parse the diamond results
    :param m8_table:
    :return:
    """

    results = defaultdict(int)

    for line in open(m8_table, 'r'):
        if not line.strip():
            continue

        line = line.rstrip()

        elements = line.split("\t")
        read = elements[0]
        hit = elements[1]

        # Just taking the top hit
        if read in results:
            continue

        else:
            if hit[-4:] == "recA":
                results["recA"] += 1
            else:
                results[hit] += 1

    return results


def get_gene_lengths(input_fasta):
    """
    Get the gene lengths for each gene in the fasta file. This is used to get the
    gene lenghts for the antibiotic database

    :param input_fasta:
    :return:
    """
    lengths = defaultdict(int)

    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        lengths[seq_record.id] = len(seq_record)

    return lengths


# ---------------- #
program_description = "Script that runs the read mapping analysis, using BWA"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-i", "--input_table", type=str,
                    help="Input tabular file with the Tara oceans sample", required=True)
parser.add_argument("-o", "--output_folder", type=str, help="Output folder", required=True)
parser.add_argument("-d", "--database", type=str, help="Diamond database", required=True)
parser.add_argument("-c", "--cpus", type=int, help="Number of CPUs to use", required=True)
#parser.add_argument("-g", "--arg_fasta", type=str, help="Fasta file with the ARGs used for search", required=True)
parser.add_argument("-a", "--arg_categories", type=str, help="File with the ARG categories", required=True)

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

# Prepare the input dataset of antibiotic genes
#length_arg = get_gene_lengths(args.arg_fasta)

# Read the file with the ARGs categories
arg_categories = defaultdict(list)

for line in open(args.arg_categories, 'r'):
    if not line.strip():
        continue

    line = line.rstrip()
    gene, function, aro = line.split("\t")
    gene = gene.rstrip(".p01")
    gene = gene.rstrip(".p1")
    arg_categories[gene].append((function, aro))

# Process each file

summary_table = open(args.output_folder + "/summary_table.txt", 'w')
file_sample_results = open(args.output_folder + "/sample_results.txt", 'w')

# Run Diamond, and parse results

for sample in tara_data:

    sample_results = defaultdict(int)

    print "\n# Processing sample: %s" % sample

    # Create output folder
    sample_folder = args.output_folder + "/" + sample

    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)

    # Get the xml file for each entry, get the sequences and run diamond
    for entry in read_xml(tara_data[sample]):
        seq_type, checksum, seq_file = entry

        output_seq_file = sample_folder + "/" + seq_file.split("/")[-1]

        file_url = "ftp://ftp.sra.ebi.ac.uk/vol1/" + seq_file

        # Download seqs file
        print "## Downloading file %s \n" % file_url

        subprocess.call(["wget", file_url, "-O", output_seq_file])
        # print seq_type

        # Files could be either Fastq or SFF
        if seq_type == "fastq":
            base = os.path.basename(output_seq_file)
            file_prefix = os.path.splitext(base)[0]

            fastq_file = sample_folder + "/" + os.path.basename(output_seq_file)
            output_diamond = sample_folder + "/" + file_prefix + ".daa"
            output_tab = sample_folder + "/" + file_prefix + ".m8"

            # Get the number of reads on the file
            print "### Counting records\n"

            proc = subprocess.Popen(["wc", "-l", fastq_file], stdout=subprocess.PIPE)
            line_count = proc.stdout.read()
            read_count = int(line_count.split()[0])/4

            summary_table.write(sample + "\t" + file_prefix + "\t" + str(read_count) + "\t")

            # Run diamond

            print "###### Running diamond\n"

            subprocess.call(["diamond", "blastx", "-p", "32", "-q", fastq_file, "-e", "0.00001", "--sensitive",
                             "-a", output_diamond, "-d", args.database])

            subprocess.call(["diamond", "view", "-a", output_diamond, "-o", output_tab])

            # Delete temporal files
            os.remove(fastq_file)
            os.remove(output_diamond)

            for gene, count in parse_diamond_results(output_tab).items():
                sample_results[gene] += count

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

            # Run Diamond
            print "###### Running diamond\n"
            subprocess.call(["diamond", "blastx", "-p", "32", "-q", sff_fastq_name, "-e", "0.00001", "--sensitive",
                             "-a", output_diamond, "-d", args.database])

            subprocess.call(["diamond", "view", "-a", output_diamond, "-o", output_tab])

        # Delete the files
            os.remove(output_seq_file)
            os.remove(sff_fastq_name)
            os.remove(output_diamond)

            # Parse the diamond results
            for gene, count in parse_diamond_results(output_tab).items():
                sample_results[gene] += count

    # Process sample results
    raw_sample_results = defaultdict(int)
    normalized_sample_results = defaultdict(float)

    recA_count = sample_results["recA"]

    for entry in sample_results:
        arg_functions = [value[0] for value in arg_categories[entry]]

        for function in arg_functions:
            raw_sample_results[function] += int(sample_results[entry])

    for entry in raw_sample_results:
        norm_count = int(raw_sample_results[entry]) / float(recA_count)
        normalized_sample_results[entry] = norm_count

    file_sample_results.write("#" + sample + "\n")
    for result in normalized_sample_results:
        file_sample_results.write(result + "\t" + str(normalized_sample_results[result]) + "\n")

    file_sample_results.flush()
    summary_table.flush()

    os.fsync(file_sample_results.fileno())
    os.fsync(summary_table.fileno())

summary_table.close()
file_sample_results.close()