from MA import *
from random import choice
import tempfile
import random

# create a reference string
reference = ""
for _ in range(1000):
    reference += choice(['C', 'T', 'G', 'A'])
# create pack and fmd index from that string
reference_pack = Pack()
reference_pack.append(
    "random pack", "sequence description", NucSeq(reference))


for _ in range(100):
    query = ""
    qual = ""
    for _ in range(100):
        query += choice(['C', 'T', 'G', 'A'])
        qual += str(chr(random.randint(33, 126))) # random quality

    # write a fasta file
    fasta_file_name = tempfile.gettempdir() + "/.tmp.fasta"
    sam_file_name = tempfile.gettempdir() + "/.tmp.sam"
    with open(fasta_file_name, "w") as out_file:
        out_file.write("@name desc\n")
        out_file.write(query)
        out_file.write("\n+\n")
        out_file.write(qual)
        out_file.write("\n")

    file_reader = FileReader(ParameterSetManager(), libMA.path(fasta_file_name))
    query_nuc_seq = file_reader.execute()
    alignment_vec = libMA.AlignmentVector()
    file_writer = FileWriter(ParameterSetManager(), sam_file_name, reference_pack)
    file_writer.execute(query_nuc_seq, alignment_vec, reference_pack)
    del file_writer

    with open(sam_file_name, "r") as in_file:
        for line in in_file:
            if line[0] != '@':
                #print("|", line, "|", sep="")
                fields = line[:-1].split("\t")
                #print("|" + fields[10] + "|", "|" + qual + "|", sep='\n')
                assert fields[9] == query
                assert fields[10] == qual

for _ in range(100):
    query = ""
    for _ in range(100):
        query += choice(['C', 'T', 'G', 'A'])

    # write a fasta file
    fasta_file_name = tempfile.gettempdir() + "/.tmp.fasta"
    sam_file_name = tempfile.gettempdir() + "/.tmp.sam"
    with open(fasta_file_name, "w") as out_file:
        out_file.write(">name desc\n")
        out_file.write(query)
        out_file.write("\n")

    file_reader = FileReader(ParameterSetManager(), libMA.path(fasta_file_name))
    query_nuc_seq = file_reader.execute()
    alignment_vec = libMA.AlignmentVector()
    file_writer = FileWriter(ParameterSetManager(), sam_file_name, reference_pack)
    file_writer.execute(query_nuc_seq, alignment_vec, reference_pack)
    del file_writer

    with open(sam_file_name, "r") as in_file:
        for line in in_file:
            if line[0] != '@':
                #print("|", line, "|", sep="")
                fields = line[:-1].split("\t")
                #print("|" + fields[10] + "|", "|" + qual + "|", sep='\n')
                assert fields[9] == query
                assert fields[10] == "*"