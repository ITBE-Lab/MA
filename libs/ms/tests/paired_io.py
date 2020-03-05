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
    query1 = ""
    qual1 = ""
    for _ in range(100):
        query1 += choice(['C', 'T', 'G', 'A'])
        qual1 += str(chr(random.randint(33, 126))) # random quality
        
    query2 = ""
    qual2 = ""
    for _ in range(100):
        query2 += choice(['C', 'T', 'G', 'A'])
        qual2 += str(chr(random.randint(33, 126))) # random quality

    # write a fasta file
    fasta_1_file_name = tempfile.gettempdir() + "/.tmp1.fasta"
    fasta_2_file_name = tempfile.gettempdir() + "/.tmp2.fasta"
    sam_file_name = tempfile.gettempdir() + "/.tmp.p.sam"
    with open(fasta_1_file_name, "w") as out_file:
        out_file.write("@name desc\n")
        out_file.write(query1)
        out_file.write("\n+\n")
        out_file.write(qual1)
        out_file.write("\n")

    with open(fasta_2_file_name, "w") as out_file:
        out_file.write("@name desc\n")
        out_file.write(query2)
        out_file.write("\n+\n")
        out_file.write(qual2)
        out_file.write("\n")

    file_reader = PairedFileReader(ParameterSetManager(),
                                   libMA.filePathVector([libMA.path(fasta_1_file_name)]),
                                   libMA.filePathVector([libMA.path(fasta_2_file_name)]))
    query_vec = file_reader.execute()
    alignment_vec = libMA.AlignmentVector()
    file_writer = PairedFileWriter(ParameterSetManager(), sam_file_name, reference_pack)
    file_writer.execute(query_vec[0], query_vec[1], alignment_vec, reference_pack)
    del file_writer

    fst = True
    with open(sam_file_name, "r") as in_file:
        for line in in_file:
            if line[0] != '@':
                #print("|", line, "|", sep="")
                fields = line[:-1].split("\t")
                #print("|" + fields[10] + "|", "|" + qual + "|", sep='\n')
                assert fields[9] == query1 if fst else query2
                assert fields[10] == qual1 if fst else query2
                fst = False

for _ in range(100):
    query1 = ""
    for _ in range(100):
        query1 += choice(['C', 'T', 'G', 'A'])
        
    query2 = ""
    for _ in range(100):
        query2 += choice(['C', 'T', 'G', 'A'])

    # write a fasta file
    fasta_1_file_name = tempfile.gettempdir() + "/.tmp.1.fasta"
    fasta_2_file_name = tempfile.gettempdir() + "/.tmp.2.fasta"
    sam_file_name = tempfile.gettempdir() + "/.tmp.p.sam"
    with open(fasta_1_file_name, "w") as out_file:
        out_file.write(">name desc\n")
        out_file.write(query1)
        out_file.write("\n")

    with open(fasta_2_file_name, "w") as out_file:
        out_file.write(">name desc\n")
        out_file.write(query2)
        out_file.write("\n")

    file_reader = PairedFileReader(ParameterSetManager(),
                                   libMA.filePathVector([libMA.path(fasta_1_file_name)]),
                                   libMA.filePathVector([libMA.path(fasta_2_file_name)]))
    query_vec = file_reader.execute()
    alignment_vec = libMA.AlignmentVector()
    file_writer = PairedFileWriter(ParameterSetManager(), sam_file_name, reference_pack)
    file_writer.execute(query_vec[0], query_vec[1], alignment_vec, reference_pack)
    del file_writer

    fst = True
    with open(sam_file_name, "r") as in_file:
        for line in in_file:
            if line[0] != '@':
                #print("|", line, "|", sep="")
                fields = line[:-1].split("\t")
                #print("|" + fields[10] + "|", "|" + qual + "|", sep='\n')
                assert fields[9] == query1 if fst else query2
                assert fields[10] == "*"
                fst = False