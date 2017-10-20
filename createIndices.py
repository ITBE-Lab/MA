from LAuS import *
import random
    
file_all = open("/mnt/ssd0/chrom/human/all", "w")
ref_seq = Pack()

"""
file_all.write(">blub/n")
seq = ""
for index in range(1000000):
    num = random.randint(0,3)
    if num == 0:
        seq += 'a'
        file_all.write("a")
    elif num == 2:
        seq += 'c'
        file_all.write("c")
    elif num == 3:
        seq += 't'
        file_all.write("t")
    else:
        seq += 'g'
        file_all.write("g")
    if index % 20 == 0:
        file_all.write("\n")
ref = NucSeq(seq)
rev_ref = NucSeq(seq)
"""

for index in range(1,25):
    name = str(index)
    if index == 23:
        name = "X"
    if index == 24:
        name = "Y"
    ref_seq.append_fasta_file("/mnt/ssd0/chrom/human/chr" + name + ".fna")
    file = open("/mnt/ssd0/chrom/human/chr" + name + ".fna", "r")
    line = file.readline()
    while line:
        file_all.write(line)
        line = file.readline()
    file.close()
file_all.close()

ref_seq.store("/mnt/ssd0/chrom/human/pack")

fm_index = FMIndex(ref_seq)
fm_index.store("/mnt/ssd0/chrom/human/index")