from LAuS import *
import random
    
file_all = open("/mnt/ssd0/chrom/random/all", "w")

ref = NucSeq("")

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
    file = open("/mnt/ssd0/chrom/human/chr" + name + ".fna", "r")
    line = file.readline()#ignore first line
    #file_all.write(line)
    line = file.readline()
    while line:
        ref.append(line)
        #file_all.write(line)
        line = file.readline()
"""
file_all.close()



print("b")
ref_seq = BWAPack()
ref_seq.append("name", "no comment",ref)
ref_seq.store("/mnt/ssd0/chrom/random/pack")

print("c")
fm_index = FMIndex(ref_seq)
print("d")
fm_index.store("/mnt/ssd0/chrom/random/index")
print("e")
rev_ref.reverse()
rev_ref_seq = BWAPack()
rev_ref_seq.append("name", "no comment",rev_ref)
print("a")
rev_fm_index = FMIndex(rev_ref_seq)
print("f")
rev_fm_index.store("/mnt/ssd0/chrom/random/rev_index")