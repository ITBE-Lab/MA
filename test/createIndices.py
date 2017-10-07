from LAuS import *
import random
    
#file_all = open("/mnt/ssd0/chrom/human/all", "w")

ref = NucSeq("")
"""
seq = ""
for _ in range(10000):
    num = random.randint(0,3)
    if num == 0:
        seq += 'a'
    elif num == 2:
        seq += 'c'
    elif num == 3:
        seq += 't'
    else:
        seq += 'g'
ref = NucSeq(seq)

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

#file_all.write(">blub")
#file_all.write(seq)
#file_all.close()



print("b")
#ref_seq = BWAPack()
#ref_seq.append("name", "no comment",ref)
#ref_seq.store("/mnt/ssd0/chrom/human/pack")

print("c")
#fm_index = FMIndex(ref_seq)
print("d")
#fm_index.store("/mnt/ssd0/chrom/human/index")
print("e")
ref.reverse()
rev_ref_seq = BWAPack()
rev_ref_seq.append("name", "no comment",ref)
print("a")
rev_fm_index = FMIndex(rev_ref_seq)
print("f")
rev_fm_index.store("/mnt/ssd0/chrom/human/rev_index")