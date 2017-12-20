from MABS import *
import random

def make(filenames, out_name):
    file_all = open(out_name + "bwa", "w")
    ref_seq = Pack()

    for name in filenames:
        ref_seq.append_fasta_file( name )
        file = open(name , "r")
        line = file.readline()
        while line:
            file_all.write(line)
            line = file.readline()
        file.close()
    file_all.close()

    ref_seq.store(out_name)

    fm_index = FMIndex(ref_seq)
    fm_index.store(out_name)

def makeRandom(out_name, size):
    file_all = open(out_name + "bwa", "w")
    ref_seq = Pack()

    file_all.write(">random sequence/n")
    seq = ""
    for index in range(size):
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
    ref_seq.append( "random sequence", "no desc", NucSeq(seq) )
    file_all.close()

    ref_seq.store(out_name)

    fm_index = FMIndex(ref_seq)
    fm_index.store(out_name)

def chrNames(prefix, num, suffix):
    ret = []
    for index in range(0,num):
        name = str(index+1)# files count from 1 to num (inclusive)
        ret.append(prefix + name + suffix)
    ret.append(prefix + "X" + suffix)
    ret.append(prefix + "Y" + suffix)
    return ret


make(chrNames("/mnt/ssd0/chrom/human/chr", 22, ".fna"), "/mnt/ssd0/genome/human")
#make(
#    ["/mnt/ssd0/chrom/human/GCF_000001405.37_GRCh38.p11_genomic.fna"], "/mnt/ssd0/genome/human")
#make(chrNames("/mnt/ssd0/chrom/mouse/chr", 21, ".fna"), "/mnt/ssd0/chrom/mouse/all")
#makeRandom("/mnt/ssd0/chrom/random/random", 1000000)