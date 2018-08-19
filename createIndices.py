from MA import *
import random
import os


def replace_n(filenames, out_file):
    file_out = open(out_file, "w")
    for name in filenames:
        file = open(name , "r")
        line = file.readline()
        while line:
            line = line[:-1]
            if line[0] != '>':
                i = 0
                while i < len(line):
                    if not line[i] in ["a", "c", "g", "t", "A", "C", "G", "T"]:
                        line = line[:i] + line[i+1:] #remove the N
                    else:
                        i += 1
                assert(not 'N' in line)
                assert(not 'n' in line)
            else:
                print(line[0:-1])
            if len(line) > 1:
                file_out.write(line + "\n")
            line = file.readline()
        file.close()
    file_out.close()
    print("self check...")
    file = open(out_file , "r")
    line = file.readline()
    while line:
        if line[0] != '>':
            if "N" in line or 'n' in line:
                print(line)
            #assert(len(line) > 1)
        line = file.readline()
    file.close()
    print("[OK]")

def concat(filenames, out_file):
    file_out = open(out_file, "w")
    for name in filenames:
        file = open(name , "r")
        line = file.readline()
        while line:
            line = line[:-1]
            if len(line) > 1:
                file_out.write(line + "\n")
            line = file.readline()
        file.close()
    file_out.close()
    print("self check...")
    file = open(out_file , "r")
    line = file.readline()
    while line:
        if line[0] != '>':
            if "N" in line or 'n' in line:
                print(line)
            #assert(len(line) > 1)
        line = file.readline()
    file.close()
    print("[OK]")

bowtie2_home = "/usr/home/markus/workspace/bowtie2/bowtie2-2.3.3.1/bowtie2-build "
blasr_home = "/usr/home/markus/workspace/blasr/build/bin/sawriter "
bwa_home = "/usr/home/markus/workspace/bwa/bwa index "
minimap_home = "/usr/home/markus/workspace/minimap2/minimap2 -d "
gem_home = "/usr/home/markus/workspace/gemtools/GEMTools/bin/"
def make(filename, out_name):
    # my pack + fmd Index
    ref_seq = Pack()
    ref_seq.append_fasta_file( filename )
    ref_seq.store(out_name)
    fm_index = FMIndex(ref_seq)
    fm_index.store(out_name)
    # BWA fmd index
    os.system(bwa_home + "-p " + out_name + "bwa " + filename)
    # bowtie index
    os.system(bowtie2_home + filename + " " + out_name + "bowtie2")
    # blasr index
    os.system(blasr_home + out_name + "blasr " + filename)
    # minimap2 index
    os.system(minimap_home + out_name + ".mmi " + filename)

    # # GEM index
    # os.system(gem_home + "gt.construct -i " + filename + " -o " + out_name + ".gem")

""" # DEPRECATED #
def make_hash(name):
    print("loading pack...")
    ref_seq = Pack()
    ref_seq.load(name)
    print("done loading pack\nextracting sequence...")
    #seq = ref_seq.extract_from_to(500, 600)
    #ref_seq = Pack()
    #ref_seq.append("a", "b", seq)
    #ref_seq.store("/MAdata/genome/minimal")
    #exit()
    seq = ref_seq.extract_forward_strand()
    print("done extracting sequence")
    #seq = NucSeq("ACCCCTGTGTTGTCACATCGATACGACTACGACACATCAGCACTACGACTACACCCCTGTGTTGTCACATCGATACGACTACGACACATCAGCACTACGACTACACCCCTGTGTTGTCACATCGATACGACTACGACACATCAGCACTACGACTACACCCCTGTGTTGTCACATCGATACGACTACGACACATCAGCACTACGACTACACCCCTGTGTTGTCACATCGATACGACTACGACACATCAGCACTACGACTACACCCCTGTGTTGTCACATCGATACGACTACGACACATCAGCACTACGACTACACCCCTGTGTTGTCACATCGATACGACTACGACACATCAGCACTACGACTAC")
    minimizer = Minimizers()
    minimizer.print = True
    hash_index = minimizer.execute(seq).toHash(ref_seq.unpacked_size(), ref_seq)
    print("saving...")
    hash_index.to_file(name + ".maRef")
    print("done saving")
"""

def makeRandom(out_name, size):
    file_all = open(out_name, "w")

    file_all.write(">randomSequence0 description/n")
    for index in range(size):
        num = random.randint(0,3)
        if num == 0:
            file_all.write("A")
        elif num == 2:
            file_all.write("C")
        elif num == 3:
            file_all.write("T")
        else:
            file_all.write("G")
        if index % 50 == 0:
            file_all.write("\n")
            if index % 50000 == 0:
                file_all.write(">randomSequence" + str(index) + " description/n")
    file_all.close()

def chrNames(prefix, num, suffix):
    ret = []
    for index in range(0,num):
        name = str(index+1)# files count from 1 to num (inclusive)
        ret.append(prefix + name + suffix)
    ret.append(prefix + "X" + suffix)
    ret.append(prefix + "Y" + suffix)
    return ret

#replace_n(["/MAdata/chrom/zebrafish/GCF_000002035.6_GRCz11_genomic.fna"], "/MAdata/chrom/zebrafish/n_free.fasta")
#make("/MAdata/chrom/zebrafish/n_free.fasta", "/MAdata/genome/zebrafish")

#replace_n(["/MAdata/chrom/wheat/full_genome.fasta"], "/MAdata/chrom/wheat/n_free.fasta")
#make("/MAdata/chrom/wheat/n_free.fasta", "/MAdata/genome/wheat")

#replace_n(
#        chrNames("/MAdata/chrom/human_hg37/chr", 22, ".fa"),
#        "/MAdata/chrom/human_hg37/n_free.fasta"
#    )
#concat(
#        chrNames("/MAdata/chrom/human_hg37/chr", 22, ".fa"),
#        "/MAdata/chrom/human_hg37/full_sequence.fasta"
#    )
#make("/MAdata/chrom/eColi/GCA_000005845.2_ASM584v2_genomic.fna", "/MAdata/genome/eColi_full")
make("/MAdata/chrom/zebrafish/GCF_000002035.6_GRCz11_genomic.fna", "/MAdata/genome/zebrafish_full")

#replace_n(chrNames("/MAdata/chrom/mouse/chr", 19, ".fna"), "/MAdata/chrom/mouse/n_free.fasta")
#make("/MAdata/chrom/mouse/n_free.fasta", "/MAdata/genome/mouse")

#replace_n(["/MAdata/chrom/zebrafish/genome.fasta"], "/MAdata/chrom/zebrafish/n_free.fasta")
#make("/MAdata/chrom/plasmodium/n_free.fasta", "/MAdata/genome/plasmodium")

#makeRandom("/MAdata/chrom/random.fasta", 3 * 10**9)
#make("/MAdata/chrom/zebrafish/n_free.fasta", "/MAdata/genome/zebrafish")
#make_hash("/MAdata/genome/human")
#makeRandom("/MAdata/genome/random_3_10_7", 3 * 10**7)

#replace_n(["/MAdata/chrom/eColi/GCA_000005845.2_ASM584v2_genomic.fna"], "/MAdata/chrom/eColi/n_free.fasta")
#make("/MAdata/chrom/eColi/n_free.fasta", "/MAdata/genome/eColi")