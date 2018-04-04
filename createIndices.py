from MA import *
import random

def replace_n(filenames, out_file):
    file_out = open(out_file, "w")
    for name in filenames:
        file = open(name , "r")
        line = file.readline()
        while line:
            if line[0] != '>':
                i = 0
                while i < len(line):
                    if line[i] == 'N' or line[i] == 'n':
                        line = line[:i] + line[i+1:] #remove the N
                    else:
                        i += 1
            else:
                print(line[0:-1])
            if len(line) > 1:
                file_out.write(line)
            line = file.readline()
        file.close()
    file_out.close()
    print("self check...")
    file = open(out_file , "r")
    line = file.readline()
    while line:
        if line[0] != '>':
            assert(not 'N' in line)
            assert(not 'n' in line)
            assert(len(line) > 1)
        line = file.readline()
    file.close()
    print("[OK]")

bowtie2_home = "/usr/home/markus/workspace/bowtie2/bowtie2-2.3.3.1/bowtie2-build "
blasr_home = "/usr/home/markus/workspace/blasr/build/bin/sawriter "
bwa_home = "/usr/home/markus/workspace/bwa/bwa index "
def make(filename, out_name):
    # my pack + fmd Index
    #ref_seq = Pack()
    #ref_seq.append_fasta_file( filename )
    #ref_seq.store(out_name)

    #fm_index = FMIndex(ref_seq)
    #fm_index.store(out_name)

    #BWA fmd index
    os.system(bwa_home + "-p " + out_name + "bwa " + filename)
    # #bowtie index
    #os.system(bowtie2_home + filename + " " + out_name + "bowtie2")
    # #blasr index
    #os.system(blasr_home + out_name + "blasr " + filename)

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


def makeRandom(out_name, size):
    file_all = open(out_name + "bwa", "w")
    ref_seq = Pack()

    file_all.write(">randomSequence/n")
    seq = ""
    for index in range(size):
        num = random.randint(0,3)
        if num == 0:
            seq += 'A'
            file_all.write("A")
        elif num == 2:
            seq += 'C'
            file_all.write("C")
        elif num == 3:
            seq += 'T'
            file_all.write("T")
        else:
            seq += 'G'
            file_all.write("G")
        if index % 20 == 0:
            file_all.write("\n")
        if index % 10000 == 0:
            ref_seq.append( "randomSequence" + str(index), "noDesc", NucSeq(seq) )
            seq = ""
    if seq != "":
        ref_seq.append( "randomSequenceLast", "noDesc", NucSeq(seq) )
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


#replace_n(chrNames("/MAdata/chrom/human/chr", 22, ".fna"), "/MAdata/chrom/human/n_free.fasta")
make("/MAdata/chrom/human/n_free.fasta", "/MAdata/genome/human")

#replace_n(chrNames("/MAdata/chrom/mouse/chr", 19, ".fna"), "/MAdata/chrom/mouse/n_free.fasta")
make("/MAdata/chrom/mouse/n_free.fasta", "/MAdata/genome/mouse")

#replace_n(["/MAdata/chrom/plasmodium/genome.fasta"], "/MAdata/chrom/plasmodium/n_free.fasta")
make("/MAdata/chrom/plasmodium/n_free.fasta", "/MAdata/genome/plasmodium")

#makeRandom("/MAdata/genome/random", 3 * 10**9)
#make(["/MAdata/chrom/human/chr1.fna" ], "/MAdata/genome/humanchr1_bugged")
#make_hash("/MAdata/genome/human")
#makeRandom("/MAdata/genome/random_3_10_7", 3 * 10**7)