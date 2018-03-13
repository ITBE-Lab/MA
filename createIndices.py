from MA import *
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

def make_hash(name):
    print("loading pack...")
    ref_seq = Pack()
    ref_seq.load(name)
    print("done loading pack\nextracting sequence...")
    #seq = ref_seq.extract_from_to(500, 600)
    #ref_seq = Pack()
    #ref_seq.append("a", "b", seq)
    #ref_seq.store("/mnt/ssd0/genome/minimal")
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
        if index % 50000 == 0:
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


#make(chrNames("/mnt/ssd0/chrom/human/chr", 22, ".fna"), "/mnt/ssd0/genome/human")
#make(["/mnt/ssd0/chrom/human/chr1.fna" ], "/mnt/ssd0/genome/humanchr1")
#make_hash("/mnt/ssd0/genome/human")
#make(["/mnt/ssd0/chrom/plasmodium/genome.fasta"], "/mnt/ssd0/genome/plasmodium")
#make(chrNames("/mnt/ssd0/chrom/mouse/chr", 21, ".fna"), "/mnt/ssd0/chrom/mouse/all")
#makeRandom("/mnt/ssd0/genome/random", 3 * 10**9)
#makeRandom("/mnt/ssd0/genome/random_3_10_7", 3 * 10**7)