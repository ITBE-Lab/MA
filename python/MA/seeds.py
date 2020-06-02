from MA import *
import argparse
import time

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_parser():
    parser = argparse.ArgumentParser(description='Compute and print seeds', epilog="Use 'seeds.py createIndex ...' to compute the Pack, FMD-index and Minimizer-Index for a reference genome (fasta file given by read_or_genome_fasta). Then use that index_prefix for running the various seeding techniques")
    parser.add_argument('seeding_technique', help="Seeding technique to be used",
                        choices=["createIndex", "MEMs", "SMEMs", "maxSpanning", "Minimizers", "MiniToMEM", "MiniToSMEM", "MiniToMaxSpanning", "ExtendPurge"])
    parser.add_argument('index_prefix', help="prefix of the FMD-index, Minimizer-index and Pack")
    parser.add_argument('read_or_genome_fasta', help="path of the fasta file with reads. If 'createIndex' is used, this parameter shall be the path to the genome")
    parser.add_argument('-l', help='minimal length of variable-sized seeds and seed length for minimizers [19]',
                        default="19", type=int)
    parser.add_argument('-w', help='window size for minimizers [10]', default="10", type=int)
    parser.add_argument('-o',
                        help='occurrence filter: Filter out all seeds that occur more than O times on the reference [200]',
                        default="200", type=int)
    parser.add_argument('-p', help='print the computed seeds to stdout [True]', default=True, type=str2bool)
    return parser

def load_reads(p_m, path):
    print("# loading reads...")
    start = time.time()
    reads = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(path)])).cpp_module.read_all()
    # the minimizer index needs the reads as strings 
    # we convert here so that the conversion time is not included in the mm seeding runtime
    str_reads = libMA.StringVector()
    for nuc_seq in reads:
        str_reads.append(str(nuc_seq))
    end = time.time()
    print("# loading reads took", end-start, "seconds")
    return reads, str_reads

def load_indices(p_m, path):
    print("# loading indices...")
    start = time.time()
    pack = Pack()
    pack.load(path)
    mm_index = libMA.MinimizerIndex(p_m, path + ".mmi")
    fm_index = FMIndex()
    fm_index.load(path)
    end = time.time()
    print("# loading indices took", end-start, "seconds")
    return fm_index, mm_index, pack

def print_seeds(reads, seeds_l):
    print("# printing seeds as query pos\tref pos\tlength. Reverse complement strand locations are in the interval "
          "(genome_size,2*genome_size], where the location 2*genome_size corresponds to the first nucleotide "
          "of the forward strand.")
    for read, seeds in zip(reads, seeds_l):
        print(">", read.name)
        for seed in seeds:
            print(seed.start, seed.start_ref, seed.size, sep='\t')

def compute_mems(p_m, reads, str_reads, fm_index, mm_index, pack):
    p_m.by_name("Seeding Technique").set(2)
    return BinarySeeding(p_m).cpp_module.seed(fm_index, reads)

def compute_smems(p_m, reads, str_reads, fm_index, mm_index, pack):
    p_m.by_name("Seeding Technique").set(1)
    return BinarySeeding(p_m).cpp_module.seed(fm_index, reads)

def compute_mss(p_m, reads, str_reads, fm_index, mm_index, pack):
    p_m.by_name("Seeding Technique").set(0)
    return BinarySeeding(p_m).cpp_module.seed(fm_index, reads)

def compute_mm(p_m, reads, str_reads, fm_index, mm_index, pack):
    return mm_index.seed(str_reads, pack)

def compute_mm_to_mem(p_m, reads, str_reads, fm_index, mm_index, pack):
    mm_seeds = mm_index.seed(str_reads, pack)
    return libMA.SeedLumping(p_m).cpp_module.lump(mm_seeds, reads, pack)

def compute_mm_to_smem(p_m, reads, str_reads, fm_index, mm_index, pack):
    mm_seeds = mm_index.seed(str_reads, pack)
    mems = libMA.SeedLumping(p_m).cpp_module.lump(mm_seeds, reads, pack)
    return libMA.MaxExtendedToSMEM(p_m).cpp_module.filter(mems)

def compute_mm_to_max_sp(p_m, reads, str_reads, fm_index, mm_index, pack):
    mm_seeds = mm_index.seed(str_reads, pack)
    mems = libMA.SeedLumping(p_m).cpp_module.lump(mm_seeds, reads, pack)
    return libMA.MaxExtendedToMaxSpanning(p_m).cpp_module.filter(mems)

def compute_mm_extend_purge(p_m, reads, str_reads, fm_index, mm_index, pack):
    mm_seeds = mm_index.seed(str_reads, pack)
    mems = libMA.SeedExtender(p_m).cpp_module.extend(mm_seeds, reads, pack)
    return libMA.SortRemoveDuplicates(p_m).cpp_module.filter(mems)

seeders = {
    "MEMs": compute_mems,
    "SMEMs": compute_smems,
    "maxSpanning": compute_mss,
    "Minimizers": compute_mm,
    "MiniToMEM": compute_mm_to_mem,
    "MiniToSMEM": compute_mm_to_smem,
    "MiniToMaxSpanning": compute_mm_to_max_sp,
    "ExtendPurge": compute_mm_extend_purge,
}

if __name__ == "__main__":
    args = get_parser().parse_args()
    p_m = ParameterSetManager()

    p_m.by_name("Number of Threads").set(1)
    p_m.by_name("Use all Processor Cores").set(False)
    p_m.by_name("Minimal Seed Length").set(args.l)
    p_m.by_name("Maximal Ambiguity").set(args.o) # occurrence filter setting for fmd-index
    p_m.by_name("Minimizer Size").set(args.l) # k
    p_m.by_name("Minimizer Window Size").set(args.w) # w

    if args.seeding_technique != "createIndex":
        fm_index, mm_index, pack = load_indices(p_m, args.index_prefix)
        mm_index.set_mid_occ(args.o) # occurrence filter setting for minimizer index
        reads, str_reads = load_reads(p_m, args.read_or_genome_fasta)
        print("# seeding...")
        start = time.time()
        seeds = seeders[args.seeding_technique](p_m, reads, str_reads, fm_index, mm_index, pack)
        end = time.time()
        print("# seeding took", end-start, "seconds")
        if args.p:
            print_seeds(reads, seeds)
    else:
        print("# creating indices and pack...")
        start = time.time()
        pack = Pack()
        pack.append_fasta_file(args.read_or_genome_fasta)
        pack.store(args.index_prefix)
        FMIndex(pack).store(args.index_prefix)
        libMA.MinimizerIndex(p_m, pack.contigSeqs(), pack.contigNames()).dump(args.index_prefix + ".mmi")
        end = time.time()
        print("# creating indices and pack took", end-start, "seconds")

