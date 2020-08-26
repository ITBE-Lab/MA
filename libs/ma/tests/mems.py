from MA import *
import random

for _ in range(1000):
    def random_nuc_seq(l):
        ret = ""
        for _ in range(l):
            ret += random.choice(['a', 'c', 'g', 't'])
        return NucSeq(ret)

    min_size = 3

    p_m = ParameterSetManager()
    p_m.by_name("Seeding Technique").set(2)
    p_m.by_name("Minimizers - k").set(min_size)
    p_m.by_name("Minimizers - w").set(1)
    p_m.by_name("Minimal Seed Length").set(min_size-2)


    query = random_nuc_seq(20)

    pack = Pack()
    pack.append("chr1", "chr1-desc", random_nuc_seq(10))
    fm_index = FMIndex(pack)

    fm_seeds = BinarySeeding(p_m).execute(fm_index, query).extract_seeds(fm_index, 1000, min_size - 1, len(query), True)


    mm_index = MinimizerIndex(p_m, pack.contigSeqs(), pack.contigNames())

    reads = StringVector()
    reads.append(str(query))
    queries = VecRetNuc()
    queries.append(query)
    mm_seeds = SeedLumping(p_m).cpp_module.lump(mm_index.seed(reads, True, pack), queries, pack)

    for i in range(10):
        print(i, end="")
    for i in range(10):
        print(i, end="")
    print()
    print(str(query))
    print(str(pack.extract_complete()))



    seed_order = lambda x: (x.start, x.start_ref, x.size)
    def sorted_n_filtered(seeds):
        def filter(seeds):
            for seed in seeds:
                if seed.on_forward_strand:
                    if not pack.is_bridging(seed.start_ref, seed.size + 1):
                        yield seed
                else:
                    if not pack.is_bridging(seed.start_ref + 1 - seed.size, seed.size + 1):
                        yield seed
        return sorted(filter(seeds), key=seed_order)

    mm_filterd = list(sorted_n_filtered(mm_seeds[0]))
    fm_filterd = list(sorted_n_filtered(fm_seeds))

    print("mm_seed")
    for seed in mm_filterd:
        print("qrl:", seed.start, seed.start_ref, seed.size)

    print("fm_seeds")
    for seed in fm_filterd:
        print("qrl:", seed.start, seed.start_ref, seed.size)

    for a, b in zip(fm_filterd, fm_filterd):
        if not seed_order(a) == seed_order(b):
            print("seeds do not match!")
            print("a = qrl:", a.start, a.start_ref, a.size)
            print("b = qrl:", b.start, b.start_ref, b.size)
            assert False

    assert len(fm_filterd) == len(mm_filterd)
