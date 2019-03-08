from .aligner import *
from .alignmentPrinter import *
from random import choice


def test_aligner():
    parameter_manager = ParameterSetManager()
    parameter_manager.set_selected("PacBio")

    # create a reference string
    reference = ""
    for _ in range(1000000):
        reference += choice(['C', 'T', 'G'])

    # create pack and fmd index from that string
    reference_pack = Pack()
    reference_pack.append(
        "sequence name", "sequence description", NucSeq(reference))
    reference_index = FMIndex(reference_pack)

    query = NucSeq(reference[100:125] + "A" + reference[126:150] +
                   "A" + reference[151:175] + "A" + reference[176:200])

    # initialize modules
    seeding_module = BinarySeeding(parameter_manager)
    soc_module = StripOfConsideration(parameter_manager)
    harm_module = Harmonization(parameter_manager)
    nw_module = NeedlemanWunsch(parameter_manager)
    printer_module = AlignmentPrinter(parameter_manager)

    # execute modules
    seeds = seeding_module.execute(reference_index, query)
    print("seed set (q,r,l):")
    for seed in seeds.extract_seeds(reference_index, 100, 0, len(query), True):
        print(seed.start, seed.start_ref, seed.size)
    socs = soc_module.execute(seeds, query, reference_pack, reference_index)
    harm_seeds = harm_module.execute(socs, query, reference_index)
    for i, seeds in enumerate(harm_seeds):
        print("harmonization seed set ", i, ": (q,r,l):")
        for seed in seeds:
            print(seed.start, seed.start_ref, seed.size)
    alignments = nw_module.execute(harm_seeds, query, reference_pack)
    for index, alignment in enumerate(alignments):
        print("Alignment ", index, ":")
        printer_module.execute(alignment, query, reference_pack)

    print("[Test Successful]")


def quick_align(parameter_set, pack, fm_index, queries=None, output = None, paired_queries=None, paired_output = None):
    fm_index_pledge = Pledge()
    fm_index_pledge.set(fm_index)

    pack_pledge = Pledge()
    pack_pledge.set(pack)

    module_lock = Lock(parameter_set)
    module_seeding = BinarySeeding(parameter_set)
    module_soc = StripOfConsideration(parameter_set)
    module_harm = Harmonization(parameter_set)
    module_dp = NeedlemanWunsch(parameter_set)
    module_mapping_qual = MappingQuality(parameter_set)

    res = VectorPledge()
    if not queries is None:
        assert(not output is None)
        for _ in range(parameter_set.get_num_threads()):
            locked_query = promise_me(module_lock, queries)
            seeds = promise_me(module_seeding, fm_index_pledge, locked_query)
            socs = promise_me(module_soc, seeds, locked_query,
                              pack_pledge, fm_index_pledge)
            harm = promise_me(module_harm, socs, locked_query, fm_index_pledge)
            alignments = promise_me(module_dp, harm, locked_query, pack_pledge)
            alignments_w_map_q = promise_me(
                module_mapping_qual, locked_query, alignments)
            empty = promise_me(output, locked_query,
                               alignments_w_map_q, pack_pledge)
            unlock = promise_me(UnLock(parameter_set, locked_query), empty)
            res.append(unlock)
    if not paired_queries is None:
        assert(not paired_output is None)
        module_paired = PairedReads(parameter_set)
        module_get_first = GetFirstQuery(parameter_set)
        module_get_second = GetSecondQuery(parameter_set)
        for _ in range(parameter_set.get_num_threads()):
            locked_query_pair = promise_me(module_lock, paired_queries)
            # primary read
            locked_query = promise_me(module_get_first, locked_query_pair)
            seeds = promise_me(module_seeding, fm_index_pledge, locked_query)
            socs = promise_me(module_soc, seeds, locked_query,
                              pack_pledge, fm_index_pledge)
            harm = promise_me(module_harm, socs, locked_query, fm_index_pledge)
            alignments = promise_me(module_dp, harm, locked_query, pack_pledge)
            alignments_w_map_q = promise_me(
                module_mapping_qual, locked_query, alignments)

            # mate read
            locked_query_mate = promise_me(
                module_get_second, locked_query_pair)
            seeds_mate = promise_me(
                module_seeding, fm_index_pledge, locked_query_mate)
            socs_mate = promise_me(module_soc, seeds_mate, locked_query_mate,
                                   pack_pledge, fm_index_pledge)
            harm_mate = promise_me(module_harm, socs_mate,
                                   locked_query_mate, fm_index_pledge)
            alignments_mate = promise_me(
                module_dp, harm_mate, locked_query_mate, pack_pledge)
            alignments_w_map_q_mate = promise_me(
                module_mapping_qual, locked_query_mate, alignments_mate)

            # combine & output
            combined_alignments = promise_me(
                module_paired, locked_query, locked_query_mate, alignments_w_map_q, alignments_w_map_q_mate, pack_pledge)
            empty = promise_me(paired_output, locked_query,
                               locked_query_mate, combined_alignments, pack_pledge)

            # unlock
            unlock = promise_me(
                UnLock(parameter_set, locked_query_pair), empty)
            res.append(unlock)

    res.simultaneous_get(parameter_set.get_num_threads())


def quick_align_paths(queries, genome_prefix, parameter_set_manager=ParameterSetManager(),
                      output_path="alignments.sam"):
    pack = Pack()
    pack.load(genome_prefix)
    fm_index = FMIndex()
    fm_index.load(genome_prefix)

    query_vec_pledge = Pledge()
    query_vec_pledge.set(ContainerVectorNucSeq(queries))

    module_splitter = libMA.NucSeqSplitter(parameter_set_manager)
    queries_pledge = promise_me(module_splitter, query_vec_pledge)

    module_output = FileWriter(parameter_set_manager, output_path, pack)

    quick_align(parameter_set_manager, queries_pledge,
                pack, fm_index, module_output)

# queries = [MA.NucSeq("ACCCGTGTGTACGACTACGGCATCAGACTACGAC")]
# MA.quick_align_paths(queries, "/mnt/ssd0/genome/GRCh38.p12")
