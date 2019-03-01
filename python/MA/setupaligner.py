from .aligner import *
from .alignmentPrinter import *


def test_aligner():
    configureFast()

    reference_pack = Pack()
    reference_pack.load("C:\\Users\\Markus\\Desktop\\MAdata\\eColi")
    reference_index = FMIndex()
    reference_index.load("C:\\Users\\Markus\\Desktop\\MAdata\\eColi")

    query = reference_pack.extract_from_to(10000, 10200)

    # initialize modules
    seeding_module = BinarySeeding()
    soc_module = StripOfConsideration()
    harm_module = Harmonization()
    nw_module = NeedlemanWunsch()
    printer_module = AlignmentPrinter()

    # execute modules
    seeds = seeding_module.execute(reference_index, query)
    print("seed set (q,r,l):")
    for seed in seeds.extract_seeds(reference_index, 100, 0, len(query), True):
        print(seed.start, seed.start_ref, seed.size)
    socs = soc_module.execute(seeds, query, reference_pack, reference_index)
    harm_seeds = harm_module.execute(socs, query, reference_index)
    alignments = nw_module.execute(harm_seeds, query, reference_pack)

    for index, alignment in enumerate(alignments):
        print("Alignment ", index, ":")
        printer_module.execute(alignment, query, reference_pack)

    print("[Test Successful]")


def quick_align(parameter_set, queries, pack, fm_index, output):
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
