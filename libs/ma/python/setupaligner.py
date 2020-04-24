from .alignmentPrinter import *
from random import choice
from MA import *

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


def quick_align(parameter_set, pack, fm_index, output_graph, file_queue, file_queue_2=None):
    fm_index_pledge = Pledge()
    fm_index_pledge.set(fm_index)

    pack_pledge = Pledge()
    pack_pledge.set(pack)

    module_seeding = BinarySeeding(parameter_set)
    module_soc = StripOfConsideration(parameter_set)
    module_harm = Harmonization(parameter_set)
    module_dp = NeedlemanWunsch(parameter_set)
    module_mapping_qual = MappingQuality(parameter_set)
    module_inv = SmallInversions(parameter_set)

    combined_queue = file_queue
    lock = Lock(parameter_set)
    if not file_queue_2 is None:
        combined_queue = combine_file_streams(file_queue, file_queue_2)
        module_get_first = GetFirstQuery(parameter_set)
        module_get_second = GetSecondQuery(parameter_set)
        module_paired = PairedReads(parameter_set)

        #set up modules for paired read insert
        queue_picker = PairedFilePicker(parameter_set)
        queue_placer = PairedFileNucSeqPlacer(parameter_set)
        file_reader = PairedFileReader(parameter_set)
    else:
        #set up modules for single read insert
        queue_picker = FilePicker(parameter_set)
        queue_placer = FileNucSeqPlacer(parameter_set)
        file_reader = FileReader(parameter_set)
    
    queue_pledge = Pledge()
    queue_pledge.set(combined_queue)

    def query_to_alignment(query):
        seeds = promise_me(module_seeding, fm_index_pledge, query)
        socs = promise_me(module_soc, seeds, query,
                            pack_pledge, fm_index_pledge)
        harm = promise_me(module_harm, socs, query, fm_index_pledge)
        alignments = promise_me(module_dp, harm, query, pack_pledge)
        alignments_w_map_q = promise_me(
            module_mapping_qual, query, alignments)
        if parameter_set.by_name("Detect Small Inversions").get():
            alignments_w_map_q = promise_me(
                module_inv, alignments_w_map_q, query, pack_pledge)
        return alignments_w_map_q

    res = VectorPledge()
    for _ in range(parameter_set.get_num_threads()):
        picked_file = promise_me(queue_picker, queue_pledge)
        locked_file = promise_me(lock, picked_file)
        query_ = promise_me(file_reader, locked_file)

        locked_query = promise_me(queue_placer, query_, locked_file, queue_pledge)
        if not file_queue_2 is None:
            query_prim = promise_me(module_get_first, locked_query)
            alignment_pledge = query_to_alignment(query_prim)

            query_mate = promise_me(module_get_second, locked_query)
            alignment_pledge_mate = query_to_alignment(query_mate)

            alignment_pledge_combined = promise_me(
                module_paired, query_prim, query_mate, alignment_pledge, alignment_pledge_mate, pack_pledge)
            empty = output_graph(alignment_pledge_combined, pack_pledge, query_prim, query_mate)
        else:
            alignment_pledge = query_to_alignment(locked_query)
            empty = output_graph(alignment_pledge, pack_pledge, locked_query)

        unlock = promise_me(UnLock(parameter_set, locked_file), empty)
        res.append(unlock)

    res.simultaneous_get(parameter_set.get_num_threads())


def quick_align_paths(queries, genome_prefix, parameter_set_manager=ParameterSetManager(),
                      output_path="alignments.sam"):
    def to_file_queue(string_vec):
        if string_vec is None:
            return None
        file_queue = FileQueue()
        for string in string_vec:
            file_queue.add(FileStreamFromPath(string))
        return file_queue
    pack = Pack()
    pack.load(genome_prefix)
    fm_index = FMIndex()
    fm_index.load(genome_prefix)

    file_writer = FileWriter(parameter_set_manager, output_path, pack)
    def output_graph(alignment_pledge, pack_pledge, locked_query):
        return promise_me(file_writer,
                          locked_query, alignment_pledge, pack_pledge)

    quick_align(parameter_set_manager, pack, fm_index, output_graph, to_file_queue(queries))
