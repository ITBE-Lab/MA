from .aligner import *

def insert_paired_reads(parameter_set, database_name, sequencer_name, filename_vec1, filename_vec2):

    get_read_inserter = GetPairedReadInserter(parameter_set, sequencer_name)
    inserter_module = ReadInserterModule(PairedReadInserterModule)

    return


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
    modlue_inv = SmallInversions(parameter_set)

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

            if parameter_set.by_name("Detect Small Inversions").get():
                alignments_w_map_q = promise_me(
                    modlue_inv, alignments_w_map_q, locked_query, pack_pledge)

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

            if parameter_set.by_name("Detect Small Inversions").get():
                alignments_w_map_q = promise_me(
                    modlue_inv, alignments_w_map_q, locked_query, pack_pledge)
                alignments_w_map_q_mate = promise_me(
                    modlue_inv, alignments_w_map_q_mate, locked_query_mate, pack_pledge)

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