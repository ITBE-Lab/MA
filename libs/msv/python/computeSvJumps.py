from MSV import *
from MS import *
from MA import *
import datetime

def compute_sv_jumps(parameter_set_manager, mm_index, pack, dataset_name, seq_ids=0, runtime_file=None):
    #parameter_set_manager.by_name("Number of Threads").set(1)
    #parameter_set_manager.by_name("Use all Processor Cores").set(False)
    #assert parameter_set_manager.get_num_threads() == 1

    parameter_set_manager.by_name("Fixed SoC Width").set(50)
    parameter_set_manager.by_name("Max Size Reseed").set(2000)
    parameter_set_manager.by_name("Maximal Ambiguity").set(1)
    parameter_set_manager.by_name("Min Size Edge").set(10) # runtime optimization...
    parameter_set_manager.by_name("Min NT in SoC").set(50)
    parameter_set_manager.by_name("Rectangular SoC").set(False)
    mm_index.set_max_occ(2)
    def scope():
        single_con = DbConn(dataset_name)

        nuc_seq_fetcher = NucSeqFetcher(parameter_set_manager)
        lock_module = Lock(parameter_set_manager)
        seeding_module = MMFilteredSeeding(parameter_set_manager)
        seed_lumper = SeedLumping(parameter_set_manager)
        contig_filter = FilterContigBorder(parameter_set_manager)
        soc_module = StripOfConsiderationSeeds(parameter_set_manager)
        soc_filter = GetAllFeasibleSoCs(parameter_set_manager)
        jumps_from_seeds = SvJumpsFromExtractedSeeds(parameter_set_manager, pack)
        filter_by_ambiguity = FilterJumpsByRefAmbiguity(parameter_set_manager)
        get_jump_inserter = GetJumpInserter(parameter_set_manager, single_con, "MS-SV",
                                            "python built comp graph")
        jump_inserter_module = JumpInserterModule(parameter_set_manager)

        mm_pledge = Pledge()
        mm_pledge.set(mm_index)
        pack_pledge = Pledge()
        pack_pledge.set(pack)
        pool_pledge = Pledge()
        pool_pledge.set(PoolContainer(parameter_set_manager.get_num_threads() + 1, dataset_name))


        if not isinstance(seq_ids, list):
            _seq_ids = [seq_ids]
        else:
            _seq_ids = seq_ids

        jobs = parameter_set_manager.get_num_threads()
        for seq_id in _seq_ids:
            mm_counter_container = HashFilterTable(single_con).get_counter(seq_id)
            mm_counter = Pledge()
            mm_counter.set(mm_counter_container)

            analyze = AnalyzeRuntimes()
            res = VectorPledge()
            inserter_vec = []
            # graph for single reads
            for idx in parallel_graph(jobs):
                # @todo there should be a set of modules distributing reads
                # (only a problem if threads finish at different times)...
                nuc_seq_query_getter = GetNucSeqFromSqlQuery(parameter_set_manager, seq_id, idx, jobs, True, True)
                nuc_seq_query = promise_me(nuc_seq_query_getter, pool_pledge)
                analyze.register("GetNucSeqFromSqlQuery", nuc_seq_query, True)
                queries_pledge = promise_me(nuc_seq_fetcher, pool_pledge, nuc_seq_query)
                analyze.register("NucSeqFetcher", queries_pledge, True)
                query_pledge = promise_me(lock_module, queries_pledge)
                analyze.register("Lock", query_pledge, True)
                seeds_pledge = promise_me(seeding_module, mm_pledge, query_pledge, pack_pledge, mm_counter)
                analyze.register("MinimizerSeeding", seeds_pledge, True)
                lumped_seeds = promise_me(seed_lumper, seeds_pledge, query_pledge, pack_pledge)
                analyze.register("SeedLumping", lumped_seeds, True)
                c_filter_seeds = promise_me(contig_filter, lumped_seeds, pack_pledge)
                analyze.register("FilterContigBorder", c_filter_seeds, True)
                socs = promise_me(soc_module, c_filter_seeds, query_pledge, pack_pledge)
                analyze.register("SoC", socs, True)
                filtered_seeds_pledge_2 = promise_me(soc_filter, socs)
                analyze.register("SoCFilter", filtered_seeds_pledge_2, True)
                jumps_pledge = promise_me(jumps_from_seeds, filtered_seeds_pledge_2, pack_pledge, query_pledge)
                analyze.register("SvJumpsFromSeeds", jumps_pledge, True)
                filtered_jumps_pledge = promise_me(filter_by_ambiguity, jumps_pledge, pack_pledge)
                analyze.register("FilterJumpsByRefAmbiguity", filtered_jumps_pledge, True)
                jump_inserter = promise_me(get_jump_inserter, pool_pledge)
                inserter_vec.append(jump_inserter)
                analyze.register("GetJumpInserter", jump_inserter, True)
                write_to_db_pledge = promise_me(jump_inserter_module, jump_inserter, pool_pledge, filtered_jumps_pledge,
                                                query_pledge)
                analyze.register("JumpInserterModule", write_to_db_pledge, True)
                unlock_pledge = promise_me(UnLock(parameter_set_manager, query_pledge), write_to_db_pledge)
                analyze.register("UnLock", unlock_pledge, True)
                res.append(unlock_pledge)

            # drain all sources
            res.simultaneous_get(parameter_set_manager.get_num_threads())
            for inserter in inserter_vec:
                inserter.get().close(pool_pledge.get()) # @todo for some reason the destructor does not trigger automatically :(

            analyze.analyze(runtime_file)

        return get_jump_inserter.cpp_module.id

    db_conn = DbConn(dataset_name)
    JumpRunTable(db_conn)

    jump_table = SvJumpTable(db_conn)
    jump_table.drop_indices(0) # number does nothing at the moment

    jump_id = scope()

    analyze = AnalyzeRuntimes()

    print("num jumps:", jump_table.num_jumps(jump_id))

    print("creating index...")
    start = datetime.datetime.now()
    jump_table.create_indices( jump_id )
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("create_indices", delta.total_seconds(), False, lambda x: x)
    print("created index")
    analyze.analyze(runtime_file)

    if not runtime_file is None:
        runtime_file.write("sv_jump_run_id is " + str(jump_id) + "\n")

    # return the run_id
    return jump_id

