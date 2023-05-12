from MSV import *
from MS import *
from MA import *
import datetime

def compute_sv_jumps(parameter_set_manager, mm_index, pack, dataset_name, seq_ids=0, runtime_file=None):
    #parameter_set_manager.by_name("Number of Threads").set(1)
    #parameter_set_manager.by_name("Use all Processor Cores").set(False)
    #assert parameter_set_manager.get_num_threads() == 1

    mm_index.set_max_occ(2)
    def scope():
        single_con = DbConn(dataset_name)

        nuc_seq_fetcher = NucSeqFetcher(parameter_set_manager)
        lock_module = Lock(parameter_set_manager)
        seeding_module = MMFilteredSeeding(parameter_set_manager)
        seed_lumper = SeedLumping(parameter_set_manager)
        soc_module = StripOfConsiderationSeeds(parameter_set_manager)
        soc_filter = GetAllFeasibleSoCsAsSet(parameter_set_manager)
        reseeding = RecursiveReseedingSoCs(parameter_set_manager, pack)
        jumps_from_seeds = SvJumpsFromExtractedSeeds(parameter_set_manager, pack)
        contig_filter = JumpsFilterContigBorder(parameter_set_manager)
        #filter_by_ambiguity = FilterJumpsByRefAmbiguity(parameter_set_manager)
        get_jump_inserter = GetJumpInserter(parameter_set_manager, single_con, "MS-SV",
                                            "python built comp graph")
        jump_inserter_module = JumpInserterModule(parameter_set_manager)

        #cov_table.init_coverage(get_jump_inserter.cpp_module.id, 
        #            pack.unpacked_size_single_strand // parameter_set_manager.by_name("Coveragebin size").get())
        #get_coverage_updater = GetCoverageUpdater(parameter_set_manager, get_jump_inserter.cpp_module.id)
        #coverage_updater_module = CoverageUpdaterModule(parameter_set_manager)
        cc_module = CollectSeedCoverage(parameter_set_manager)

        join_module = ContainerJoin(parameter_set_manager)

        cc_col = Pledge()
        cc_col.set(CoverageCollector(pack.unpacked_size_single_strand))
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
                unlocked_queries_pledge = promise_me(nuc_seq_fetcher, pool_pledge, nuc_seq_query)
                analyze.register("NucSeqFetcher", unlocked_queries_pledge, True)
                query_pledge = promise_me(lock_module, unlocked_queries_pledge)
                analyze.register("Lock", query_pledge, True)
                seeds_pledge = promise_me(seeding_module, mm_pledge, query_pledge, pack_pledge, mm_counter)
                analyze.register("MinimizerSeeding", seeds_pledge, True)
                lumped_seeds = promise_me(seed_lumper, seeds_pledge, query_pledge, pack_pledge)
                analyze.register("SeedLumping", lumped_seeds, True)
                socs = promise_me(soc_module, lumped_seeds, query_pledge, pack_pledge)
                analyze.register("SoC", socs, True)
                filtered_seeds_pledge_2 = promise_me(soc_filter, socs)
                analyze.register("SoCFilter", filtered_seeds_pledge_2, True)
                filtered_seeds_pledge_3 = promise_me(reseeding, filtered_seeds_pledge_2, pack_pledge, query_pledge)
                analyze.register("RecursiveReseedingSoCs", filtered_seeds_pledge_3, True)
                jumps_pledge = promise_me(jumps_from_seeds, filtered_seeds_pledge_3, pack_pledge, query_pledge)
                analyze.register("SvJumpsFromSeeds", jumps_pledge, True)
                filtered_jumps = promise_me(contig_filter, jumps_pledge, pack_pledge)
                analyze.register("FilterContigBorder", filtered_jumps, True)
                #filtered_jumps_pledge = promise_me(filter_by_ambiguity, jumps_pledge, pack_pledge)
                #analyze.register("FilterJumpsByRefAmbiguity", filtered_jumps_pledge, True)
                jump_inserter = promise_me(get_jump_inserter, pool_pledge)
                inserter_vec.append(jump_inserter)
                analyze.register("GetJumpInserter", jump_inserter, True)
                write_to_db_pledge = promise_me(jump_inserter_module, jump_inserter, pool_pledge, filtered_jumps,
                                                query_pledge)
                analyze.register("JumpInserterModule", write_to_db_pledge, True)
                #cov_updater = promise_me(get_coverage_updater, pool_pledge)
                #inserter_vec.append(cov_updater)
                #analyze.register("GetCovUpdater", cov_updater, True)
                #write_to_db_pledge_2 = promise_me(coverage_updater_module, cov_updater, pool_pledge, 
                #                                  filtered_seeds_pledge_3)
                #analyze.register("CoverageUpdaterModule", write_to_db_pledge_2, True)
                write_to_db_pledge_2 = promise_me(cc_module, filtered_seeds_pledge_3, cc_col)
                join_pledge = promise_me(join_module, write_to_db_pledge, write_to_db_pledge_2)
                analyze.register("Join", join_pledge, True)
                unlock_pledge = promise_me(UnLock(parameter_set_manager, query_pledge), join_pledge)
                analyze.register("UnLock", unlock_pledge, True)
                res.append(unlock_pledge)

            # drain all sources
            res.simultaneous_get(parameter_set_manager.get_num_threads())
            for inserter in inserter_vec:
                # @todo for some reason the destructor does not trigger automatically :(
                inserter.get().close(pool_pledge.get()) 

            cov_table.init_coverage_from_col(get_jump_inserter.cpp_module.id, cc_col.get())

            analyze.analyze(runtime_file)

        return get_jump_inserter.cpp_module.id

    db_conn = DbConn(dataset_name)
    JumpRunTable(db_conn)

    jump_table = SvJumpTable(db_conn)
    jump_table.drop_indices(0) # number does nothing at the moment

    cov_table = CoverageTable(db_conn)
    cov_table.drop_indices()

    jump_id = scope()

    analyze = AnalyzeRuntimes()

    print("num jumps:", jump_table.num_jumps(jump_id))

    print("creating index...")
    start = datetime.datetime.now()
    jump_table.create_indices( jump_id )
    cov_table.create_indices()
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("create_indices", delta.total_seconds(), False, lambda x: x)
    print("created index")
    analyze.analyze(runtime_file)

    if not runtime_file is None:
        runtime_file.write("sv_jump_run_id is " + str(jump_id) + "\n")

    # return the run_id
    return jump_id

