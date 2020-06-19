from MSV import *
from MS import *
from MA import *
import datetime

def compute_sv_jumps(parameter_set_manager, fm_index, pack, dataset_name, seq_ids=0, runtime_file=None):
    #parameter_set_manager.by_name("Number of Threads").set(1)
    #parameter_set_manager.by_name("Use all Processor Cores").set(False)
    #assert parameter_set_manager.get_num_threads() == 1

    def scope():
        nuc_seq_fetcher = NucSeqFetcher(parameter_set_manager)
        lock_module = Lock(parameter_set_manager)
        seeding_module = BinarySeeding(parameter_set_manager)
        jumps_from_seeds = SvJumpsFromSeeds(parameter_set_manager, pack)
        get_jump_inserter = GetJumpInserter(parameter_set_manager, DbConn(dataset_name), "MS-SV",
                                            "python built comp graph")
        jump_inserter_module = JumpInserterModule(parameter_set_manager)

        fm_pledge = Pledge()
        fm_pledge.set(fm_index)
        pack_pledge = Pledge()
        pack_pledge.set(pack)
        pool_pledge = Pledge()
        pool_pledge.set(PoolContainer(parameter_set_manager.get_num_threads() + 1, dataset_name))

        if not isinstance(seq_ids, list):
            seq_ids = [seq_ids]

        jobs = parameter_set_manager.get_num_threads()
        for seq_id in seq_ids:
            analyze = AnalyzeRuntimes()
            res = VectorPledge()
            inserter_vec = []
            # graph for single reads
            for idx in range(jobs):
                # @todo there should be a set of modules distributing reads
                # (only a problem if threads finish at different times)...
                nuc_seq_query_getter = GetNucSeqFromSqlQuery(parameter_set_manager, seq_id, idx, jobs, True, True)
                nuc_seq_query = promise_me(nuc_seq_query_getter, pool_pledge)
                analyze.register("GetNucSeqFromSqlQuery", nuc_seq_query, True)
                queries_pledge = promise_me(nuc_seq_fetcher, pool_pledge, nuc_seq_query)
                analyze.register("NucSeqFetcher", queries_pledge, True)
                query_pledge = promise_me(lock_module, queries_pledge)
                analyze.register("Lock", query_pledge, True)
                segments_pledge = promise_me(seeding_module, fm_pledge, query_pledge)
                analyze.register("BinarySeeding", segments_pledge, True)
                jumps_pledge = promise_me(jumps_from_seeds, segments_pledge, pack_pledge, fm_pledge, query_pledge)
                analyze.register("SvJumpsFromSeeds", jumps_pledge, True)

                jump_inserter = promise_me(get_jump_inserter, pool_pledge)
                inserter_vec.append(jump_inserter)
                analyze.register("GetJumpInserter", jump_inserter, True)
                write_to_db_pledge = promise_me(jump_inserter_module, jump_inserter, pool_pledge, jumps_pledge,
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

    jump_table = SvJumpTable(DbConn(dataset_name))
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

