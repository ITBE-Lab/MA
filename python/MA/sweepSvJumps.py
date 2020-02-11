from MA import *
from .analyzeRuntimes import AnalyzeRuntimes
import datetime


def sweep_sv_jumps(parameter_set_manager, dataset_name, run_id, name, desc, sequencer_ids, pack,
                       out_file=None):
    parameter_set_manager.by_name("Number of Threads").set(2)
    parameter_set_manager.by_name("Use all Processor Cores").set(False)
    assert parameter_set_manager.get_num_threads() == 2

    # creates scope so that deconstructor of call inserter is triggered (commits insert transaction)
    def graph():
        analyze = AnalyzeRuntimes()
        print("\tsetting up graph...")
        
        pack_pledge = Pledge()
        pack_pledge.set(pack)
        pool_pledge = Pledge()
        pool_pledge.set(PoolContainer(parameter_set_manager.get_num_threads() + 1, dataset_name))

        assert len(sequencer_ids) == 1

        section_fac = libMA.GenomeSectionFactory(parameter_set_manager, pack)
        lock_module = Lock(parameter_set_manager)
        sweep1 = libMA.CompleteBipartiteSubgraphSweep(parameter_set_manager, sequencer_ids[0])
        sweep2 = libMA.ExactCompleteBipartiteSubgraphSweep(parameter_set_manager)

        filter1 = libMA.FilterLowSupportShortCalls(parameter_set_manager)
        filter2 = libMA.FilterFuzzyCalls(parameter_set_manager)
        filter5 = libMA.FilterDiagonalLineCalls(parameter_set_manager)
        filter6 = libMA.FilterLowScoreCalls(parameter_set_manager)
        call_ambiguity = libMA.ComputeCallAmbiguity(parameter_set_manager)

        get_call_inserter = GetCallVectorInserter(parameter_set_manager, DbConn(dataset_name), name, desc, run_id)
        call_inserter_module = CallVectorInserterModule(parameter_set_manager)

        res = VectorPledge()
        sections_pledge = promise_me(section_fac) # @note this cannot be in the loop (synchronization!)
        # graph for single reads
        for _ in range(parameter_set_manager.get_num_threads()):
            section_pledge = promise_me(lock_module, sections_pledge)

            sweep1_pledge = promise_me(sweep1, pool_pledge, section_pledge, pack_pledge)
            #analyze.register("CompleteBipartiteSubgraphSweep", sweep1_pledge) this would cause duplicate time
            analyze.register("CompleteBipartiteSubgraphSweep::init", sweep1, True, lambda x: x.cpp_module.time_init)
            analyze.register("CompleteBipartiteSubgraphSweep::outer_while", sweep1, True,
                             lambda x: x.cpp_module.time_complete_while - x.cpp_module.time_inner_while)
            analyze.register("CompleteBipartiteSubgraphSweep::inner_while", sweep1, True,
                             lambda x: x.cpp_module.time_inner_while)
            sweep2_pledge = promise_me(sweep2, sweep1_pledge, pack_pledge)
            analyze.register("ExactCompleteBipartiteSubgraphSweep", sweep2_pledge, True)

            #### FILTERS ####

            filter1_pledge = promise_me(filter1, sweep2_pledge)
            analyze.register("FilterLowSupportShortCalls", filter1_pledge, True)
            filter2_pledge = promise_me(filter2, filter1_pledge)
            analyze.register("FilterFuzzyCalls", filter2_pledge, True)

            #filter3 = libMA.ConnectorPatternFilter(parameter_set_manager, sv_db)
            #filter3_pledge = promise_me(filter3, filter2_pledge, pack_pledge) # this filter was off already
            #analyze.register("[4] ConnectorPatternFilter", filter3_pledge)

            filter3_pledge = promise_me(filter5, filter2_pledge)
            analyze.register("FilterDiagonalLineCalls", filter3_pledge, True)

            call_ambiguity_pledge = promise_me(call_ambiguity, filter3_pledge, pack_pledge)
            analyze.register("ComputeCallAmbiguity", call_ambiguity_pledge, True)

            filter6_pledge = promise_me(filter6, call_ambiguity_pledge)
            analyze.register("FilterLowScoreCalls", filter6_pledge, True)

            call_inserter = promise_me(get_call_inserter, pool_pledge)

            write_to_db_pledge = promise_me(call_inserter_module, call_inserter, pool_pledge, filter6_pledge)
            analyze.register("CallInserterModule", write_to_db_pledge, True)

            unlock_pledge = promise_me(UnLock(parameter_set_manager, section_pledge), write_to_db_pledge)
            analyze.register("UnLock", unlock_pledge, True)
            res.append(unlock_pledge)

        # drain all sources
        print("\texecuting graph...")
        res.simultaneous_get( parameter_set_manager.get_num_threads() )
        print("\tdone executing")
        analyze.analyze(out_file)

        return get_call_inserter.cpp_module.id

    sv_caller_run_id = graph()
    print("done sweeping")
    analyze = AnalyzeRuntimes()

    conn = DbConn(dataset_name)
    call_table = SvCallTable(conn)

    print("num calls:", call_table.num_calls(sv_caller_run_id, 0))

    print("overlapping...")
    start = datetime.datetime.now()
    num_combined = libMA.combine_overlapping_calls(parameter_set_manager, conn, sv_caller_run_id)
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("combine_overlapping_calls", delta.total_seconds(), False, lambda x: x)
    print("done overlapping; combined", num_combined, "calls")

    print("computing score index...")
    start = datetime.datetime.now()
    call_table.add_score_index(sv_caller_run_id) # does nothing at the moment
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("compute_score_index", delta.total_seconds(), False, lambda x: x)
    print("done computing score index")

    analyze.analyze(out_file)
    if not out_file is None:
        out_file.write("run_id is " + str(sv_caller_run_id) + "\n")