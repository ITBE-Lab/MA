from MA import *
from .analyzeRuntimes import AnalyzeRuntimes
import datetime


def sweep_sv_jumps(parameter_set_manager, sv_db, run_id, name, desc, sequencer_ids, pack,
                       out_file=None):
    analyze = AnalyzeRuntimes()
    # creates scope so that deconstructor of call inserter is triggered (commits insert transaction)
    def graph():
        print("\tsetting up graph...")
        
        pack_pledge = Pledge()
        pack_pledge.set(pack)

        section_fac = libMA.GenomeSectionFactory(parameter_set_manager, pack)
        lock_module = Lock(parameter_set_manager)
        sweep2 = libMA.ExactCompleteBipartiteSubgraphSweep(parameter_set_manager, sv_db, pack, sequencer_ids[0])

        sv_caller_run_id = sv_db.insert_sv_caller_run(name, desc, run_id)
        filter1 = libMA.FilterLowSupportShortCalls(parameter_set_manager)
        filter2 = libMA.FilterFuzzyCalls(parameter_set_manager)
        filter5 = libMA.FilterDiagonalLineCalls(parameter_set_manager)
        filter6 = libMA.FilterLowScoreCalls(parameter_set_manager)
        call_ambiguity = libMA.ComputeCallAmbiguity(parameter_set_manager)
        assert len(sequencer_ids) == 1

        call_inserter = libMA.SvCallInserter(sv_db, sv_caller_run_id)

        res = VectorPledge()
        sections_pledge = promise_me(section_fac) # @note this cannot be in the loop (synchronization!)
        sinks = []
        # graph for single reads
        for _ in range(1) : # Arne reduced to 1 from parameter_set_manager.get_num_threads()
            # in order to allow multithreading this module needs individual db connections for each thread
            sweep1 = libMA.CompleteBipartiteSubgraphSweep(parameter_set_manager, sv_db, pack, run_id, sequencer_ids[0])
            sink = libMA.BufferedSvCallSink(parameter_set_manager, call_inserter)
            filter3 = libMA.ConnectorPatternFilter(parameter_set_manager, sv_db)
            sinks.append(sink)

            section_pledge = promise_me(lock_module, sections_pledge)
            sweep1_pledge = promise_me(sweep1, section_pledge)
            #analyze.register("CompleteBipartiteSubgraphSweep", sweep1_pledge)
            analyze.register("CompleteBipartiteSubgraphSweep::init", sweep1, True, lambda x: x.cpp_module.time_init)
            analyze.register("CompleteBipartiteSubgraphSweep::outer_while", sweep1, True,
                             lambda x: x.cpp_module.time_complete_while - x.cpp_module.time_inner_while,)
            analyze.register("CompleteBipartiteSubgraphSweep::inner_while", sweep1, True,
                             lambda x: x.cpp_module.time_inner_while)
            sweep2_pledge = promise_me(sweep2, sweep1_pledge)
            analyze.register("ExactCompleteBipartiteSubgraphSweep", sweep2_pledge, True)
            #filters

            filter1_pledge = promise_me(filter1, sweep2_pledge)
            analyze.register("FilterLowSupportShortCalls", filter1_pledge, True)
            filter2_pledge = promise_me(filter2, filter1_pledge)
            analyze.register("FilterFuzzyCalls", filter2_pledge, True)

            #filter3_pledge = promise_me(filter3, filter2_pledge, pack_pledge) # this filter was off already
            #analyze.register("[4] ConnectorPatternFilter", filter3_pledge)

            filter3_pledge = promise_me(filter5, filter2_pledge)
            analyze.register("FilterDiagonalLineCalls", filter3_pledge, True)

            call_ambiguity_pledge = promise_me(call_ambiguity, filter3_pledge, pack_pledge)
            analyze.register("ComputeCallAmbiguity", call_ambiguity_pledge, True)

            filter6_pledge = promise_me(filter6, call_ambiguity_pledge)
            analyze.register("FilterLowScoreCalls", filter6_pledge, True)

            write_to_db_pledge = promise_me(sink, filter6_pledge)
            analyze.register("SvCallSink", write_to_db_pledge, True)
            unlock_pledge = promise_me(UnLock(parameter_set_manager, section_pledge), write_to_db_pledge)
            res.append(unlock_pledge)

        # drain all sources
        print("\texecuting graph...")
        res.simultaneous_get( parameter_set_manager.get_num_threads() )
        print("\tcommiting calls...")
        for sink in sinks:
            sink.cpp_module.commit( True )
        print("\tdone")

        call_inserter.end_transaction()

        return sv_caller_run_id

    sv_caller_run_id = graph()
    print("done sweeping")
    print("num calls:", sv_db.get_num_calls(sv_caller_run_id, 0))

    print("overlapping...")
    start = datetime.datetime.now()
    num_combined = libMA.combine_overlapping_calls(parameter_set_manager, sv_db, sv_caller_run_id)
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("combine_overlapping_calls", delta.total_seconds(), False, lambda x: x)
    print("done overlapping; combined", num_combined, "calls")

    print("computing score index...")
    start = datetime.datetime.now()
    sv_db.add_score_index(sv_caller_run_id)
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("compute_score_index", delta.total_seconds(), False, lambda x: x)
    print("done computing score index")

    analyze.analyze(out_file)
    if not out_file is None:
        out_file.write("run_id is " + str(sv_caller_run_id) + "\n")