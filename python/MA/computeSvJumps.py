from MA import *
from .analyzeRuntimes import AnalyzeRuntimes

def compute_sv_jumps(parameter_set_manager, fm_index, pack, sv_db, seq_id=0, runtime_file=None):
    #parameter_set_manager.by_name("Number of Threads").set(1)
    #parameter_set_manager.by_name("Use all Processor Cores").set(False)
    #assert parameter_set_manager.get_num_threads() == 1

    lock_module = Lock(parameter_set_manager)
    seeding_module = BinarySeeding(parameter_set_manager)
    run_id = sv_db.insert_sv_jump_run("MS-SV", "python built compt graph")
    jumps_from_seeds = SvJumpsFromSeeds(parameter_set_manager, seq_id, sv_db, pack)

    fm_pledge = Pledge()
    fm_pledge.set(fm_index)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    analyze = AnalyzeRuntimes()

    res = VectorPledge()
    jump_to_dbs = []
    # graph for single reads
    jobs = parameter_set_manager.get_num_threads() * 2
    for idx in range(jobs):
        # @todo there should be a set of modules distributing reads...
        nuc_seq_getter = AllNucSeqFromSql(parameter_set_manager, sv_db, seq_id, idx, jobs)
        jumps_to_db = libMA.BufferedSvDbInserter(parameter_set_manager, sv_db, run_id)
        jump_to_dbs.append(jumps_to_db)
        queries_pledge = promise_me(nuc_seq_getter)
        analyze.register("AllNucSeqFromSql", queries_pledge)
        query_pledge = promise_me(lock_module, queries_pledge)
        segments_pledge = promise_me(seeding_module, fm_pledge, query_pledge)
        analyze.register("BinarySeeding", segments_pledge)
        jumps_pledge = promise_me(jumps_from_seeds, segments_pledge, pack_pledge, fm_pledge, query_pledge)
        analyze.register("SvJumpsFromSeeds", jumps_pledge)
        write_to_db_pledge = promise_me(jumps_to_db, jumps_pledge, query_pledge)
        analyze.register("SvDbInserter", write_to_db_pledge)
        unlock_pledge = promise_me(UnLock(parameter_set_manager, query_pledge), write_to_db_pledge)
        res.append(unlock_pledge)

    # drain all sources
    res.simultaneous_get( parameter_set_manager.get_num_threads() )

    for jump_to_db in jump_to_dbs:
        jump_to_db.cpp_module.commit(False, True)

    analyze.analyze(runtime_file)
    if not runtime_file is None:
        runtime_file.write("sv_jump_run_id is " + str(run_id) + "\n")

    sv_db.create_jump_indices( run_id )

    # return the run_id
    return run_id

