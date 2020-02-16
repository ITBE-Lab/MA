from .aligner import *
from .analyzeRuntimes import *

def to_file_path_vec(string_vec):
    return libMA.filePathVector([libMA.path(x) for x in string_vec])

def insert_paired_reads(parameter_set, dataset_name, sequencer_name, filename_vec1, filename_vec2, runtime_file=None):
    file_reader = PairedFileReader(parameter_set, to_file_path_vec(filename_vec1), to_file_path_vec(filename_vec2))
    get_read_inserter = GetPairedReadInserter(parameter_set, DbConn(dataset_name), sequencer_name)
    read_inserter_module = PairedReadInserterModule(parameter_set)
    module_lock = Lock(parameter_set)
    module_get_first = GetFirstQuery(parameter_set)
    module_get_second = GetSecondQuery(parameter_set)

    pool_pledge = Pledge()
    pool_pledge.set(PoolContainer(parameter_set.get_num_threads() + 1, dataset_name))

    analyze = AnalyzeRuntimes()

    res = VectorPledge()
    inserter_vec = []
    queries_pledge = promise_me(file_reader)
    analyze.register("PairedFileReader", queries_pledge, False)
    for _ in range(parameter_set.get_num_threads()):
        locked_query_pair = promise_me(module_lock, queries_pledge)
        locked_query_primary = promise_me(module_get_first, locked_query_pair)
        locked_query_mate = promise_me(module_get_second, locked_query_pair)
        read_inserter = promise_me(get_read_inserter, pool_pledge)
        inserter_vec.append(read_inserter)
        # insert the queries
        empty = promise_me(read_inserter_module, read_inserter, pool_pledge, locked_query_primary, locked_query_mate)
        analyze.register("PairedReadInserterModule", empty, True)
        # unlock the locked query
        unlock = promise_me(UnLock(parameter_set, locked_query_pair), empty)
        res.append(unlock)

    res.simultaneous_get(parameter_set.get_num_threads())

    for inserter in inserter_vec:
        inserter.get().close() # @todo for some reason the destructor does not trigger automatically :(

    analyze.analyze(runtime_file)

    return get_read_inserter.cpp_module.id

def insert_reads_from_pledge(parameter_set, dataset_name, sequencer_name, queries_pledge, runtime_file=None):
    get_read_inserter = GetReadInserter(parameter_set, DbConn(dataset_name), sequencer_name)
    read_inserter_module = ReadInserterModule(parameter_set)
    module_lock = Lock(parameter_set)

    pool_pledge = Pledge()
    pool_pledge.set(PoolContainer(parameter_set.get_num_threads() + 1, dataset_name))

    analyze = AnalyzeRuntimes()

    res = VectorPledge()
    inserter_vec = []
    analyze.register("GetQueries", queries_pledge, False)
    for _ in range(parameter_set.get_num_threads()):
        locked_query = promise_me(module_lock, queries_pledge)
        read_inserter = promise_me(get_read_inserter, pool_pledge)
        inserter_vec.append(read_inserter)
        # insert the queries
        empty = promise_me(read_inserter_module, read_inserter, pool_pledge, locked_query)
        analyze.register("ReadInserterModule", empty, True)
        unlock = promise_me(UnLock(parameter_set, locked_query), empty)
        res.append(unlock)

    res.simultaneous_get(parameter_set.get_num_threads())

    for inserter in inserter_vec:
        inserter.get().close() # @todo for some reason the destructor does not trigger automatically :(

    analyze.analyze(runtime_file)

    return get_read_inserter.cpp_module.id

def insert_reads(parameter_set, dataset_name, sequencer_name, filename_vec1, runtime_file=None):
    file_reader = FileReader(parameter_set, to_file_path_vec(filename_vec1))
    return insert_reads_from_pledge(parameter_set, dataset_name, sequencer_name, promise_me(file_reader), runtime_file)

def insert_reads_vec(parameter_set, dataset_name, sequencer_name, query_vec, runtime_file=None):
    splitter = StaticNucSeqSplitter(parameter_set, query_vec)
    return insert_reads_from_pledge(parameter_set, dataset_name, sequencer_name,
                                    promise_me(splitter), runtime_file)
