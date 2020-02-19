from .aligner import *
from .analyzeRuntimes import *


def insert_reads(parameter_set, dataset_name, sequencer_name, file_queue, file_queue_2=None,
                 runtime_file=None):
    combined_queue = file_queue
    if not file_queue_2 is None:
        combined_queue = libMA.combine_file_streams(file_queue, file_queue_2)
        module_get_first = GetFirstQuery(parameter_set)
        module_get_second = GetSecondQuery(parameter_set)

        #set up modules for paired read insert
        queue_picker = libMA.PairedFilePicker(parameter_set)
        queue_placer = libMA.PairedFilePlacer(parameter_set)
        file_reader = PairedFileReader(parameter_set)
        get_read_inserter = GetPairedReadInserter(parameter_set, DbConn(dataset_name), sequencer_name)
        read_inserter_module = PairedReadInserterModule(parameter_set)
    else:
        #set up modules for single read insert
        queue_picker = libMA.FilePicker(parameter_set)
        queue_placer = libMA.FilePlacer(parameter_set)
        file_reader = FileReader(parameter_set)
        get_read_inserter = GetReadInserter(parameter_set, DbConn(dataset_name), sequencer_name)
        read_inserter_module = ReadInserterModule(parameter_set)

    queue_pledge = Pledge()
    queue_pledge.set(combined_queue)

    pool_pledge = Pledge()
    pool_pledge.set(PoolContainer(parameter_set.get_num_threads() + 1, dataset_name))

    analyze = AnalyzeRuntimes()

    # put together the graph
    res = VectorPledge()
    inserter_vec = []
    for _ in range(parameter_set.get_num_threads()):
        picked_file = promise_me(queue_picker, queue_pledge)
        analyze.register("queue_picker", picked_file, True)

        locked_file = promise_me(queue_placer, picked_file, queue_pledge)

        query = promise_me(file_reader, locked_file)
        analyze.register("file_reader", query, True)

        read_inserter = promise_me(get_read_inserter, pool_pledge)
        inserter_vec.append(read_inserter)

        if not file_queue_2 is None:
            query_primary = promise_me(module_get_first, query)
            query_mate = promise_me(module_get_second, query)

            empty = promise_me(read_inserter_module, read_inserter, pool_pledge, query_primary, query_mate)
        else:
            empty = promise_me(read_inserter_module, read_inserter, pool_pledge, query)
        analyze.register("read_inserter", empty, True)

        unlock = promise_me(UnLock(parameter_set, locked_file), empty)
        res.append(unlock)

    # run the graph
    res.simultaneous_get(parameter_set.get_num_threads())

    for inserter in inserter_vec:
        inserter.get().close(pool_pledge.get()) # @todo for some reason the destructor does not trigger automatically :(

    analyze.analyze(runtime_file)

    return get_read_inserter.cpp_module.id

def insert_reads_path_string_vec(parameter_set, dataset_name, sequencer_name, string_vec, string_vec_2=None,
                                 runtime_file=None):
    def to_file_queue(string_vec):
        if string_vec is None:
            return None
        file_queue = libMA.FileQueue()
        for string in string_vec:
            file_queue.add(libMA.FileStreamFromPath(string))
        return file_queue

    return insert_reads(parameter_set, dataset_name, sequencer_name, to_file_queue(string_vec),
                        to_file_queue(string_vec_2), runtime_file)
