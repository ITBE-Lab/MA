from MS import *
from MA import *
from MSV import *


def insert_reads(parameter_set, dataset_name, sequencer_name, file_queue, file_queue_2=None,
                 runtime_file=None, coverage=50):
    #parameter_set.by_name("Number of Threads").set(16)
    #parameter_set.by_name("Use all Processor Cores").set(False)
    #assert parameter_set.get_num_threads() == 16

    combined_queue = file_queue
    lock = Lock(parameter_set)
    single_con = DbConn(dataset_name)
    if not file_queue_2 is None:
        combined_queue = combine_file_streams(file_queue, file_queue_2)
        module_get_first = GetFirstQuery(parameter_set)
        module_get_second = GetSecondQuery(parameter_set)

        #set up modules for paired read insert
        queue_picker = PairedFilePicker(parameter_set)
        queue_placer = PairedFileNucSeqPlacer(parameter_set)
        file_reader = PairedFileReader(parameter_set)
        get_read_inserter = GetPairedReadInserter(parameter_set, single_con, sequencer_name)
        read_inserter_module = PairedReadInserterModule(parameter_set)
        printer = ProgressPrinterPairedFileStreamQueue(parameter_set)
    else:
        #set up modules for single read insert
        queue_picker = FilePicker(parameter_set)
        queue_placer = FileNucSeqPlacer(parameter_set)
        file_reader = FileReader(parameter_set)
        get_read_inserter = GetReadInserter(parameter_set, single_con, sequencer_name)
        read_inserter_module = ReadInserterModule(parameter_set)
        printer = ProgressPrinterFileStreamQueue(parameter_set)

    hash_counter_container = HashCounter()
    hash_counter = Pledge()
    hash_counter.set(hash_counter_container)
    mm_counter_module = MMCounterModule(parameter_set)

    queue_pledge = Pledge()
    queue_pledge.set(combined_queue)

    pool_pledge = Pledge()
    pool_pledge.set(PoolContainer(parameter_set.get_num_threads() + 1, dataset_name))

    analyze = AnalyzeRuntimes()

    # put together the graph
    res = VectorPledge()
    inserter_vec = []
    for _ in parallel_graph(parameter_set.get_num_threads()):
        picked_file = promise_me(queue_picker, queue_pledge)
        analyze.register("queue_picker", picked_file, True)
        locked_file = promise_me(lock, picked_file)
        analyze.register("lock", locked_file, True)

        query_ = promise_me(file_reader, locked_file)
        analyze.register("file_reader", query_, True)

        query = promise_me(queue_placer, query_, locked_file, queue_pledge)
        analyze.register("queue_placer", query, True)

        read_inserter = promise_me(get_read_inserter, pool_pledge)
        analyze.register("get_read_inserter", read_inserter, True)
        inserter_vec.append(read_inserter)

        if not file_queue_2 is None:
            query_primary = promise_me(module_get_first, query)
            analyze.register("get_first", query_primary, True)

            counted_primary = promise_me(mm_counter_module, query_primary, hash_counter)
            analyze.register("mm_counter", counted_primary , True)

            query_mate = promise_me(module_get_second, query)
            analyze.register("get_second", query_mate, True)

            counted_mate = promise_me(mm_counter_module, query_mate, hash_counter)
            analyze.register("mm_counter", counted_mate , True)

            empty = promise_me(read_inserter_module, read_inserter, pool_pledge, counted_primary, counted_mate)
            analyze.register("read_inserter", empty, True)

        else:
            counted_query = promise_me(mm_counter_module, query, hash_counter)
            analyze.register("mm_counter", counted_query, True)

            empty = promise_me(read_inserter_module, read_inserter, pool_pledge, counted_query)
            analyze.register("read_inserter", empty, True)

        empty_print = promise_me(printer, empty, queue_pledge)
        analyze.register("printer", empty_print, True)

        unlock = promise_me(UnLock(parameter_set, locked_file), empty_print)
        res.append(unlock)

    # run the graph
    res.simultaneous_get(parameter_set.get_num_threads())

    for inserter in inserter_vec:
        inserter.get().close(pool_pledge.get()) # @todo for some reason the destructor does not trigger automatically :(

    HashFilterTable(single_con).insert_counter_set(get_read_inserter.cpp_module.id, hash_counter_container, coverage)

    analyze.analyze(runtime_file)

    return get_read_inserter.cpp_module.id

def insert_reads_path_string_vec(parameter_set, dataset_name, sequencer_name, string_vec, string_vec_2=None,
                                 runtime_file=None, coverage=50):
    def to_file_queue(string_vec):
        if string_vec is None:
            return None
        file_queue = FileQueue()
        for string in string_vec:
            file_queue.add(FileStreamFromPath(string))
        return file_queue

    return insert_reads(parameter_set, dataset_name, sequencer_name, to_file_queue(string_vec),
                        to_file_queue(string_vec_2), runtime_file, coverage)

def write_used_reads_to_file(dataset_name, file_name):
    reads = ReadTable(DbConn(dataset_name)).get_used_reads()
    with open(file_name + ".fasta", "w") as fasta_out:
        for read in reads:
            fasta_out.write(">")
            fasta_out.write(read.name)
            fasta_out.write("\n")
            fasta_out.write(str(read))
            fasta_out.write("\n")

def iterate_reads(parameter_set, db_name, seq_id):
    db_conn_pool = PoolContainer(1, db_name)
    nuc_seq_fetcher = NucSeqFetcher(parameter_set)
    nuc_seq_query_getter = GetNucSeqFromSqlQuery(parameter_set, seq_id, 0, 1, True, True)
    nuc_seq_query = nuc_seq_query_getter.execute(db_conn_pool)
    while True:
        read = nuc_seq_fetcher.execute(db_conn_pool, nuc_seq_query)
        if read is None:
            break
        yield read
