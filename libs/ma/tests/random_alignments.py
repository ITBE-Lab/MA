from MA import *
import random

def random_nuc_seq(l):
    ret = ""
    for _ in range(l):
        ret += random.choice(['a', 'c', 'g', 't'])
    return ret

pack = Pack()
pack.append("chr1", "chr1-desc", NucSeq(random_nuc_seq(65536)))
fm_index = FMIndex(pack)

file_string = ""
for idx in range(1000):
    file_string += ">short_query"
    file_string += str(idx)
    file_string += "\n"
    file_string += random_nuc_seq(1000)
    file_string += "\n"
for idx in range(100):
    file_string += ">long_query"
    file_string += str(idx)
    file_string += "\n"
    file_string += random_nuc_seq(10000)
    file_string += "\n"

parameter_set = ParameterSetManager()
parameter_set.by_name("Detect Small Inversions").set(True) # actually run all modules

alignment_collector = AlignmentCollector(parameter_set)
def output_graph(alignment_pledge, pack_pledge, locked_query):
    return promise_me(alignment_collector,
                        locked_query, alignment_pledge, pack_pledge)

file_queue = FileQueue()
file_queue.add(StringStream(file_string))
quick_align(parameter_set, pack, fm_index, output_graph, file_queue)

exit(0)