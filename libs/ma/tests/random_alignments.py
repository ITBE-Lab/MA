from MA import *
import random

def random_nuc_seq(l):
    ret = ""
    for _ in range(l):
        ret += random.choice(['a', 'c', 'g', 't'])
    return NucSeq(ret)

pack = Pack()
pack.append("chr1", "chr1-desc", random_nuc_seq(65536))
fm_index = FMIndex(pack)

query_vec = QueryVector()
for _ in range(1000):
    query_vec.append(random_nuc_seq(1000))
for _ in range(100):
    query_vec.append(random_nuc_seq(10000))

parameter_set = ParameterSetManager()
parameter_set.by_name("Detect Small Inversions").set(True) # actually run all modules

splitter = NucSeqSplitter(parameter_set)
query_vec_pledge = Pledge()
query_vec_pledge.set(query_vec)
queries_pledge = promise_me(splitter, query_vec_pledge)

alignment_collector = AlignmentCollector(parameter_set)

quick_align(parameter_set, pack, fm_index, queries_pledge, alignment_collector)

exit(0)