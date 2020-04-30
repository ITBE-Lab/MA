from MSV import *
from MA import *
import tempfile

#
# there is a ppt file that describes the example used in this test
#

def insert_calls(db_conn, dataset_name):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "simulated_sv",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    # deletion
    sv_inserter.insert(SvCall(4, 7, 0, 0, False, 1000))  # a

    # inversion
    sv_inserter.insert(SvCall(9, 14, 0, 0, True, 1000))  # b
    sv_inserter.insert(SvCall(15, 10, 0, 0, True, 1000)) # c

    # insertion
    insertion = SvCall(16, 17, 0, 0, False, 1000)
    insertion.inserted_sequence = NucSeq("TGTT")
    sv_inserter.insert(insertion) # d

    # translocation
    sv_inserter.insert(SvCall(0, 19, 0, 0, False, 1000))  # e
    sv_inserter.insert(SvCall(19, 1, 0, 0, False, 1000))  # f
    sv_inserter.insert(SvCall(18, 20, 0, 0, False, 1000))  # g

    sv_inserter.close(pool)

    return get_inserter.cpp_module.id

def get_reference():
    reference = Pack()
    reference.append("chr1", "chr1-desc", NucSeq("GATCGTATC"))
    reference.append("chr2", "chr2-desc", NucSeq("CTCGTCAACAG"))
    return reference

if __name__ == "__main__":
    db_conn = DbConn({"SCHEMA": {"NAME": "tmp_2", "FLAGS": ["DROP_ON_CLOSURE"]}})

    run_id = insert_calls(db_conn, "tmp_2")
    reference = get_reference()

    seeds, inserts = SvCallTable(db_conn).calls_to_seeds(reference, run_id)

    seed_printer = SeedPrinter(ParameterSetManager(), "seed")
    #seed_printer.execute(seeds) # uncomment to draw the seeds with bokeh

    seeds_to_expect = [
        # q,r,l,forw
        (0,0,1,True), # M
        (1,19,1,True), # P
        (2,1,4,True), # U
        (6,7,2,True), # W - chr1
        (8,9,1,True), # W - chr2
        (9,15,5,False), # ~X
        (14,15,2,True), # Y
        (20,17,2,True), # Z
    ]
    assert len(seeds) == len(inserts)
    assert len(seeds) == len(seeds_to_expect)

    for seed, (q,r,l,f) in zip(seeds, seeds_to_expect):
        assert seed.start == q
        assert seed.start_ref == r
        assert seed.size == l
        assert seed.on_forward_strand == f

    for idx, insert in enumerate(inserts):
        if idx != 6:
            assert len(insert) == 0
        else:
            assert str(insert) == "TGTT"