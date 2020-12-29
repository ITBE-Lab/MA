from MSV import *
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
    sv_inserter.insert(SvCall(4, 7, 0, 0, True, True, 1000))  # a

    # inversion
    sv_inserter.insert(SvCall(9, 14, 0, 0, True, False, 1000))  # b
    sv_inserter.insert(SvCall(10, 15, 0, 0, False, True, 1000)) # c

    # insertion
    insertion = SvCall(16, 17, 0, 0, True, True, 1000)
    insertion.inserted_sequence = NucSeq("TGTT")
    sv_inserter.insert(insertion) # d

    # translocation
    sv_inserter.insert(SvCall(0, 19, 0, 0, True, True, 1000))  # e
    sv_inserter.insert(SvCall(1, 19, 0, 0, False, False, 1000))  # f
    sv_inserter.insert(SvCall(18, 20, 0, 0, True, True, 1000))  # g

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

    expected_sequence = "GGATCGTCCGACGAAATGTTCA"

    reconstr = reconstruct_sequenced_genome(reference, run_id)

    if str(reconstr.extract_forward_strand_n()) != expected_sequence:
        print("original sequence     ", reference.extract_forward_strand_n())
        print("expected sequence     ", expected_sequence)
        print("reconstructed sequence", reconstr.extract_forward_strand_n())
        for i, l in enumerate(reconstr.contigLengths()):
            print("contig", i,"length =", l)
    assert str(reconstr.extract_forward_strand_n()) == expected_sequence

    if [*reconstr.contigLengths()] != [8, 14]:
        print("contig lengths don't match; expects 8 14, got ", *reconstr.contigLengths())
    assert [*reconstr.contigLengths()] == [8, 14]