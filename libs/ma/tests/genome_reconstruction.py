from MA import *
import tempfile

#
# there is a ppt file that describes the example used in this test
#

def insert_calls(database):
    sv_inserter = SvCallInserter(database, "simulated sv", "the sv's that were simulated", -1)

    # deletion
    sv_inserter.insert_call(SvCall(4, 7, 0, 0, False, 1000))  # a

    # inversion
    sv_inserter.insert_call(SvCall(9, 14, 0, 0, True, 1000))  # b
    sv_inserter.insert_call(SvCall(15, 10, 0, 0, True, 1000)) # c

    # insertion
    insertion = SvCall(16, 17, 0, 0, False, 1000)
    insertion.inserted_sequence = NucSeq("TGTT")
    sv_inserter.insert_call(insertion) # d

    # translocation
    sv_inserter.insert_call(SvCall(0, 19, 0, 0, False, 1000))  # e
    sv_inserter.insert_call(SvCall(19, 1, 0, 0, False, 1000))  # f
    sv_inserter.insert_call(SvCall(18, 20, 0, 0, False, 1000))  # g

    sv_inserter.end_transaction()
    return sv_inserter.sv_caller_run_id

def get_reference():
    reference = Pack()
    reference.append("chr1", "chr1-desc", NucSeq("GATCGTATC"))
    reference.append("chr2", "chr2-desc", NucSeq("CTCGTCAACAG"))
    return reference

if __name__ == "__main__":
    db_name = tempfile.gettempdir() + "/.tmp_2.db"

    print(db_name)

    database = SV_DB(db_name, "create")

    run_id = insert_calls(database)
    reference = get_reference()

    expected_sequence = "GGATCGTCCGACGAAATGTTCA"

    reconstr = database.reconstruct_sequenced_genome(reference, run_id)

    if str(reconstr.extract_forward_strand_n()) != expected_sequence:
        print("original sequence     ", reference.extract_forward_strand_n())
        print("expected sequence     ", expected_sequence)
        print("reconstructed sequence", reconstr.extract_forward_strand_n())
        for i, l in enumerate(reconstr.contigLengths()):
            print("contig", i,"length =", l)

    assert str(reconstr.extract_forward_strand_n()) == expected_sequence
    assert reconstr.contigLengths()[0] == 8
    assert reconstr.contigLengths()[1] == 14