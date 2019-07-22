from MA import *
import tempfile

db_name = tempfile.gettempdir() + "/.tmp_2.db"

print(db_name)

database = SV_DB(db_name, "create")
def insert(database):
    sv_inserter = SvCallInserter(database, "simulated sv", "the sv's that were simulated", -1)

    # insertion
    insertion = SvCall(10, 11, 1, 1, False, 1000)
    insertion.inserted_sequence = NucSeq("GGG")
    sv_inserter.insert_call(insertion)                        # 1

    # deletion 2-4
    sv_inserter.insert_call(SvCall(11, 14, 1, 1, False, 1000))  # 2

    # translocation between chromosomes
    sv_inserter.insert_call(SvCall(14, 26, 1, 1, False, 1000)) # 3
    sv_inserter.insert_call(SvCall(28, 18, 1, 1, False, 1000)) # 4

    # inversion
    sv_inserter.insert_call(SvCall(19, 22, 1, 1, True, 1000))  # 5
    sv_inserter.insert_call(SvCall(23, 20, 1, 1, True, 1000)) # 6

    # translocation (continued)
    sv_inserter.insert_call(SvCall(25, 15, 1, 1, False, 1000)) # 7
    sv_inserter.insert_call(SvCall(17, 29, 1, 1, False, 1000)) # 8

    return sv_inserter.sv_caller_run_id

run_id = insert(database)

reference = Pack()
reference.append("chr1", "chr1-desc", NucSeq("AGCTC"))
reference.append("chr2", "chr2-desc", NucSeq("GCT"))
reference.append("chr3", "chr3-desc", NucSeq("CA"))
reference.append("chr4", "chr4-desc", NucSeq("ACGTACGNACGGCAT"))
reference.append("chr5", "chr5-desc", NucSeq("GCGCG"))
reconstr = database.reconstruct_sequenced_genome(reference, run_id)

print(str(reconstr.extract_forward_strand_n())," ?= AGCTCGCTCAAGGGCACGCACGCCATGCGNG")
print("contig lengths:")
for x in reconstr.contigLengths():
    print(x)

assert str(reconstr.extract_forward_strand_n()) == "AGCTCGCTCAAGGGCACGCACGCCATGCGNG"
assert reconstr.contigLengths()[0] == 5
assert reconstr.contigLengths()[1] == 3
assert reconstr.contigLengths()[2] == 2
assert reconstr.contigLengths()[3] == 16
assert reconstr.contigLengths()[4] == 5

exit(0)