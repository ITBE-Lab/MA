from MA import *
import tempfile
from genome_reconstruction import *

def insert_reads(database, reference, parameter_set):
    print("inserting calls...")
    ground_truth_id = insert_calls(database)
    print("done")

    print("creating fmd index...")
    fm_index = FMIndex(reference)
    print("done")

    print("reconstructing genome...")
    reconstr = database.reconstruct_sequenced_genome(reference, ground_truth_id)
    print("done")

    print("inserting reads...")
    reconstr_nuc_seq = reconstr.extract_forward_strand_n()
    reads = libMA.VecRetNuc()
    # create 10x coverage
    for x in range(5):
        reads.append(reconstr_nuc_seq)
    # also insert reads from rev comp strand
    reconstr_nuc_seq.complement()
    reconstr_nuc_seq.reverse()
    for x in range(5):
        reads.append(reconstr_nuc_seq)
    sequencer_id = libMA.insert_reads(reads, database, "perfect whole genome read", reference)
    print("done")

    return ground_truth_id, fm_index, sequencer_id


def fetch_call_rectangles(database, parameter_set, run_id):
    class CallRectangle:
        def __init__(self, call):
            self.from_start = call.from_start
            self.to_start = call.to_start
            self.from_size = call.from_size
            self.to_size = call.to_size

        # so that we can use jump in list to check wether each jump overlaps at least one ground truth call
        def __eq__(self, other):
            def interval_overlap(start_a, len_a, start_b, len_b):
                return start_a + len_a >= start_b and start_b + len_b >= start_a
            if isinstance(other, libMA.SvJump):
                return interval_overlap(self.from_start, self.from_size,
                                        other.from_start_same_strand(), other.from_size()) and \
                       interval_overlap(self.to_start, self.to_size, other.to_start(), other.to_size())
            return False
        
        def __str__(self):
            return str(("call", self.from_start, self.to_start, self.from_size, self.to_size))

    calls_from_db = SvCallsFromDb(parameter_set, database, run_id)
    ret = []
    while calls_from_db.hasNext():
        ret.append(CallRectangle(calls_from_db.next()))
    return ret

if __name__ == "__main__":
    parameter_set = ParameterSetManager()
    # @todo create example where there are only unique seeds...
    parameter_set.by_name("Minimal Seed Size SV").set(2)

    db_name = tempfile.gettempdir() + "/.tmp_3.db"
    database = SV_DB(db_name, "create")

    reference = get_reference()

    ground_truth_id, fm_index, sequencer_id = insert_reads(database, reference, parameter_set)
    del database

    print("computing jumps...")
    # @todo i should not have to reopen the DB but somehow it is locked if i don't
    database = SV_DB(db_name, "open")
    jump_id = compute_sv_jumps(parameter_set, fm_index, reference, database, sequencer_id)
    print("done")

    call_rectangles = fetch_call_rectangles(database, parameter_set, ground_truth_id)
    for x in call_rectangles:
        print(x)

    sweeper = SortedSvJumpFromSql(parameter_set, database, jump_id)
    while sweeper.has_next_start():
        jump = sweeper.get_next_start()
        print(("jump", jump.from_start_same_strand(), jump.to_start(), jump.from_size(), jump.to_size()))
        #assert jump in call_rectangles


