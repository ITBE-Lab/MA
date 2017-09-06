from LAuS import *
import random

seq = ""
q = ""

for _ in range(1000):
    char = random.randint(1,4)
    if char == 1:
        seq += "a"
    elif char == 2:
        seq += "c"
    elif char == 3:
        seq += "t"
    else:
        seq += "g"

for _ in range(10):
    char = random.randint(1,4)
    if char == 1:
        q += "a"
    elif char == 2:
        q += "c"
    elif char == 3:
        q += "t"
    else:
        q += "g"


ref = NucSeq(seq)

rev_ref = NucSeq(seq)
rev_ref.reverse()

query = NucSeq(q)

fm_index = FM_index(ref)
rev_fm_index = FM_index(rev_ref)

ref_seq = BWAPack()
ref_seq.append("name", "no comment",ref)

seg = Segmentation(True, True, 10, 100000)
seg.bSkipLongBWTIntervals = False

segments = seg.execute([fm_index, rev_fm_index, query, ref_seq])

print("done")