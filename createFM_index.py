from LAuS import *
import random

seq = ""
q = ""

for _ in range(100000):
    char = random.randint(1,4)
    if char == 1:
        seq += "a"
    elif char == 2:
        seq += "c"
    elif char == 3:
        seq += "t"
    else:
        seq += "g"

for _ in range(100):
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

fm_index = FMIndex(ref)
rev_fm_index = FMIndex(rev_ref)

ref_seq = BWAPack()
ref_seq.append("name", "no comment",ref)

seg = Segmentation(True, True, 10, 100000)
seg.bSkipLongBWTIntervals = False

segments = seg.execute((fm_index, rev_fm_index, query, ref_seq))

iterator = segments.begin()
while iterator.exits():
    start = iterator.get().start()
    end = iterator.get().end()
    sequence = ""
    for i in range(start, end):
        sequence = sequence + query.at(i)
    print "segment: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence
    iterator.next()



anc = NlongestIntervalsAsAnchors(2)
anchors = anc.execute((segments,))
    
iterator = anchors.begin()
while iterator.exits():
    start = iterator.get().start()
    end = iterator.get().end()
    sequence = ""
    for i in range(start, end):
        sequence = sequence + query.at(i)
    print "anchor: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence
    iterator.next()

bucketing = Bucketing()
strips_of_consideration = bucketing.execute((
    segments,
    anchors,
    query,
    ref_seq,
    fm_index,
    rev_fm_index)).x

print "found " + str(len(strips_of_consideration)) + " strips of consideration with the scores:"
for strip in strips_of_consideration:
    print strip.get_score()
