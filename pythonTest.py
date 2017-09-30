from LAuS import *
from alignmentPrinter import AlignmentPrinter
import random

run()
exit()

q = ""

#random.seed(1)

def create_and_store():
    seq = ""
    for _ in range(1000000):
        char = random.randint(1,4)
        if char == 1:
            seq += "a"
        elif char == 2:
            seq += "c"
        elif char == 3:
            seq += "t"
        else:
            seq += "g"

            
    ref = NucSeq(seq)

    rev_ref = NucSeq(seq)
    rev_ref.reverse()

    query = NucSeq(q)

    fm_index = FMIndex(ref)
    fm_index.store("test")
    rev_fm_index = FMIndex(rev_ref)
    rev_fm_index.store("rev_test")

    ref_seq = BWAPack()
    ref_seq.append("name", "no comment",ref)
    ref_seq.store("test_pack")


ref_seq = BWAPack()
ref_seq.load("test_pack")
fm_index = FMIndex()
fm_index.load("test")
rev_fm_index = FMIndex()
rev_fm_index.load("rev_test")

q_from = random.randint(0, 1000)
q_to = q_from + 1000
q = str(ref_seq.extract_from_to(q_from, q_to))
for _ in range(10):
    pos = random.randint(1,len(q)-1)
    char = random.randint(1,4)
    if char == 1:
        q = q[:pos-1] + "a" + q[pos:]
    elif char == 2:
        q = q[:pos-1] + "c" + q[pos:]
    elif char == 3:
        q = q[:pos-1] + "t" + q[pos:]
    else:
        q = q[:pos-1] + "g" + q[pos:]

for _ in range(10):
    pos = random.randint(1,len(q)-11)
    l = random.randint(1,10)
    q = q[:pos-1] + q[pos + l:]

for _ in range(10):
    pos = random.randint(1,len(q)-1)
    l = random.randint(1,10)
    for _ in range(l):
        char = random.randint(1,4)
        if char == 1:
            q = q[:pos] + "a" + q[pos:]
        elif char == 2:
            q = q[:pos] + "c" + q[pos:]
        elif char == 3:
            q = q[:pos] + "t" + q[pos:]
        else:
            q = q[:pos] + "g" + q[pos:]



query = NucSeq(q)

seg = Segmentation(True, True, 10, 100000)
seg.bSkipLongBWTIntervals = False

segments = seg.execute((fm_index, rev_fm_index, query, ref_seq))

seeds = Seeds()

for segment in segments:
    start = segment.start()
    end = segment.end()
    sequence = ""
    for i in range(start, end):
        sequence = sequence + query[i]
    print("segment: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence)
    seeds.append(segment.get_seeds(fm_index, rev_fm_index))

#===================== unused code from here on =====================

anc = NlongestIntervalsAsAnchors(2)
anchors = anc.execute((segments,))
    
for anchor in anchors:
    start = anchor.start()
    end = anchor.end()
    sequence = ""
    for i in range(start, end):
        sequence = sequence + query[i]
    print("anchor: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence)

bucketing = Bucketing()
bucketing.strip_size = 50
strips_of_consideration = bucketing.execute((
    segments,
    anchors,
    query,
    ref_seq,
    fm_index,
    rev_fm_index))

if len(strips_of_consideration) == 0:
    print("no match found")
    exit()

print("found " + str(len(strips_of_consideration)) + " strips of consideration")
print("scores before linesweep:")
for strip in strips_of_consideration:
    print(strip.get_score())

best_strip = []
liesweep = LineSweep()
for strip in strips_of_consideration:
    best_strip.append(liesweep.execute((query, ref_seq, strip)))


print("scores after linesweep:")
for strip in best_strip:
    print(strip.get_score())

best = 0
for index, strip in enumerate(best_strip):
    if strip.get_score() > best_strip[best].get_score():
        best = index

print("best score: " + str(best_strip[best].get_score()))
nmw = NeedlemanWunsch()

align = nmw.execute((best_strip[best], query, ref_seq))

printer = AlignmentPrinter()

printer.execute((align, query, ref_seq))

