from LAuS import *
from setupaligner import set_up_aligner
import random
"""
before = "AATTTTATGACTTTTTATTTTCAG"
actual = "GTGTAAACGTCCTACTTTGGGGCAACTTGCCAGAGATAGAAGAGAGTACAGATGAAGATGTGTTAAATATCTCAGCAGAGGAGTGTATTAGATAA"
after = "ATGGAATTATGATATATATGATATACAA"

sequence = before + actual + after

start = len(before) - 20
end = len(before) + len(actual) + 20

for index in range(start, end):
    if sequence[index] == sequence[index+1] == 'G':
        print( "\"" + sequence[index-21:index-1] + "\"," )
    if sequence[index] == sequence[index+1] == 'C':
        s = sequence[index+2:index+23]
        s2 = ""
        for index in range(20):
            if s[index] == "G":
                s2 = "C" + s2
            if s[index] == "C":
                s2 = "G" + s2
            if s[index] == "A":
                s2 = "T" + s2
            if s[index] == "T":
                s2 = "A" + s2
        print( "\"" + s2 + "\"," )

exit()


analyse_crisper()
exit()
"""

q = ""

query_pledge = Pledge(ContainerType.nucSeq)
reference_pledge = Pledge(ContainerType.packedNucSeq)
fm_index_pledge = Pledge(ContainerType.fM_index)
rev_fm_index_pledge = Pledge(ContainerType.fM_index)

result_pledge = set_up_aligner(
    query_pledge,
    reference_pledge,
    fm_index_pledge,
    rev_fm_index_pledge
)

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

q_from = random.randint(0, 2000)
q_to = q_from + 2000
q = str(ref_seq.extract_from_to(q_from, q_to))
for _ in range(50):
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

for _ in range(50):
    pos = random.randint(1,len(q)-11)
    l = random.randint(1,10)
    q = q[:pos-1] + q[pos + l:]

for _ in range(50):
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

query_pledge.set(query)
reference_pledge.set(ref_seq)
fm_index_pledge.set(fm_index)
rev_fm_index_pledge.set(rev_fm_index)

result_pledge.next()
