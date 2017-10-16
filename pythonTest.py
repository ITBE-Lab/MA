from LAuS import *
import random
import gc


def mutate(char):
    num = 0
    if char == "c":
        num = 1
    elif char == "t":
        num = 2
    else:
        num = 3
    num += random.randint(1,3)
    num %= 4
    if num == 0:
        return 'a'
    elif num == 2:
        return 'c'
    elif num == 3:
        return 't'
    else:
        return 'g'

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

query_pledge = [Pledge(ContainerType.nucSeq)] * 10
reference_pledge = Pledge(ContainerType.packedNucSeq)
fm_index_pledge = Pledge(ContainerType.fM_index)

result_pledges = set_up_aligner(
    query_pledge,
    reference_pledge,
    fm_index_pledge
)



ref_seq = Pack()
#ref_seq.append("no name", "no desc", NucSeq("AACG"))
ref_seq.load("/mnt/ssd0/chrom/random/pack")
#fm_index = FMIndex(ref_seq)
fm_index = FMIndex()
fm_index.load("/mnt/ssd0/chrom/random/index")



del_ins_size = 10
q_len = 1000
mutation_amount = 50


reference_pledge.set(ref_seq)
fm_index_pledge.set(fm_index)

for _ in range(100):
    for i in range(10):
        q = ""


        q_from = random.randint(0, ref_seq.unpacked_size_single_strand - q_len)
        q_to = q_from + q_len
        q = str(ref_seq.extract_from_to(q_from, q_to))
        for _ in range(mutation_amount):
            pos = random.randint(1,len(q)-1)
            q = q[:pos-1] + mutate(q[pos]) + q[pos:]

        for _ in range(mutation_amount):
            pos = random.randint(1,len(q)-11)
            l = random.randint(1,del_ins_size)
            q = q[:pos-1] + q[pos + l:]

        for _ in range(mutation_amount):
            pos = random.randint(1,len(q)-1)
            l = random.randint(1,del_ins_size)
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

        query_pledge[i].set(query)

    Pledge.simultaneous_get(result_pledges, 2)
print("done")
gc.collect()
print(gc.garbage)
del(result_pledges)
print(gc.garbage)
gc.collect()
print(gc.garbage)
print("gc.collect - done")