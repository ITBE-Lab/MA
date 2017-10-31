from LAuS import *
import random
import gc
import os

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]

def get_memory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since

def near(index, index_2):
    return index + 100 > index_2 and index - 100 < index_2


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

num_test = 100
query_pledge = []
for _ in range(num_test):
    query_pledge.append(Pledge(ContainerType.nucSeq))

reference_pledge = Pledge(ContainerType.packedNucSeq)
fm_index_pledge = Pledge(ContainerType.fM_index)

result_pledges = set_up_aligner_lr_segmentation(
    query_pledge,
    reference_pledge,
    fm_index_pledge
)



ref_seq = Pack()
#ref_seq.append("no name", "no desc", NucSeq("AACG"))
ref_seq.load("/mnt/ssd0/chrom/human/pack")
#fm_index = FMIndex(ref_seq)
fm_index = FMIndex()
fm_index.load("/mnt/ssd0/chrom/human/index")


memory = open("memory_usage.txt", "w")
memory.close()

del_ins_size = 10
q_len = 100


reference_pledge.set(ref_seq)
fm_index_pledge.set(fm_index)

while True:
    mutation_amount = 0#random.randint(0, 30)
    starts = []
    ends = []
    hits = 0
    memory = open("memory_usage.txt", "a")
    gc.collect()
    memory.write(str(get_memory()) + "\n")
    memory.close()
    for i in range(num_test):
        q = ""

        q_from = random.randint(0, ref_seq.unpacked_size_single_strand - q_len)
        starts.append(q_from)
        q_to = q_from + q_len
        ends.append(q_to)
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
        print(q)

        query = NucSeq(q)

        query_pledge[i].set(query)

    results = Pledge.simultaneous_get(result_pledges, 48)

    for i, alignment in enumerate(results):
        if alignment is None:
            continue
        if near(alignment.begin_on_ref(), starts[i]) and near(alignment.end_on_ref(), ends[i]):
            hits += 1
            print("hit")

    print(hits/num_test)
    break



print("done")