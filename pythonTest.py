from LAuS import *
import random
import gc
import os
import math
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column


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


def make_list_pretty(lx, ly):
    lastx = 0
    lasty = 0
    for i in range(len(lx)):
        if lx[i] is float('nan'):
            lastx = 0
            lasty = 0
            continue
        if lx[i] < lastx or ly[i] < lasty:
            add = max(lastx - lx[i], lasty - ly[i])
            lx[i] += add
            ly[i] += add
        else:
            lastx = lx[i]
            lasty = ly[i]
    return lx, ly


def get_query(ref_seq, q_len, mutation_amount, del_ins_size):
    q = ""

    q_from = random.randint(0, ref_seq.unpacked_size_single_strand - q_len)
    q_to = q_from + q_len
    q = str(ref_seq.extract_from_to(q_from, q_to))
    for _ in range(mutation_amount * 10):
        pos = random.randint(1,len(q)-1)
        q = q[:pos-1] + mutate(q[pos]) + q[pos:]

    for _ in range(mutation_amount):
        if len(q) <= del_ins_size:
            break
        pos = random.randint(1,len(q)-del_ins_size - 1)
        l = del_ins_size
        q = q[:pos-1] + q[pos + l:]

    for _ in range(mutation_amount):
        pos = random.randint(1,len(q)-1)
        l = del_ins_size
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

    return NucSeq(q)

def test_chaining(seeds, results):
    output_file("chaining_comp.html")

    p1 = figure(
            title="chaining comp", 
            x_axis_label='reference', 
            y_axis_label='query',
            plot_width=800,
            plot_height=800)

    listx = []
    listy = []
    linex = []
    liney = []
    for seed in seeds:
        listx.append( seed.end_ref() )
        listy.append( seed.end() )

        linex.append( seed.start_ref() )
        liney.append( seed.start() )
        linex.append( seed.end_ref() )
        liney.append( seed.end() )
        linex.append( float('nan') )
        liney.append( float('nan') )

    p1.circle(listx, listy, size=10, color="black")
    p1.line(linex, liney, line_width=5 ,color="black")

    for index, result_tuple in enumerate(results):
        result, color, name = result_tuple
        listx = []
        listy = []
        listx_ = []
        listy_ = []
        for seed in result:
            listx.append(seed.start_ref())
            listy.append(seed.start())
            listx.append(seed.end_ref())
            listy.append(seed.end())
            listx_.append(seed.end_ref())
            listy_.append(seed.end())
            listx.append( float('nan') )
            listy.append( float('nan') )
        listx, listy = make_list_pretty(listx, listy)
        dashes = []
        for i in range(len(results)):
            if i == index:
                dashes.append(10)
                dashes.append(0)
            else:
                dashes.append(0)
                dashes.append(10)
        p1.line(listx, listy, line_width=5, line_dash=dashes, color=color, legend=name)
        p1.circle(listx_, listy_, size=10, color=color)

    show(p1)

def make_up_seeds():
    seeds = Seeds()
    num = 15
    for _ in range(num):
        q_pos = random.randint(0,200)
        ref_pos = 300000000 * random.randint(0,num/5) + random.randint(-10,10) + q_pos
        seeds.append(Seed(
            q_pos,#query
            random.randint(30,30),#length
            ref_pos))#ref
    return seeds

def compare_chaining_linesweep_visual():
    ref_seq = Pack()
    ref_seq.load("/mnt/ssd0/chrom/human/pack")
    fm_index = FMIndex()
    fm_index.load("/mnt/ssd0/chrom/human/index")
    query = get_query(ref_seq, 1000, 5, 10)

    segments = LongestLRSegments().execute((
            fm_index,
            query
        ))

    anchors = NlongestIntervalsAsAnchors(10).execute((segments,))

    bucketing = Bucketing()
    bucketing.max_hits = 5
    strips = bucketing.execute((segments, anchors, query, ref_seq, fm_index))

    ls_res = SweepAllReturnBest(LineSweep()).execute((strips,))
    ch_res = SweepAllReturnBest(Chaining()).execute((strips,))

    seeds = segments.get_seeds(fm_index, 5)
    """

    #seeds = make_up_seeds()

    ls_res = LineSweep().execute((
            seeds,
        ))

    ch_res = Chaining().execute((
            seeds,
        ))"""

    test_chaining(seeds, [(ch_res, "green", "chaining"), (ls_res, "blue", "linesweep")])


compare_chaining_linesweep_visual()
exit()

q = ""

num_test = 1000
query_pledge = []
for _ in range(num_test):
    query_pledge.append(Pledge(ContainerType.nucSeq))

reference_pledge = Pledge(ContainerType.packedNucSeq)
fm_index_pledge = Pledge(ContainerType.fM_index)




ref_seq = Pack()
#ref_seq.append("no name", "no desc", NucSeq("AACG"))
ref_seq.load("/mnt/ssd0/chrom/human/pack")
#fm_index = FMIndex(ref_seq)
fm_index = FMIndex()
fm_index.load("/mnt/ssd0/chrom/human/index")


memory = open("memory_usage.txt", "w")
memory.close()

del_ins_size = 10
q_len = 1000

# output to static HTML file
output_file("result.html")

# create a new plot with a title and axis labels
p = figure(title="quality comparison BWA", x_axis_label='genetic distance', y_axis_label='hit rate')

reference_pledge.set(ref_seq)
fm_index_pledge.set(fm_index)

for max_h in range(100,1000, 100):
    print("max_hits = " + str(max_h))
    result_pledges = set_up_aligner(
        query_pledge,
        reference_pledge,
        fm_index_pledge,
        max_hits=max_h,
        chain=Chaining()
    )

    x = []
    y = []

    #while True:
    for mutation_amount in range(0,30):
        #mutation_amount = 0#random.randint(0, 30)
        x.append(mutation_amount)
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
            for _ in range(mutation_amount * 10):
                pos = random.randint(1,len(q)-1)
                q = q[:pos-1] + mutate(q[pos]) + q[pos:]

            for _ in range(mutation_amount):
                if len(q) <= del_ins_size:
                    break
                pos = random.randint(1,len(q)-del_ins_size - 1)
                l = del_ins_size
                q = q[:pos-1] + q[pos + l:]

            for _ in range(mutation_amount):
                pos = random.randint(1,len(q)-1)
                l = del_ins_size
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

        Pledge.simultaneous_get(result_pledges[-1], 48)

        for i, pledge in enumerate(result_pledges[-1]):
            alignment = pledge.get()
            if alignment is None:
                continue
            if near(alignment.begin_on_ref(), starts[i]) and near(alignment.end_on_ref(), ends[i]):
                hits += 1
                #print("hit")
        y.append(hits/num_test)
        print(hits/num_test)

    cdgt1 = math.floor(10*max_h/1000)
    cdgt2 = math.floor(100*max_h/1000)%10
        
    # add a line renderer with legend and line thickness
    p.line(x, y, legend=str(max_h), color="#A0A0" + str(cdgt1) + str(cdgt2), line_width=3)
    #break




# show the results
show(p)

print("done")