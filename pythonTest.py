from LAuS import *
import random
import gc
import os
import math
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.palettes import d3


def setup_multiple(module, in_list, perm):
    ret = []
    for in_ele in in_list:
        ret.append(module.promise_me((in_ele, perm)))
    return [ret]


human_genome = "/mnt/ssd0/chrom/human/all"



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




def test_chaining(seeds, strip, results):
    output_file("chaining_comp.html")

    p1 = figure(
            title="chaining comp", 
            x_axis_label='reference', 
            y_axis_label='query',
            plot_width=800,
            plot_height=800)

    listx = []
    listy = []
    stripx = []
    stripy = []
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
    for seed in strip:
        stripx.append( seed.end_ref() )
        stripy.append( seed.end() )

    p1.circle(listx, listy, size=7, color="black")
    p1.circle(stripx, stripy, size=10, legend="in strip of cons.", color="red")
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
        listx, listy = make_list_pretty(listx, listy)
        dashes = []
        for i in range(len(results)):
            if i == index:
                dashes.append(10)
                dashes.append(0)
            else:
                dashes.append(0)
                dashes.append(10)
        p1.line(listx, listy, line_width=3, line_dash=dashes, color=color, legend=name)

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

def make_up_seeds_simple():
    seeds = Seeds()
    num = 15
    for _ in range(num):
        q_pos = random.randint(0,200)
        ref_pos = random.randint(0,200)
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
    query = get_query(ref_seq, 1000, 10, 10)

    segments = LongestNonEnclosedSegments().execute((
            fm_index,
            query
        ))

    anchors = NlongestIntervalsAsAnchors(1).execute((segments,))

    bucketing = Bucketing()
    bucketing.max_hits = 5
    tail = Tail(Seeds())
    extractAll = ExtractAllSeeds()
    extractAll.max_hits = bucketing.max_hits

    strips = bucketing.execute((segments, anchors, query, ref_seq, fm_index))
    strip = tail.execute((strips,))

    #only get the best result for linesweep & bucketing
    ls_res = LineSweep().execute((strip,))

    #strip = extractAll.execute((segments, fm_index))

    ls2_res = LineSweep2().execute((strip,))

    seeds = segments.get_seeds(fm_index, 5)
    
    print(len(ls_res))
    print(len(ls2_res))


    test_chaining(seeds, strip, [(ls2_res, "green", "linesweep 2"), (ls_res, "blue", "linesweep")])

def print_runtime_breakdown(distances, runtimes_list, names, module_names):
    output_file("runtimes_breakdown.html")
    plots = []
    for index, runtimes in enumerate(runtimes_list):
        colors = d3["Category10"][10]
        p1 = figure(
                title=names[index],
                x_axis_label='distance',
                y_axis_label='seconds',
                plot_width=800,
                plot_height=800)

        for index2, runtime in enumerate(runtimes):
            p1.line(distances, runtime, legend=module_names[index][index2], line_width=5, color=colors[index2])

        plots.append(p1)

    show(row(plots))

def runtime_breakdown_helper(result_pledges):
    ret = []
    for index, step in enumerate(result_pledges):
        avg = 0
        for pledge in step:
            avg += pledge.exec_time
        avg /= len(step)
        ret.append(avg)
    return ret


def runtime_breakdown(num_test, query_pledge, result_pledges_list, names):
    distances = []

    #setup a new list
    module_names = []
    runtimes_list = []
    for result_pledges in result_pledges_list:
        runtimes = []
        names = []
        for pledge in result_pledges:
            runtimes.append([])
            names.append(pledge[0].get_pledger().get_name())
        names.append("total")
        runtimes.append([])
        runtimes_list.append(runtimes)
        module_names.append(names)

    max_dist = 15
    for distance in range(0,max_dist):
        print(str(distance) + "/" + str(max_dist))
        distances.append(distance)
        for i in range(num_test):
            query_pledge[i].set(get_query(ref_seq, 1000, distance, 10))

        for result_pledges in result_pledges_list:
            Pledge.simultaneous_get(result_pledges[-1], 10)

        for index1, result_pledges in enumerate(result_pledges_list):
            total = 0
            for index, runtime in enumerate(runtime_breakdown_helper(result_pledges)):
                runtimes_list[index1][index].append(runtime)
                total += runtime
            runtimes_list[index1][-1].append(total)

    print_runtime_breakdown(distances, runtimes_list, names, module_names)

##
# run smith waterman for all new samples within the database and save the results
#
def run_smw_for_all_samples(reference_name):
    query_pledges = []
    sample_ids = []

    #
    # collect samples and respective id's
    #
    print("collecting new queries for smw...")
    for sequence, sample_id in getNewQueries(db_name, "smw", human_genome, 100, 10):
        sample_ids.append(sample_id)
        query_pledges.append(Pledge(NucSeq()))
        query_pledges[-1].set(NucSeq(sequence))
    print("collected " + str(len(sample_ids)) + " new queries for smw")

    #
    # create pledges
    #
    reference_pledge = Pledge(Pack())
    ref_seq = Pack()
    ref_seq.load(reference_name)
    reference_pledge.set(ref_seq)
    result_pledge = setup_multiple(SMW(), query_pledges, reference_pledge)

    #
    # run smw for all new samples
    #
    print("running smw...")
    Pledge.simultaneous_get(result_pledge[-1], 16)
    print("done")

    #
    # collect the results
    #
    print("collecting and saving results...")
    results = []
    for alignment_vec, sample_id in zip(result_pledge[-1], sample_ids):
        for alignment in alignment_vec:
            result = (
                        sample_id,
                        alignment.get_score(),
                        alignment.begin_on_ref(),
                        alignment.length(),
                        "smw"
                    )
            results.append(result)

    #
    # save the results
    #
    submitResults(db_name, results)
    print("done")

#function

#exit()

def make_runtime_breakdown():
    print("setup...")
    num_test = 500
    query_pledge = []
    for _ in range(num_test):
        query_pledge.append(Pledge(NucSeq()))

    reference_pledge = Pledge(Pack())
    fm_index_pledge = Pledge(FMIndex())


    ref_seq = Pack()
    ref_seq.load("/mnt/ssd0/chrom/human/pack")
    fm_index = FMIndex()
    fm_index.load("/mnt/ssd0/chrom/human/index")

    reference_pledge.set(ref_seq)
    fm_index_pledge.set(fm_index)


    result_pledges = set_up_aligner(
        query_pledge,
        reference_pledge,
        fm_index_pledge,
        chain=LineSweep2(),
        seg=LongestLRSegments()
    )

    result_pledges_2 = set_up_aligner(
        query_pledge,
        reference_pledge,
        fm_index_pledge,
        LongestLRSegments()
    )

    print("done")
    """
    while True:
        for i in range(num_test):
            query_pledge[i].set(get_query(ref_seq, 100, random.randint(0,5), 10))

        Pledge.simultaneous_get(result_pledges[-1], 32)
        print("iteration!")
    """

    runtime_breakdown(num_test, query_pledge, 
        [result_pledges, result_pledges_2], ["Bs SoC sLs", "Bs SoC Ls"])
#function

#run_smw_for_all_samples(human_genome)

#creating samples int the database
createSampleQueries(human_genome, db_name, 1000, 100, 256)


#exit()

#compare_chaining_linesweep_visual()

