from LAuS import *
import random
import gc
import os
import math
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, layout
from bokeh.palettes import d3
from bokeh.models import LinearAxis, Range1d, LogColorMapper
from bokeh.models.formatters import FuncTickFormatter
from bokeh.io import save


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
    ref_seq.load(human_genome)
    fm_index = FMIndex()
    fm_index.load(human_genome)
    query_seq, query_id = getQueriesFor(db_name, human_genome, 50, 0, 100)[0]
    query = NucSeq(query_seq)

    query_orig = getOriginOf(db_name, query_id)

    print(query_seq)
    print(str(ref_seq.extract(query_orig, query_orig+100)))

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
#createSampleQueries(human_genome, db_name, 1000, 100, 256)

def test_my_approach(
            db_name,
            reference,
            name,
            seg=LongestNonEnclosedSegments(),
            chain=LineSweep(), 
            max_hits=5,
            num_anchors=10, 
            strips_of_consideration=True,
            low_res = False,
            re_seed = None,
            max_sweep = None,
            strip_size = 500
        ):
    print("collecting samples (" + name + ") ...")
    reference_pledge = Pledge(Pack())
    fm_index_pledge = Pledge(FMIndex())

    res1 = 1
    if low_res:
        res1 = 100
    res2 = 1
    if low_res:
        res2 = 5

    queries = getNewQueries(db_name, name, reference, res1, res2)
    #queries = getQueriesFor(db_name, reference, 40, 0, size)

    print("extracting " + str(len(queries)) + " samples (" + name + ") ...")
    #setup the query pledges
    query_pledge = []
    ids = []
    for sequence, sample_id in queries:
        query_pledge.append(Pledge(NucSeq()))
        query_pledge[-1].set(NucSeq(sequence))
        ids.append(sample_id)


    print("setting up (" + name + ") ...")
    ref_pack = Pack()
    ref_pack.load(reference)
    reference_pledge.set(ref_pack)

    fm_index = FMIndex()
    fm_index.load(reference)
    fm_index_pledge.set(fm_index)

    result_pledge = set_up_aligner(
        query_pledge,
        reference_pledge,
        fm_index_pledge,
        seg=seg,
        chain=LineSweep2(),
        max_hits=max_hits,
        num_anchors=num_anchors,
        strips_of_consideration=strips_of_consideration,
        re_seed=re_seed,
        max_sweep=max_sweep,
        strip_size=strip_size
    )

    #temp code
    if False:
        num = 0
        #result_pledge[-1][num].get()

        printer = AlignmentPrinter()

        print("origin: " + str(ids[num]))
        print("num results: " + str(len(result_pledge[-2][num].get().get())))

        for alignment in result_pledge[-2][num].get().get():
            printer.execute((alignment, query_pledge[num].get(), ref_pack))

        exit()

    print("computing (" + name + ") ...")
    Pledge.simultaneous_get(result_pledge[0], 32)
    Pledge.simultaneous_get(result_pledge[-1], 32)

    print("extracting results (" + name + ") ...")
    result = []
    for index in range(len(queries)):
        alignment = result_pledge[-1][index].get()
        segment_list = result_pledge[0][index].get()
        total_time = 0
        for step_index in range(len(result_pledge)):
            total_time += result_pledge[step_index][index].exec_time
        sample_id = queries[index][1]
        result.append(
                (
                    sample_id,
                    alignment.get_score(),
                    alignment.begin_on_ref(),
                    segment_list.num_seeds(fm_index, max_hits),
                    total_time,
                    name
                )
            )
    print("submitting results (" + name + ") ...")
    submitResults(db_name, result)
    print("done")


def test_my_approaches():
    #clearResults(db_name, human_genome, "pBs:50000 SoC:5 sLs")
    #test_my_approach(db_name, human_genome, "pBs:50000 SoC:5 sLs", max_hits=50000)

    #clearResults(db_name, human_genome, "pBs:5 SoC:10 sLs")
    #test_my_approach(db_name, human_genome, "pBs:5 SoC:10 sLs", num_anchors=10)

    #clearResults(db_name, human_genome, "Bs:50000 SoC:5 sLs")
    #test_my_approach(db_name, human_genome, "Bs:50000 SoC:5 sLs", max_hits=50000, seg=LongestLRSegments())

    clearResults(db_name, human_genome, "pBs:5 SoC:100,500nt sLs:5")
    test_my_approach(db_name, human_genome, "pBs:5 SoC:100,500nt sLs:5", num_anchors=100, max_sweep=5)

    clearResults(db_name, human_genome, "pBs:5 SoC:100,500nt sLs:1")
    test_my_approach(db_name, human_genome, "pBs:5 SoC:100,500nt sLs:1", num_anchors=100, max_sweep=1)

    clearResults(db_name, human_genome, "pBs:5 SoC:1,500nt sLs:inf")
    test_my_approach(db_name, human_genome, "pBs:5 SoC:1,500nt sLs:inf", num_anchors=1)

    clearResults(db_name, human_genome, "pBs:5 SoC:1,100nt sLs:inf")
    test_my_approach(db_name, human_genome, "pBs:5 SoC:1,100nt sLs:inf", num_anchors=1, strip_size=100)

    clearResults(db_name, human_genome, "pBs:5 SoC:1,10000nt sLs:inf")
    test_my_approach(db_name, human_genome, "pBs:5 SoC:1,10000nt sLs:inf", num_anchors=1, strip_size=10000)

    clearResults(db_name, human_genome, "pBs:5 SoC:1,1000000nt sLs:inf")
    test_my_approach(db_name, human_genome, "pBs:5 SoC:1,1000000nt sLs:inf", num_anchors=1, strip_size=1000000)

    #clearResults(db_name, human_genome, "pBsr:1000 SoC:10 sLs")
    #test_my_approach(db_name, human_genome, "pBsr:1000 SoC:10 sLs", max_hits=1000, num_anchors=10, re_seed=0)

    #clearResults(db_name, human_genome, "Bs:5 SoC:5 sLs:inf")
    #test_my_approach(db_name, human_genome, "Bs:5 SoC:5 sLs:inf", seg=LongestLRSegments())

    #clearResults(db_name, human_genome, "pBs:5 SoC:1000 sLs:5")
    #test_my_approach(db_name, human_genome, "pBs:5 SoC:1000 sLs:5", num_anchors=1000, max_sweep=5)

    #clearResults(db_name, human_genome, "pBs:5 Ch")
    #test_my_approach(db_name, human_genome, "pBs:5 Ch", chain=Chaining(), strips_of_consideration=False)

    #clearResults(db_name, human_genome, "Bs:10000 SoC:10000 sLs")
    #test_my_approach(db_name, human_genome, "Bs:10000 SoC:10000 sLs", seg=LongestLRSegments(),
    #                num_anchors=10000, max_hits=10000)



def analyse_all_approaches(query_size = 100, indel_size = 10):
    output_file("results_log.html")
    approaches = getApproachesWithData(db_name)
    plots = [ [], [], [], [] ]

    for approach_ in approaches:
        approach = approach_[0]
        results = getResults(db_name, approach, query_size, indel_size)
        hits = {}
        tries = {}
        run_times = {}
        scores = {}
        nums_seeds = {}

        max_indel = 2*query_size/indel_size
        max_mut = query_size

        def init(d, x, y):
            if x not in d:
                d[x] = {}
            if y not in d[x]:
                d[x][y] = 0
            return d

        for score, result_start, original_start, num_mutation, num_indels, num_seeds, run_time in results:
            hits = init(hits, num_mutation, num_indels)
            tries = init(tries, num_mutation, num_indels)
            run_times = init(run_times, num_mutation, num_indels)
            if not score is None:
                scores = init(scores, num_mutation, num_indels)
            nums_seeds = init(nums_seeds, num_mutation, num_indels)

            if near(result_start, original_start):
                hits[num_mutation][num_indels] += 1
            tries[num_mutation][num_indels] += 1
            run_times[num_mutation][num_indels] += run_time
            nums_seeds[num_mutation][num_indels] += num_seeds
            if not score is None:
                scores[num_mutation][num_indels] += score


        def makePicFromDict(d, w, h, divideBy, title, ignore_max_n = 0, log = False, set_max = None, set_min=None):
            pic = []
            min_ = 10000.0
            max_ = 0.0
            if len(d.keys()) == 0:
                return None
            else:
                w_keys = sorted(d.keys())
                h_keys = sorted(d[0].keys())

                for x in w_keys:
                    pic.append( [] )
                    for y in h_keys:
                        if x not in d or y not in d[x]:
                            pic[-1].append( float("nan") )
                        else:
                            pic[-1].append( d[x][y] / divideBy[x][y] )
                            if pic[-1][-1] < min_:
                                min_ = pic[-1][-1]
                            if pic[-1][-1] > max_:
                                max_ = pic[-1][-1]

                for _ in range(ignore_max_n):
                    max_x = 0
                    max_y = 0
                    for x, row in enumerate(pic):
                        for y, p in enumerate(row):
                            if p > pic[max_x][max_y]:
                                max_x = x
                                max_y = y
                    pic[max_x][max_y] = float('nan')
                if ignore_max_n > 0:
                    max_ = 0
                    for row in pic:
                        for p in row:
                            if p > max_:
                                max_ = p

            if set_max is not None:
                max_ = set_max
            if set_min is not None:
                min_ = set_min

            color_mapper = LinearColorMapper(palette="Viridis256", low=min_, high=max_)
            if log:
                color_mapper = LogColorMapper(palette="Viridis256", low=min_, high=max_)

            tick_formater = FuncTickFormatter(code="""
                return Math.max(Math.floor( (tick+1)/2),0) + '; ' +
                        Math.max(Math.floor( (tick)/2),0)"""
                )
            #tick_formater = FuncTickFormatter(code="return 'a')

            plot = figure(title=title,
                    x_range=(0,h), y_range=(0,w),
                    x_axis_label='num ' + str(indel_size) + ' nt insertions; num ' + str(indel_size) + ' nt deletions', y_axis_label='num mutations'
                )
            plot.xaxis.formatter = tick_formater
            plot.image(image=[pic], color_mapper=color_mapper,
                    dh=[w], dw=[h], x=[0], y=[0])

            color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))

            plot.add_layout(color_bar, 'left')

            return plot

        avg_hits = makePicFromDict(hits, max_mut, max_indel, tries, "accuracy " + approach, set_max=1, set_min=0)
        avg_runtime = makePicFromDict(run_times, max_mut, max_indel, tries, "runtime " + approach, 0, True)
        avg_score = makePicFromDict(scores, max_mut, max_indel, tries, "score " + approach)
        avg_seeds = makePicFromDict(nums_seeds, max_mut, max_indel, tries, "num seeds " + approach)

        if not avg_hits is None:
            plots[0].append(avg_hits)
        if not avg_runtime is None:
            plots[1].append(avg_runtime)
        if not avg_score is None:
            plots[2].append(avg_score)
        if not avg_seeds is None:
            plots[3].append(avg_seeds)

    save(layout(plots))

#exit()
def manualCheckSequences():
    ref_seq = Pack()
    ref_seq.load(human_genome)
    for query_seq, query_id in getQueriesFor(db_name, human_genome, 40, 0, 100):
        query = NucSeq(query_seq)

        query_orig = getOriginOf(db_name, query_id)

        print(query_seq)
        print(str(ref_seq.extract_from_to(query_orig, query_orig+100)))
        print("")

#manualCheckSequences()
#exit()

#exit()
#createSampleQueries(human_genome, db_name, 1000, 100, 64, True)
#test_my_approaches()
analyse_all_approaches(1000, 100)
#compare_chaining_linesweep_visual()
