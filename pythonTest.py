from MA import *
import random
import gc
import os
import math
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column, layout
from bokeh.palettes import d3
from bokeh.models import LinearAxis, Range1d, LogColorMapper, FixedTicker, BasicTicker, Grid
from bokeh.models.formatters import FuncTickFormatter
from bokeh.layouts import gridplot
from bokeh.io import save
import numpy as np
import colorsys
from statistics import mean, median, stdev
from sys import stderr

# Computes a single rainbow-color. (Like in C++ and Java-code)
def rainbow(value):
    minimum = 0
    maximum = 1
    """
    Local helper function for rgb-value computation.
    """
    def clamp(x):
        """
        Local helper function for min-max encapsulation.
        """
        return max(0, min(x, 255))

    value = min(max(minimum, value), maximum)
    difference = maximum - minimum
    red = green = blue = 1.0
    if value < minimum + 0.25 * difference:
        red = 0.0
        green = 4.0 * (value - minimum) / difference
    else:
        if value < minimum + 0.5 * difference:
            red = 0.0
            blue = 1.0 + 4.0 * (minimum + 0.25 * difference - value) / difference
        else:
            if value < minimum + 0.75 * difference:
                red = 4.0 * (value - minimum - 0.5 * difference) / difference
                blue = 0.0
            else:
                green = 1.0 + 4.0 * (minimum + 0.75 * difference - value) / difference
                blue = 0.0
    return (red, green, blue)


def light_spec_approximation(x):
    #map input [0, 1] to wavelength [350, 645]
    w = 370 + x * (645-370)
    r = 0.0
    g = 0.0
    b = 0.0
    if w < 440:
        r = -(w - 440.) / (440. - 380.)
        b = 1.0
    elif w >= 440 and w < 490:
        g = (w - 440.) / (490. - 440.)
        b = 1.0
    elif w >= 490 and w < 510:
        g = 1.0
        b = -(w - 510.) / (510. - 490.)
    elif w >= 510 and w < 580:
        r = (w - 510.) / (580. - 510.)
        g = 1.0
    elif w >= 580 and w < 645:
        r = 1.0
        g = -(w - 645.) / (645. - 580.)
    elif w >= 645:
        r = 1.0

    #intensity
    i = 1.0
    if w > 650:
        i = .3 + .7*(780-w)/(780-650)
    elif w < 420:
        i = .3 + .7*(w-380)/(420-380)

    #gamma
    m = .8

    return (i*r**m,i*g**m,i*b**m)

def hsv_walk(x):
    return colorsys.hsv_to_rgb(-x*.8+.8, 1, x*.5+.5)

def light_spec_approximation_2(x):
    #map input [0, 1] to wavelength [350, 645]
    w = 370 + x * (645-370)
    r = 0.0
    g = 0.0
    b = 0.0
    if w < 440:
        r = -(w - 440.) / (440. - 380.)
        b = 1.0
    elif w >= 440 and w < 490:
        g = (w - 440.) / (490. - 440.)
        b = 1.0
    elif w >= 490 and w < 580:
        g = 1.0
        b = max(0.0, -(w - 510.) / (510. - 490.))
        r = max(0.0, (w - 510.) / (580. - 510.))
    elif w >= 580 and w < 645:
        r = 1.0
        g = -(w - 645.) / (645. - 580.)
    elif w >= 645:
        r = 1.0

    #intensity
    i = 1.0
    if w > 650:
        i = .3 + .7*(780-w)/(780-650)
    elif w < 420:
        i = .3 + .7*(w-380)/(420-380)

    #gamma
    m = .8

    return (i*r**m,i*g**m,i*b**m)

def heatmap_palette(scheme, num_colors):
    def format(rgb):
        def clamp(x):
            return max(0, min(x, 255))
        red, green, blue = rgb
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                               clamp(int(blue * 255)))
    return [format(scheme(x)) for x in np.linspace(0, 1, num_colors)]


def setup_multiple(module, in_list, perm):
    ret = []
    for in_ele in in_list:
        ret.append(module.promise_me((in_ele, perm)))
    return [ret]


human_genome = "/mnt/ssd0/genome/human"



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
    gc.collect()
    return _VmB('VmSize:') - since

def near(start, start_2, end, end_2):
    max_d = 0#10000
    return end_2 + max_d >= start and start_2 - max_d <= start


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

    segments = BinarySeeding().execute((
            fm_index,
            query
        ))

    anchors = GetAnchors(1).execute((segments,))

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


def memory_test(reference, test_index):
    reference_pledge = Pledge(Pack())
    fm_index_pledge = Pledge(FMIndex())
    query_pledge = Pledge(NucSeq())
    
    ref_pack = Pack()
    ref_pack.load(reference)
    reference_pledge.set(ref_pack)

    fm_index = FMIndex()
    fm_index.load(reference)
    fm_index_pledge.set(fm_index)

    result_pledge = set_up_aligner(
        query_pledge,
        reference_pledge,
        fm_index_pledge
    )
    module = BinarySeeding(True)
    
    mem = None

    num = 0
    print("running")
    while True:
        q_from, q, original_nuc_dist, modified_nuc_dist = get_query(ref_pack, 200, 0, 0, 1)
        query_pledge.set(NucSeq(q))

        #module.execute(fm_index, NucSeq(q))

        result_pledge[test_index].get()
        #Pledge.simultaneous_get( [result_pledge[test_index]] , 1)

        #query_pledge.get()

        if not mem is None and not get_memory(mem) == 0:
            print(get_memory(mem), "\t", get_memory())
            #print(q)
        mem = get_memory()

        num += 1
        if num >= 1000:
            print("marker")
            num = 0

## @brief Yield successive n-sized chunks from l.
def chunks(l, n):
    if n >= len(l):
        yield l
    else:
        for i in range(0, len(l), n):
            yield l[i:i + n]

#run_smw_for_all_samples(human_genome)

#creating samples int the database
#createSampleQueries(human_genome, db_name, 1000, 100, 256)

def test_my_approach(
            db_name,
            reference,
            name,
            max_hits=100,
            num_strips=10, 
            complete_seeds = False
        ):
    print("collecting samples (" + name + ") ...")

    all_queries = getNewQueries(db_name, name, reference)
    #queries = getQueriesFor(db_name, reference, 40, 0, size)

    print("having ", len(all_queries), " samples total (", name, ") ...")
    if len(all_queries) == 0:
        print("no new queries ... done")
        return

    aligner = Aligner(max_hits, num_strips, complete_seeds)

    runtimes = {}

    extract_size = len(all_queries) #2**15
    # break samples into chunks of 2^14
    for index, queries in enumerate(chunks(all_queries, extract_size)):
        print("extracting", len(queries), "samples", index, "/",
            len(all_queries)/extract_size, "(", name, ") ...")
        #setup the query pledges
        query_container = ContainerVector(NucSeq())
        del query_container[:]
        ids = []
        optimal_alignment_in = []
        for position, data in enumerate(queries):
            sequence, sample_id = data
            query_container.append(NucSeq(sequence))
            query_container[-1].name = str(position)
            ids.append(sample_id)
            optimal_alignment_in.append( (Pledge(NucSeq()),Pledge(NucSeq())) )

        print("setting up (", name, ") ...")

        ref_pack = Pack()
        ref_pack.load(reference)
        aligner.setRef(ref_pack)

        fm_index = FMIndex()
        fm_index.load(reference)
        aligner.setInd(fm_index)

        print("computing (", name, ") ...")
        results = aligner.align(query_container)

        do_optimal = False
        optimal_alignment = None
        if do_optimal:
            smw = SMW()
            optimal_alignment_out = []
            for a, b in optimal_alignment_in:
                optimal_alignment_out.append(smw.promise_me(a, b))

            print("computing optimal (", name, ") ...")
            for i in range(len(queries)):
                alignment = results[i]
                index = int(alignment.stats.name)
                query = queries[index][0]
                optimal_alignment_in[index][0].set(NucSeq(query[alignment.begin_on_query:alignment.end_on_query]))
                if alignment.begin_on_ref() != alignment.end_on_ref():
                    optimal_alignment_in[index][1].set(
                        ref_pack.extract_from_to(alignment.begin_on_ref(),alignment.end_on_ref()))
                else:
                    optimal_alignment_in[index][1].set(NucSeq(""))
            Pledge.simultaneous_get(optimal_alignment_out, 32)

        print("extracting results (", name, ") ...")
        result = []
        for i in range(len(queries)):
            # @todo change sorting order
            alignment = results[i]
            alignment2 = None
            #print(alignment.stats.name)
            if alignment.stats.name == "unknown":
                continue
            index = int(alignment.stats.name)

            total_time = 0
            for key, value in aligner.get_runtimes().items():
                if not key in runtimes:
                    runtimes[key] = 0.0
                runtimes[key] += value
                total_time += value

            sample_id = queries[index][1]

            max_nmw_area=0
            max_diag_deviation = 0.0
            curr_diag_deviation = 0
            curr_max_diag_deviation = 0
            gap_size = 0
            cur_nmw_h = 0
            cur_nmw_w = 0

            if optimal_alignment != None and alignment.get_score() > optimal_alignment.get_score():
                print("WARNING: alignment computed better than optimal score")
                query = queries[index][0]
                AlignmentPrinter().execute(
                        alignment, 
                        NucSeq(query[alignment.begin_on_query:alignment.end_on_query]),
                        ref_pack
                    )

            """
            for pos in range(len(alignment)):
                match_type = alignment[pos]
                if match_type == MatchType.seed:
                    # @todo really filter out all gaps smaller than 10?
                    # NO: we need to filter out all gaps so small that the pure missmatchscore 
                    #     must be smaller than the gap score 
                    if gap_size > 10 and float(abs(curr_max_diag_deviation)) / gap_size > max_diag_deviation:
                        max_diag_deviation = float(abs(curr_max_diag_deviation)) / gap_size
                    curr_diag_deviation = 0
                    curr_max_diag_deviation = 0
                    gap_size = 0
                    cur_nmw_h = 0
                    cur_nmw_w = 0
                else:
                    gap_size += 1
                    if match_type == MatchType.insertion:
                        curr_diag_deviation += 1
                        cur_nmw_h += 1
                    elif match_type == MatchType.deletion:
                        curr_diag_deviation -= 1
                        cur_nmw_w += 1
                    elif match_type == MatchType.missmatch:
                        cur_nmw_w += 1
                        cur_nmw_h += 1
                    if curr_diag_deviation > curr_max_diag_deviation:
                        curr_max_diag_deviation = curr_diag_deviation
                    if cur_nmw_h * cur_nmw_w > max_nmw_area:
                        max_nmw_area = cur_nmw_h * cur_nmw_w
            """

            sc2 = None
            if not optimal_alignment is None:
                sc2 = optimal_alignment.get_score()
            score2 = 0
            if not alignment2 is None:
                score2 = alignment2.get_score()

            result.append(
                    (
                        sample_id,
                        alignment.get_score(),
                        score2,
                        sc2,
                        alignment.begin_on_ref(),
                        alignment.end_on_ref(),
                        0,#segment_list.num_seeds(fm_index, max_hits),
                        alignment.stats.index_of_strip,
                        alignment.stats.seed_coverage,
                        alignment.stats.num_seeds_in_strip,
                        alignment.stats.anchor_size,
                        alignment.stats.anchor_ambiguity,
                        max_diag_deviation,
                        max_nmw_area,
                        alignment.mapping_quality,
                        total_time,
                        name
                    )
                )
        print("submitting results (", name, ") ...")
        submitResults(db_name, result)
    print("total runtimes:")
    for key, value in runtimes.items():
        print(value, "(", key, ")")
    print("done")


def test_my_approaches_rele(db_name):
    test_my_approach(db_name, human_genome, "non-enclosed pairs", seg=BinarySeeding(True), num_anchors=200, nmw_give_up=20000)

    test_my_approach(db_name, human_genome, "non-enclosed", seg=BinarySeeding(False), num_anchors=200, nmw_give_up=20000)

    #@todo BWA-MEM seeding technique

    seg = BinarySeeding(False)
    seg.do16ntevery10ntExtension = True
    test_my_approach(db_name, human_genome, "BOWTIE 2", seg=seg, num_anchors=200, nmw_give_up=20000)

    seg2 = BinarySeeding(False)
    seg2.blasrExtension = True
    test_my_approach(db_name, human_genome, "BLASR", seg=seg2, num_anchors=200, nmw_give_up=20000)

def test_my_approaches(db_name):
    # [p]Bs:<max_ambiguity>,*<max_seeds_fac> 
    # # SoC:<num_SoC>,<SoC_size>nt,\<<max_ambiguity>,\><min_seeds_in_strip>,*<min_coverage>
    # sLs:<num_nmw>

    #
    # this is the un optimized hammer method
    #

    clearResults(db_name, human_genome, "MA 1")
    clearResults(db_name, human_genome, "MA 2")
    #clearResults(db_name, human_genome, "MA 3")
    #clearResults(db_name, human_genome, "MA NMW-Band = 100")

    #test_my_approach(db_name, human_genome, "MA MAX QUALITY", num_anchors=10000, seg=BinarySeeding(False), max_sweep=500, min_seeds=2, min_seed_length=0.01, max_seeds=6.0, max_seeds_2=7.0, nmw_give_up=0)
    #
    # optimized in a way that speed is maximal without reducing accuracy by filters (hopefully)
    #
    #test_my_approach(db_name, human_genome, "MA 3", num_strips=1000, max_sweep=1000, nmw_give_up=0, max_hits=100)

    test_my_approach(db_name, human_genome, "MA 2", max_hits=100, num_strips=10, complete_seeds=True)

    #test_my_approach(db_name, human_genome, "Bs,SoC,sLs_quality&speed", num_anchors=200, max_sweep=0, seg=BinarySeeding(True), min_seeds=2, min_seed_length=0.02, max_seeds=0, max_seeds_2=0.17, nmw_give_up=7500)

    # pretty good mabs 1
    # min_seeds=2, min_seed_length=0.4, max_seeds=0, max_seeds_2=0.15,
    test_my_approach(db_name, human_genome, "MA 1", max_hits=100, num_strips=5, complete_seeds=False)
    #test_my_approach(db_name, human_genome, "MA 1", num_anchors=1000, seg=BinarySeeding(True))

    #clearResults(db_name, human_genome, "MA 2 radix")
    #test_my_approach(db_name, human_genome, "MA 2", num_anchors=200, max_sweep=0, seg=BinarySeeding(True), min_seeds=2, min_seed_length=0.02, max_seeds=0, max_seeds_2=0.15, nmw_give_up=1000)

    #test_my_approach(db_name, human_genome, "MA 2.1", num_anchors=1000, max_sweep=0, seg=BinarySeeding(True), min_seeds=2, min_seed_length=0.02, max_seeds=0, max_seeds_2=0.15, nmw_give_up=5000)

    #clearResults(db_name, human_genome, "Bs,SoC,sLs_quality")
    #test_my_approach(db_name, human_genome, "Bs,SoC,sLs_quality", num_anchors=10000, max_sweep=0, seg=BinarySeeding(True), min_seeds=0, min_seed_length=0.02, max_seeds=0, max_seeds_2=0, nmw_give_up=0)

def filter_suggestion(db_name, min_percentage_cov=0.1, max_num_seeds=5000, min_num_seed_in_strip=3):
    approaches = getApproachesWithData(db_name)
    plots = []

    for approach_ in approaches:
        true_pos = 0
        true_neg = 0
        false_pos = 0
        false_neg = 0

        approach = approach_[0]
        for result in getResults(db_name, approach):
            score, result_start, original_start, num_mutation, num_indels, num_seeds, run_time, index_of_chosen_strip, seed_coverage_chosen_strip, num_seeds_chosen_strip, anchor_size, anchor_ambiguity, max_diag_deviation, max_nmw_area, original_size = result

            is_filtered = True

            if seed_coverage_chosen_strip / original_size >= min_percentage_cov:
                is_filtered = False
            if num_seeds <= max_num_seeds:
                is_filtered = False
            if num_seeds_chosen_strip >= min_num_seed_in_strip:
                is_filtered = False

            is_kept = not is_filtered

            if near(result_start, original_start) and is_kept:
                true_pos += 1
            elif near(result_start, original_start) and not is_kept:
                false_neg += 1
            elif not near(result_start, original_start) and is_kept:
                false_pos += 1
            elif not near(result_start, original_start) and not is_kept:
                true_neg += 1

        print("======",approach,"======")
        print("true  positives:\t", true_pos, "\t*very good*")
        print("false negatives:\t", false_neg, "\t*very bad*")
        print("")
        print("false positives:\t", false_pos, "\t*bad*")
        print("true  negatives:\t", true_neg, "\t*good*")


def analyse_detailed(out_prefix, db_name):
    approaches = getApproachesWithData(db_name)
    plots = []

    for approach_ in approaches:
        strip_index = ([], [])
        seed_coverage = ([], [])
        num_seeds_tot = ([], [])
        num_seeds_strip = ([], [])
        anc_size = ([], [])
        anc_ambiguity = ([], [])
        max_diag_deviations = ([], [])
        max_nmw_areas = ([], [])
        mapping_qual = ([], [])
        score_dif = ([], [])

        approach = approach_[0]
        for result in getResults(db_name, approach):
            score, score2, optimal_score_this_region, result_start, result_end, original_start, num_mutation, num_indels, num_seeds, mapping_quality, run_time, index_of_chosen_strip, seed_coverage_chosen_strip, num_seeds_chosen_strip, anchor_size, anchor_ambiguity, max_diag_deviation, max_nmw_area, original_size = result

            index = 1
            if near(result_start, original_start, result_end, original_start+original_size):
                index = 0

            mapping_qual[index].append(mapping_quality)
            if not score is None and not score2 is None:
                score_dif[index].append(abs(score - score2))
            else:
                score_dif[index].append(score)
            strip_index[index].append(index_of_chosen_strip)
            if not seed_coverage_chosen_strip is None:
                seed_coverage[index].append(seed_coverage_chosen_strip / original_size)
            num_seeds_tot[index].append(num_seeds)
            num_seeds_strip[index].append(num_seeds_chosen_strip)
            anc_size[index].append(anchor_size)
            anc_ambiguity[index].append(anchor_ambiguity)
            max_diag_deviations[index].append(max_diag_deviation)
            max_nmw_areas[index].append(max_nmw_area)

        #print("required_nmw_area", max_diag_deviations)

        output_file(out_prefix + approach + ".html")

        def bar_plot(data, name, num_buckets=20, print_data=False):
            plot = figure(title=name,
                    x_axis_label='values',
                    y_axis_label='relative amount'
                )

            total_amount_1 = len(data[0])
            total_amount_2 = len(data[1])
            min_val = 10000
            max_val = 0
            for x in data[0]:
                if x is None:
                    return plot
                if x < min_val:
                    min_val = x
                if x > max_val:
                    max_val = x
            for x in data[1]:
                if x is None:
                    return plot
                if x < min_val:
                    min_val = x
                if x > max_val:
                    max_val = x
            bucket_width = float(max_val - min_val) / num_buckets
            #print(max_val, min_val, bucket_width)
            buckets_1 = []
            buckets_2 = []
            buckets_x = []
            for index in range(num_buckets):
                buckets_1.append(0.0)
                buckets_2.append(0.0)
                buckets_x.append(bucket_width * index + min_val + bucket_width/2)
            for x in data[0]:
                buckets_1[int( float(num_buckets*(x - min_val)) / ( (max_val+0.01)-min_val))] += 1.0 / total_amount_1
            for x in data[1]:
                buckets_2[int( float(num_buckets*(x - min_val)) / ((max_val+0.01)-min_val))] += 1.0 / total_amount_2
            if print_data:
                print(buckets_1)
                print(buckets_2)
            plot.vbar(x=buckets_x, width=bucket_width, bottom=0, top=buckets_1, color="blue", legend="accurate", fill_alpha=0.5)
            plot.vbar(x=buckets_x, width=bucket_width, bottom=0, top=buckets_2, color="red", legend="inaccurate", fill_alpha=0.5)
            return plot

        def true_pos_true_neg(data, name, amount=100):
            plot = figure(title=name,
                    y_axis_label='true positive rate',
                    x_axis_label='false positive rate'
                )

            total_amount_1 = len(data[0])
            total_amount_2 = len(data[1])
            min_val = 1
            max_val = 0
            for x in data[0]:
                if x is None:
                    return plot
                if x < min_val:
                    min_val = x
                if x > max_val:
                    max_val = x
            for x in data[1]:
                if x is None:
                    return plot
                if x < min_val:
                    min_val = x
                if x > max_val:
                    max_val = x
            adjust = max_val-min_val
            max_val += adjust*0.1
            min_val -= adjust*0.1

            values = []
            for i in range(amount):
                values.append(min_val + i*adjust/float(amount))
            line_x = []
            line_y = []
            for threshold in values:
                true_pos = 0
                false_pos = 0
                for ele in data[0]:
                    if ele > threshold:
                        true_pos += 1
                for ele in data[1]:
                    if ele > threshold:
                        false_pos += 1
                line_x.append(false_pos / total_amount_2)
                line_y.append(true_pos / total_amount_1)
            plot.line(line_x, line_y)
            #print(min_val, max_val, line_x, line_y)
            return plot

        save(column([
            true_pos_true_neg(mapping_qual, "mapping quality"),
            bar_plot(score_dif, "first and second score diff", 50),
            bar_plot(strip_index, "index of Strip of Consideration", 50),
            bar_plot(seed_coverage, "percentage of query covered by seeds in strip", 50),
            bar_plot(num_seeds_tot, "number of seeds total", 50),
            bar_plot(num_seeds_strip, "number of seeds in strip", 25),
            bar_plot(anc_size, "length of anchor seed", 50),
            bar_plot(anc_ambiguity, "ambiguity of anchor seed", 5),
            bar_plot(max_diag_deviations, "required nmw strip size to reach maximum possible score", 50, False),
            bar_plot(max_nmw_areas, "maximum area needleman wunsch", 50),
            ]))


def analyse_all_approaches_depre(out, db_name, query_size = 100, indel_size = 10):
    output_file(out)
    approaches = getApproachesWithData(db_name)
    plots = [ [], [], [], [], [], [], [], [] ]
    mapping_qual = []
    mapping_qual_illumina = []

    for approach_ in approaches:
        approach = approach_[0]
        results = getResults(db_name, approach, query_size, indel_size)
        hits = {}
        tries = {}
        run_times = {}
        scores = {}
        opt_scores = {}
        nums_seeds = {}
        nums_seeds_chosen = {}
        mapping_qual.append( ([],[]) )
        mapping_qual_illumina.append( ([],[]) )

        max_indel = 2*query_size/indel_size
        max_mut = query_size

        sub_illumina = 0.15
        indel_illumina = 0.05

        def init(d, x, y):
            if x not in d:
                d[x] = {}
            if y not in d[x]:
                d[x][y] = 0
            return d

        for score, score2, optimal_score, result_start, result_end, original_start, num_mutation, num_indels, num_seeds, num_seed_chosen_strip, mapping_quality, run_time in results:
            hits = init(hits, num_mutation, num_indels)
            tries = init(tries, num_mutation, num_indels)
            run_times = init(run_times, num_mutation, num_indels)
            if not score is None:
                scores = init(scores, num_mutation, num_indels)
            if not optimal_score is None:
                opt_scores = init(opt_scores, num_mutation, num_indels)
            nums_seeds = init(nums_seeds, num_mutation, num_indels)
            nums_seeds_chosen = init(nums_seeds_chosen, num_mutation, num_indels)

            if near(result_start, original_start, result_end, original_start+query_size):
                hits[num_mutation][num_indels] += 1
                if not optimal_score is None:
                    if optimal_score < score:
                        print("WARNING: aligner got better", score, "than optimal score", optimal_score)
                        opt_scores[num_mutation][num_indels] += 1
                    else:
                        opt_scores[num_mutation][num_indels] += score / optimal_score
                mapping_qual[-1][0].append(mapping_quality)
                if num_mutation/query_size < sub_illumina and num_indels/query_size < indel_illumina:
                    mapping_qual_illumina[-1][0].append(mapping_quality)
            else:
                mapping_qual[-1][1].append(mapping_quality)
                if num_mutation/query_size < sub_illumina and num_indels/query_size < indel_illumina:
                    mapping_qual_illumina[-1][1].append(mapping_quality)
            tries[num_mutation][num_indels] += 1
            run_times[num_mutation][num_indels] += run_time
            if not num_seeds is None:
                nums_seeds[num_mutation][num_indels] += num_seeds
            if not score is None:
                scores[num_mutation][num_indels] += score
            if not num_seed_chosen_strip is None:
                nums_seeds_chosen[num_mutation][num_indels] += num_seed_chosen_strip


        def makePicFromDict(d, w, h, divideBy, title, ignore_max_n = 0, log = False, set_max = None, set_min=None):
            pic = []
            min_ = 10000.0
            max_ = 0.0
            if len(d.keys()) == 0:
                return None
            else:
                w_keys = sorted(d.keys())
                h_keys = sorted(d[list(d.keys())[0]].keys())

                for x in w_keys:
                    pic.append( [] )
                    for y in h_keys:
                        if x not in d or y not in d[x] or divideBy[x][y] == 0:
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

            color_mapper = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 127),
                    low=min_,
                    high=max_
                )
            if log:
                color_mapper = LogColorMapper(
                        palette=heatmap_palette(light_spec_approximation, 127),
                        low=min_,
                        high=max_
                    )

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
        avg_runtime = makePicFromDict(run_times, max_mut, max_indel, tries, "runtime " + approach, 0)
        avg_score = makePicFromDict(scores, max_mut, max_indel, tries, "score " + approach)
        avg_opt_score = makePicFromDict(opt_scores, max_mut, max_indel, hits, "percent optimal score " + approach)
        avg_seeds = makePicFromDict(nums_seeds, max_mut, max_indel, tries, "num seeds " + approach)
        avg_seeds_ch = makePicFromDict(nums_seeds_chosen, max_mut, max_indel, tries, "num seeds in chosen SOC " + approach)
        seed_relevance = makePicFromDict(hits, max_mut, max_indel, nums_seeds, "seed relevance " + approach, set_max=0.01, set_min=0)

        if not avg_hits is None:
            plots[0].append(avg_hits)
        if not avg_runtime is None:
            plots[1].append(avg_runtime)
        if not avg_score is None:
            plots[2].append(avg_score)
        if not avg_seeds is None:
            plots[3].append(avg_seeds)
        if not avg_seeds_ch is None:
            plots[5].append(avg_seeds_ch)
        if not seed_relevance is None:
            plots[4].append(seed_relevance)
        if not avg_opt_score is None:
            plots[6].append(avg_opt_score)

    plot = figure(title="BWA-pic",
        x_axis_label='#wrong / #mapped', y_axis_label='#mapped / total',
        x_axis_type="log", y_range=(0,0.4),
        plot_width=650, plot_height=500,
        min_border_bottom=10, min_border_top=10,
        min_border_left=10, min_border_right=15
    )

    c_palette = heatmap_palette(light_spec_approximation, len(approaches))

    for index, approach_, in enumerate(approaches):
        approach = approach_[0]
        data = mapping_qual_illumina[index]

        total_amount_1 = len(data[0])
        total_amount_2 = len(data[1])
        all_data = []
        min_val = 1
        max_val = 0
        for x in data[0]:
            if x is None:
                continue
            all_data.append(x)
            if x < min_val:
                min_val = x
            if x > max_val:
                max_val = x
        for x in data[1]:
            if x is None:
                continue
            all_data.append(x)
            if x < min_val:
                min_val = x
            if x > max_val:
                max_val = x
                max_val = x

        values = []
        amount = 1000
        step = len(all_data) / amount
        c = 0
        for v in sorted(all_data):
            if c >= step:
                c = 0
                values.append(v)
            c += 1

        line_x = []
        line_y = []
        for val in values:
            mapped = 0
            wrong = 0
            total = total_amount_1 + total_amount_2
            for ele in data[0]:
                if not ele is None and ele > val:
                    mapped += 1
            for ele in data[1]:
                if not ele is None and ele > val:
                    mapped += 1
                    wrong += 1

            if mapped > 0:
                line_x.append( wrong/mapped )
                line_y.append( mapped/total )

        plot.line(line_x, line_y, legend=approach, color=c_palette[index], line_width=2)
        plot.x(line_x, line_y, legend=approach, color=c_palette[index], size=6)

    plot.legend.location = "top_left"
    plots[-1].append(plot)

    save(layout(plots))

def analyse_all_approaches(out, db_name, query_size = 100, indel_size = 10):
    output_file(out)
    approaches = getApproachesWithData(db_name)
    plots = [ [], [], [], [], [] ]
    mapping_qual = []
    yax = True
    for approach_ in approaches:
        approach = approach_[0]
        results = getResults(db_name, approach, query_size, indel_size)
        hits = {}
        tries = {}
        run_times = {}
        scores = {}
        nums_seeds = {}
        mapping_qual.append( ([],[]) )

        max_indel = 2*query_size/indel_size
        max_mut = query_size

        def init(d, x, y):
            if x not in d:
                d[x] = {}
            if y not in d[x]:
                d[x][y] = 0
            return d

        for score, score2, optimal_score, result_start, result_end, original_start, num_mutation, num_indels, num_seeds, num_seeds_chosen_strip, mapping_quality, run_time in results:
            hits = init(hits, num_mutation, num_indels)
            tries = init(tries, num_mutation, num_indels)
            run_times = init(run_times, num_mutation, num_indels)
            if not score is None:
                scores = init(scores, num_mutation, num_indels)
            nums_seeds = init(nums_seeds, num_mutation, num_indels)

            #print(result_start, original_start)
            if near(result_start, original_start, result_end, original_start+query_size):
                hits[num_mutation][num_indels] += 1
                mapping_qual[-1][0].append(mapping_quality)
            else:
                mapping_qual[-1][1].append(mapping_quality)
            tries[num_mutation][num_indels] += 1
            run_times[num_mutation][num_indels] += run_time
            if not num_seeds is None:
                nums_seeds[num_mutation][num_indels] += num_seeds
            if not score is None:
                scores[num_mutation][num_indels] += score


        def makePicFromDict(d, w, h, divideBy, title, ignore_max_n = 0, log = False, set_max = None, set_min=None, xax=True, yax=True):
            pic = []
            min_ = 10000.0
            max_ = 0.0
            if len(d.keys()) == 0:
                return None
            else:
                w_keys = sorted(d.keys())
                h_keys = sorted(d[list(d.keys())[0]].keys())

                for x in w_keys:
                    pic.append( [] )
                    for y in h_keys:
                        if x not in d or y not in d[x] or divideBy[x][y] == 0:
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

            color_mapper = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=min_,
                    high=max_
                )
            if log:
                color_mapper = LogColorMapper(
                        palette=heatmap_palette(light_spec_approximation, 256),
                        low=min_,
                        high=max_
                    )

            plot_width=500
            if yax:
                plot_width += 150
            plot_height=500
            if xax:
                plot_height += 40
            plot = figure(title=title,
                    x_range=(0,h), y_range=(0,w),
                    x_axis_label=str(indel_size) + 'nt indels', y_axis_label='mutations',
                    plot_width=plot_width, plot_height=plot_height,
                    min_border_bottom=10, min_border_top=10,
                    min_border_left=10, min_border_right=15,
                    tools=["save"]
                )
            #plot.xaxis.formatter = tick_formater
            plot.image(image=[pic], color_mapper=color_mapper,
                    dh=[w], dw=[h], x=[0], y=[0])

            font = "Helvetica"
            font_size = '15pt'
            if yax:
                color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
                color_bar.major_label_text_font=font
                color_bar.major_label_text_font_size=font_size
                plot.add_layout(color_bar, 'left')

            if not xax:
                plot.xaxis.visible = False
            if not yax:
                plot.yaxis.visible = False
            if not title is None:
                plot.title.text_font=font
                plot.title.text_font_size=font_size
            plot.legend.label_text_font=font
            plot.legend.label_text_font_size=font_size
            plot.axis.axis_label_text_font=font
            plot.axis.major_label_text_font=font
            plot.axis.axis_label_text_font_size=font_size
            plot.axis.major_label_text_font_size=font_size

            return plot


        avg_hits = makePicFromDict(hits, max_mut, max_indel, tries, approach, set_max=1, set_min=0, xax=False, yax=yax)
        avg_runtime = makePicFromDict(run_times, max_mut, max_indel, tries, None, 0, False, 0.01, 0, yax=yax)
        #avg_score = makePicFromDict(scores, max_mut, max_indel, tries, "score " + approach, yax=yax)
        avg_seeds = makePicFromDict(nums_seeds, max_mut, max_indel, tries, approach, yax=yax)
        avg_rel = makePicFromDict(hits, max_mut, max_indel, nums_seeds, approach, yax=yax, set_min=0, set_max=0.01)
        yax = False

        if not avg_hits is None:
            plots[0].append(avg_hits)
        if not avg_runtime is None:
            plots[1].append(avg_runtime)
        if not avg_seeds is None:
            plots[2].append(avg_seeds)
        if not avg_rel is None:
            plots[3].append(avg_rel)

    plot = figure(title="BWA-pic",
            x_axis_label='#wrong / #mapped', y_axis_label='#mapped / total',
            x_axis_type="log", y_range=(0,0.4),
            plot_width=650, plot_height=500,
            min_border_bottom=10, min_border_top=10,
            min_border_left=10, min_border_right=15
        )

    c_palette = heatmap_palette(light_spec_approximation, len(approaches))

    for index, approach_, in enumerate(approaches):
        approach = approach_[0]
        data = mapping_qual[index]

        total_amount_1 = len(data[0])
        total_amount_2 = len(data[1])
        all_data = []
        min_val = 1
        max_val = 0
        for x in data[0]:
            if x is None:
                continue
            all_data.append(x)
            if x < min_val:
                min_val = x
            if x > max_val:
                max_val = x
        for x in data[1]:
            if x is None:
                continue
            all_data.append(x)
            if x < min_val:
                min_val = x
            if x > max_val:
                max_val = x

        values = []
        amount = 1000
        step = len(all_data) / amount
        c = 0
        for v in sorted(all_data):
            if c >= step:
                c = 0
                values.append(v)
            c += 1

        line_x = []
        line_y = []
        for val in values:
            mapped = 0
            wrong = 0
            total = total_amount_1 + total_amount_2
            for ele in data[0]:
                if not ele is None and ele > val:
                    mapped += 1
            for ele in data[1]:
                if not ele is None and ele > val:
                    mapped += 1
                    wrong += 1

            if mapped > 0:
                line_x.append( wrong/mapped )
                line_y.append( mapped/total )

        plot.line(line_x, line_y, legend=approach, color=c_palette[index])
        plot.x(line_x, line_y, legend=approach, color=c_palette[index])

    plot.legend.location = "top_left"
    plots[-1].append(plot)

    save(gridplot(plots))



def compare_approaches(out, approaches, db_name, query_size = 100, indel_size = 10):
    output_file(out + ".html")
    plots = []

    max_indel = 2*query_size/indel_size
    max_mut = query_size

    def get_data(db_name, approach, query_size, indel_size):
        results = getResults(db_name, approach, query_size, indel_size)
        hits = {}
        tries = {}
        run_times = {}
        scores = {}
        nums_seeds = {}

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

        return (hits, tries, run_times, scores, nums_seeds)

    def convertToList(d, w, h, divideBy):
        pic = []
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
        return pic

    def sub(arr1, arr2):
        if arr1 is None or arr2 is None:
            return None
        for x, row in enumerate(arr2):
            for y, ele in enumerate(row):
                arr1[x][y] -= ele
        return arr1

    l_1 = get_data(db_name, approaches[0], query_size, indel_size)
    l_2 = get_data(db_name, approaches[1], query_size, indel_size)

    hits = sub(convertToList(l_1[0], max_indel, max_mut, l_1[1]),
            convertToList(l_2[0], max_indel, max_mut, l_2[1]))

    run_times = sub(convertToList(l_1[2], max_indel, max_mut, l_1[1]),
            convertToList(l_2[2], max_indel, max_mut, l_2[1]))

    scores = sub(convertToList(l_1[3], max_indel, max_mut, l_1[1]),
            convertToList(l_2[3], max_indel, max_mut, l_2[1]))

    nums_seeds = sub(convertToList(l_1[4], max_indel, max_mut, l_1[1]),
            convertToList(l_2[4], max_indel, max_mut, l_2[1]))


    def makePic(pic, w, h, title, log = False):
        if pic is None:
            return None
        def rgb(minimum, maximum, value):
            def clamp(x):
                return max(0, min(x, 255))
            value = min(max(minimum, value), maximum)
            red = green = blue = 0.0
            if value < 0:
                red = value/minimum
            else:
                green = value/maximum
            return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                                clamp(int(blue * 255)))

        def heatmap_palette(num_colors):
            """
            Color palette for visualization.
            """
            return [rgb(-1, 1, x) for x in np.linspace(-1, 1, num_colors)]

        max_ = -100000
        min_ =  100000

        for row in pic:
            for ele in row:
                if ele > max_:
                    max_ = ele
                if ele < min_:
                    min_ = ele

        if min_ < 0:
            min_ = -min_
        range_ = max_
        if min_ > range_:
            range_ = min_

        color_mapper = LinearColorMapper(palette=heatmap_palette(255), low=-range_, high=range_)
        if log:
            color_mapper = LogColorMapper(palette=heatmap_palette(255), low=-range_, high=range_)

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

    desc = approaches[0] + " - " + approaches[1]

    avg_hits = makePic(hits, max_mut, max_indel, "accuracy " + desc)
    avg_runtime = makePic(run_times, max_mut, max_indel, "runtime " + desc)
    avg_score = makePic(scores, max_mut, max_indel, "score " + desc)
    avg_seeds = makePic(nums_seeds, max_mut, max_indel, "num seeds " + desc)

    if not avg_hits is None:
        plots.append(avg_hits)
    if not avg_runtime is None:
        plots.append(avg_runtime)
    if not avg_score is None:
        plots.append(avg_score)
    if not avg_seeds is None:
        plots.append(avg_seeds)

    save(row(plots))

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


def get_ambiguity_distribution(reference, min_len=10, max_len=20):
    def get_all_queries(l):
        if l <= 0:
            yield ""
        else:
            for q in get_all_queries(l-1):
                yield q + "A"
                yield q + "C"
                yield q + "G"
                yield q + "T"
    def get_random_queries(l, amount):
        for _ in range(amount):
            q = ""
            for _ in range(l):
                char = random.randint(1,4)
                if char == 1:
                    q += "A"
                elif char == 2:
                    q += "C"
                elif char == 3:
                    q += "G"
                elif char == 4:
                    q += "T"
            yield q
    fm_index = FMIndex()
    fm_index.load(reference)
    num_queries = 100000

    r1max = 10
    r2size = 15

    indices1 = []
    indices2 = []
    for i in range(r1max):
        indices1.append(i)
    for i in range(r2size):
        indices2.append(2**i+r1max)
    data1 = []
    data2 = []
    for l in range(min_len, max_len):
        data1.append( [] )
        for _ in range(r1max):
            data1[-1].append(0.0)
        data2.append( [] )
        for _ in range(r2size):
            data2[-1].append(0.0)
        for q in get_random_queries(l, num_queries):
            ambiguity = fm_index.get_ambiguity(NucSeq(q))
            if ambiguity < r1max:
                data1[-1][int(ambiguity)] += 1.0/num_queries
            elif int(math.log2(ambiguity - r1max+1)) < r2size:
                data2[-1][int(math.log2(ambiguity - r1max+1))] += 1.0/num_queries

    print(indices1)
    print(data1)
    print(indices2)
    print(data2)

    color_mapper = LinearColorMapper(
                    palette=heatmap_palette(light_spec_approximation, 256),
                    low=0,
                    high=1
                )

    plot = figure(title="ambiguity on human genome",
            x_range=(0,r1max), y_range=(min_len, max_len),
            x_axis_label='ambiguity', y_axis_label='sequence length',
            plot_width=700, plot_height=500,
            min_border_bottom=10, min_border_top=10,
            min_border_left=10, min_border_right=15,
            tools=["save"]
        )
    plot.image(image=[data1], color_mapper=color_mapper,
            dh=[max_len - min_len], dw=[r1max], x=[0], y=[min_len])

    plot2 = figure(x_range=(r1max,2**r2size+r1max), y_range=(min_len, max_len),
            min_border_bottom=10, min_border_top=10,
            min_border_left=20, min_border_right=15,
            plot_width=500, plot_height=500,tools=[],
            x_axis_type="log"
        )
    plot2.image(image=[data2], color_mapper=color_mapper,
            dh=[max_len - min_len], dw=[2**r2size+r1max], x=[r1max], y=[min_len])

    font = "Helvetica"
    font_size = '15pt'
    color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
    color_bar.major_label_text_font=font
    color_bar.major_label_text_font_size=font_size
    plot.add_layout(color_bar, 'left')

    plot.legend.label_text_font=font
    plot.legend.label_text_font_size=font_size
    plot.axis.axis_label_text_font=font
    plot.axis.major_label_text_font=font
    plot.axis.axis_label_text_font_size=font_size
    plot.axis.major_label_text_font_size=font_size

    plot2.legend.label_text_font=font
    plot2.legend.label_text_font_size=font_size
    plot2.yaxis.visible = False
    plot2.axis.axis_label_text_font=font
    plot2.axis.major_label_text_font=font
    plot2.axis.axis_label_text_font_size=font_size
    plot2.axis.major_label_text_font_size=font_size

    show(gridplot( [[plot, plot2]] ))

"""
f = FileReader("test.fasta")
nuc = f.execute()
print(nuc.fastaq())
nuc = f.execute()
print(nuc.fastaq())

exit()
"""
#memory_test(human_genome, 1)
#get_ambiguity_distribution(human_genome)
#exit()

#createSampleQueries(human_genome, "/mnt/ssd1/shortIndels.db", 1000, 3, 128, True)
#test_my_approaches("/mnt/ssd1/shortIndels.db")
#analyse_all_approaches("shortIndels.html","/mnt/ssd1/shortIndels.db", 1000, 3)

#createSampleQueries(human_genome, "/mnt/ssd1/illumina.db", 50, 5, 64, True, True)
#test_my_approaches("/mnt/ssd1/illumina.db")
#analyse_all_approaches("illumina.html","/mnt/ssd1/illumina.db", 50, 5)
#analyse_all_approaches_depre("illumina_depre.html","/mnt/ssd1/illumina.db", 50, 5)


#high quality picture

#createSampleQueries(human_genome, "/mnt/ssd1/highQual.db", 1000, 100, 32, True, True)
test_my_approaches("/mnt/ssd1/highQual.db")
analyse_all_approaches("highQual.html","/mnt/ssd1/highQual.db", 1000, 100)
analyse_all_approaches_depre("highQual_depre.html","/mnt/ssd1/highQual.db", 1000, 100)
analyse_detailed("stats/", "/mnt/ssd1/highQual.db")
exit()

#createSampleQueries(human_genome, "/mnt/ssd1/veryHighQual.db", 1000, 100, 2**13, True, True)
#createSampleQueries(human_genome, "/mnt/ssd1/veryHighQual.db", 1000, 100, 2**7, True, True)
#createSampleQueries(human_genome, "/mnt/ssd1/veryHighQual.db", 1000, 100, 0, True, True)
test_my_approaches("/mnt/ssd1/veryHighQual.db")
analyse_all_approaches("highQual.html","/mnt/ssd1/veryHighQual.db", 1000, 100)
analyse_all_approaches_depre("highQual_depre.html","/mnt/ssd1/veryHighQual.db", 1000, 100)

#createSampleQueries(human_genome, "/mnt/ssd1/seedRelevance.db", 1000, 100, 2**7, True, True)
#test_my_approaches_rele("/mnt/ssd1/seedRelevance.db")
#analyse_all_approaches("seedRelevance.html","/mnt/ssd1/seedRelevance.db", 1000, 100)
#analyse_all_approaches_depre("seedRelevance_depre.html","/mnt/ssd1/seedRelevance.db", 1000, 100)
#analyse_detailed("seedRelevance/", "/mnt/ssd1/seedRelevance.db")

#createSampleQueries(human_genome, "/mnt/ssd1/stats.db", 1000, 100, 32, True, True)
#test_my_approaches("/mnt/ssd1/stats.db")
#analyse_all_approaches("stats.html","/mnt/ssd1/stats.db", 1000, 100)
#analyse_all_approaches_depre("stats_depre.html","/mnt/ssd1/stats.db", 1000, 100)
#analyse_detailed("stats/", "/mnt/ssd1/stats.db")
#filter_suggestion("/mnt/ssd1/stats.db", min_percentage_cov=0.02, max_num_seeds=4000, min_num_seed_in_strip=2)


#compare_approaches("results_comp_me_bwa", ["pBs:5 SoC:100,500nt sLs:1", "bwa"], db_name, 100, 10)
#compare_approaches("results_comp_me_me", ["pBs:5 SoC:100,500nt sLs:1", "pBs:5 SoC:1,500nt sLs:inf"], db_name, 1000, 100)
