from MA import *
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
import random
import gc
import os
import math
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool
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
from measure_time import *
from scipy import stats
from sklearn import linear_model
from subprocess import call
import json
import glob

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

    return (i*r**m, i*g**m, i*b**m)

def heatmap_palette(scheme, num_colors):
    def format(rgb):
        def clamp(x):
            return max(0, min(x, 255))
        red, green, blue = rgb
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                               clamp(int(blue * 255)))
    return [format(scheme(x)) for x in np.linspace(0, 1, num_colors)]

human_genome = "/MAdata/genome/human"
human_hg37_genome = "/MAdata/genome/human_hg37"
human_hg37_genome_full = "/MAdata/genome/GRCh37_full"
zebrafish_genome = "/MAdata/genome/zebrafish"
random_genome = "/MAdata/genome/random"
mouse_genome = "/MAdata/genome/mouse"
plasmodium_genome = "/MAdata/genome/plasmodium"
zebrafish_genome_n = "/MAdata/genome/zebrafish_n"
e_coli_genome = "/MAdata/genome/eColi"

## @brief Yield successive n-sized chunks from l.
def chunks(l, n):
    if n == 0:
        yield l
    elif n >= len(l):
        yield l
    else:
        for i in range(0, len(l), n):
            yield l[i:i + n]


def print_non_actg_chars(file_name):
    with open(file_name, 'r') as f:
        for line_num, line in enumerate(f):
            if line[0] == ">":
                print("sequence", line[:-1])
            else:
                for pos, char in enumerate(line):
                    if not char in ["A", "a", "C", "c", "G", "g", "T", "t", "\n"]:
                        print("char", char, "at position", pos, "in line", line_num)

def first_accurate_SoC(db_name, reference, output_file, max_span_seed_set=True):
    seeding = BinarySeeding(max_span_seed_set)
    soc = StripOfConsideration(0)

    ref_pack = Pack()
    ref_pack.load(reference)
    fm_index = FMIndex()
    fm_index.load(reference)
    
    queries = getQueriesAsASDMatrix("/MAdata/db/"+db_name)
    asd = []
    for i, row in enumerate(queries):
        asd.append( [] )
        print(i/len(queries))
        for cell in row:
            asd[-1].append( [] )
            for sequence, sample_id in cell:
                optima = getOptima("/MAdata/db/"+db_name, sample_id)
                origin = getOrigin("/MAdata/db/"+db_name, sample_id)

                pos = [ (origin[0], origin[0]+len(sequence)) ]
                pos.extend(optima)

                index = 0

                query = NucSeq(sequence)

                seed_set = seeding.execute(fm_index, query)

                soc_queue = soc.execute(seed_set, query, ref_pack, fm_index)


                found_it = False
                while not found_it:
                    if soc_queue.empty():
                        index = float("inf")
                        break
                    soc_instance = soc_queue.pop()
                    for seed in soc_instance:
                        s_start = seed.start_ref
                        s_end = seed.start_ref + seed.size

                        if s_start >= ref_pack.unpacked_size_single_strand:
                            s_start = ref_pack.unpacked_size() - s_end
                            s_end = seed.start_ref + seed.size

                        for start, end in pos:
                            if near(s_start, start, s_end, end) > 0:
                                found_it = True
                                break
                        if found_it:
                            break
                    index += 1

                asd[-1][-1].append(index)

    print("median:")
    for row in asd:
        for cell in row:
            print(median(cell), end="\t")
        print()
    print("max (of non inf values):")
    for row in asd:
        for cell in row:
            print(max([x for x in cell if x != float('inf')]), end="\t")
        print()
    print("min:")
    for row in asd:
        for cell in row:
            print(min(cell), end="\t")
        print()

    with open(output_file + ".json", "w") as f:
        json.dump([db_name, reference, max_span_seed_set, asd], f)
    print("wrote file")

#creating samples int the database
#createSampleQueries(working_genome, db_name, 1000, 100, 256)

def test_my_approach(
        db_name,
        reference,
        name,
        max_hits=1000,
        complete_seeds=False,
        use_chaining=False,
        local=True,
        reseed=False,
        full_analysis=False,
        do_minimizers=False,
        sort_after_score=True,
        max_nmw = 0,
        cheat=False,
        kMerExtension=False,
        reportN=0,
        clean_up_db=False,
        toGlobal=True,
        give_up= 0, #0.002, 16.05.18
        quitet=False,
        missed_alignments_db=None,
        min_coverage=1.1, #0.5, # 1.1 = force SW alignment SETTING DISABLED
        #optimistic_gap_estimation=False,
        specific_id=None,
        scatter_plot=False,
        specific_query=None,
        be_mean=False, # replace 10% of all symbols with N's
        analyze_heuristics=False
    ):
    if not quitet:
        print("collecting samples (" + name + ") ...")
    
    if not specific_query is None:
        db_name = None

    assert(db_name == None or specific_query == None)
    assert(not (db_name == None and specific_query == None))
    all_queries = None
    if specific_query is None: 
        all_queries = getQueries(db_name, specific_id)
        if be_mean:
            print("WARNING: i was told to be mean and replace 10% of all nucs with 'N'")
            for index in range(0, len(all_queries)):
                sequence, sample_id = all_queries[index]
                for _ in range(int(len(sequence)/10)):
                    x = random.randint(0, len(sequence))
                    sequence = sequence[:x] + "N" + sequence[x+1:]
                all_queries[index] = (sequence, sample_id)

    else:
        all_queries = [ (specific_query, 0) ] # sequence, sample_id (which is a dummy...)

    if not missed_alignments_db is None:
        conn = sqlite3.connect(missed_alignments_db)
        setUpDbTables(conn)

    if not db_name is None:
        clearApproach(db_name, name)

    if not quitet:
        print("having ", len(all_queries), " samples total (", name, ") ...")
        if len(all_queries) == 0:
            print("no new queries ... done")
            return

    runtimes = {
        'Seeding': 0,
        'SoC': 0,
        'SoC-Extraction': 0,
        'SoC-Sorting': 0,
        'SoC-Linesweep': 0,
        'Harm': 0,
        'DP': 0,
        'Map Qual': 0
    }# dict
    collect_ids = []

    num_covered_irelevant_seeds = 0
    num_irelevant_seeds = 0
    num_seeds_total = 0
    num_soc_seeds = 0
    num_coupled_seeds = 0
    warn_once = True
    picked_wrong_count = 0

    extract_size = 2**15
    # break samples into chunks of 2^15
    for index, queries, in enumerate(chunks(all_queries, extract_size)):
        if not quitet:
            if extract_size != 0:
                print("extracting", len(queries), "samples", index, "/",
                    len(all_queries)/extract_size, "(", name, ") ...")

            print("setting up (", name, ") ...")

        ref_pack = Pack()
        ref_pack.load(reference)
        ref_pledge = Pledge(Pack())
        ref_pledge.set(ref_pack)

        #ref_pledge.get().printHoles()
        #exit()

        fm_pledge = None

        if do_minimizers:
            minimizersHash = MinimizersHash().from_file(reference + ".maRef")
            minimizersHash.verify(ref_pack)# self tests...
            fm_pledge = Pledge(MinimizersHash())
            fm_pledge.set(minimizersHash)
        else:
            fm_index = FMIndex()
            fm_index.load(reference)
            fm_pledge = Pledge(FMIndex())
            fm_pledge.set(fm_index)

        missed_list = []

        # DP scoring scheme
        NeedlemanWunsch().max_gap_area = 10000 # if local else 0 #1000000
        #modules
        seeding = BinarySeeding(not complete_seeds)

        seeding.min_seed_size_drop = 0

        if kMerExtension:
            seeding = OtherSeeding(True)
        soc = StripOfConsideration(give_up) # check if 0 is okay
        ex = ExtractAllSeeds(max_hits, 0)
        ls = LinearLineSweep()

        ls.equal_score_lookahead = 5
        ls.tolerance = 0

        couple = ExecOnVec(ls, True, max_nmw)
        nmw = NeedlemanWunsch()
        optimal = ExecOnVec(nmw, sort_after_score)
        mappingQual = MappingQuality(max_nmw) # give me max_nmw alignments

        pledges = [[], [], [], [], [], []]

        for sequence, sample_id in queries:
            pledges[0].append(Pledge(NucSeq()))
            pledges[0][-1].set(NucSeq(sequence))
            pledges[0][-1].get().name = str(sample_id)
            pledges[1].append(seeding.promise_me(
                    fm_pledge, pledges[0][-1]
                ))
            pledges[2].append(soc.promise_me(
                    pledges[1][-1], pledges[0][-1], ref_pledge, fm_pledge
                ))
            pledges[3].append(ls.promise_me(
                    pledges[2][-1], pledges[0][-1]
                ))
            pledges[4].append(optimal.promise_me(
                    pledges[3][-1], pledges[0][-1], ref_pledge
                ))
            pledges[5].append(mappingQual.promise_me(
                pledges[0][-1], pledges[4][-1]
            ))

            optimal_alignment_in = []
            optimal_alignment_in.append((Pledge(NucSeq()), Pledge(NucSeq())))

        smw = SMW(full_analysis)
        optimal_alignment_out = []
        for a, b in optimal_alignment_in:
            optimal_alignment_out.append(smw.promise_me(a, b))
        
        if not quitet:
            ## print the computational graph description
            print("computational graphs: ")
            print(pledges[-1][0].get_graph_desc())
            #exit()

        if not quitet:
            print("computing (", name, ") ...")
        Pledge.simultaneous_get(pledges[-1], 32)
        print("")

        if scatter_plot:
            print("generating scatterplot...")
            plot1 = figure()
            for alignments in pledges[-1]:
                for alignment in alignments.get():
                    for scatter in alignment.vGapsScatter:
                        plot1.x(scatter.first, scatter.second, color="black", size=4)
            print("done")
            show(plot1)

        if (specific_id != None and db_name != None) or specific_query != None:
            plot1 = figure(plot_width=1200)
            x_list = []
            y_list = []
            print("Num seeds: ", pledges[1][0].get().num_seeds(0))
            print("Num SoC: ", len(pledges[2][0].get()))
            for s in pledges[2][0].get().scores:
                y_list.append(s.first)
                x_list.append(s.second)
            x_list2 = []
            y_list2 = []
            w_list2 = []
            h_list2 = []
            optima = []
            if db_name != None:
                for tup in getOptima(db_name, specific_id):
                    optima.append(tup[0])
                    optima.append(ref_pack.unpacked_size_single_strand*2 - tup[0])
                tup = getOrigin(db_name, specific_id)
                optima.append(tup[0])
                optima.append(ref_pack.unpacked_size_single_strand*2 - tup[0])
            for x in optima:
                x_list2.append(x)
                w_list2.append(2000)
                y_list2.append(100)
                h_list2.append(400)
            plot1.rect(x_list2, y_list2, w_list2, h_list2, color="green")
            plot1.line(x_list, y_list, line_width=2)
            plot1.x(x_list, y_list, color="red")
            """
            SoC scores in order
            """
            plot2 = figure(plot_width=1200)
            x_list = []
            y_list1 = []
            y_list2 = []
            y_list3 = []
            y_list4 = []
            x_list_harm_areas_start = []
            x_list_harm_areas_end = []
            x_list_soc_areas_start = []
            x_list_soc_areas_end = []
            ex = pledges[2][0].get().extract
            for alignment in pledges[-1][0].get():
                x = alignment.stats.index_of_strip
                s = ex[x]
                x_list.append(x)
                y_list1.append(s.first)
                y_list2.append(s.second)
                y_list3.append(alignment.get_score())
                y_list4.append(s.qCoverage)
                x_list_harm_areas_start.append( (s.rStart + s.rEnd) / 2 )
                x_list_harm_areas_end.append(s.rEnd - s.rStart)
                x_list_soc_areas_start.append( (s.rStartSoC + s.rEndSoC) / 2 )
                x_list_soc_areas_end.append(s.rEndSoC - s.rStartSoC)
            plot1.rect(x_list_harm_areas_start, 0, x_list_harm_areas_end, 400, color="magenta")
            plot1.rect(x_list_soc_areas_start, 0, x_list_soc_areas_end, 200, color="pink")
            plot2.line(x_list, y_list4, color="purple", legend="coverage", line_width=5)
            plot2.x(x_list, y_list4, color="purple", size=5)
            plot2.line(x_list, y_list1, color="red", legend="SoC score", line_width=5)
            plot2.x(x_list, y_list1, color="red", size=5)
            plot2.line(x_list, y_list2, color="green", legend="Harm score", line_width=2)
            plot2.x(x_list, y_list2, color="green", size=5)
            plot2.line(x_list, y_list3, color="blue", legend="SW score", line_width=1)
            plot2.x(x_list, y_list3, color="blue", size=5)
            plots = [[plot1], [plot2]]
            max_socs = 0
            for num, soc in enumerate(pledges[2][0].get().vSoCs):
                plot3 = figure(plot_width=1200, name="SoC #" + str(num))
                if max_socs <= 0:
                    break
                max_socs -= 1
                p_counter = 0
                for seed in soc:
                    plot3.line(
                            [seed.start_ref, seed.start_ref+seed.size],
                            [seed.start, seed.start+seed.size],
                            color="blue",
                            line_width=3
                        )
                    plot3.x(
                            [seed.start_ref, seed.start_ref+seed.size],
                            [seed.start, seed.start+seed.size],
                            color="black"
                        )
                    p_counter += 1
                    if p_counter % 100 == 0:
                        print(p_counter)

                plots.append([plot3])

            show(layout(plots))

        total_time = 0
        for pledge in pledges[3]:
            total_time += pledge.exec_time
        print("total Time: ", total_time)

        if db_name != None:
            putTotalRuntime(db_name, name, total_time)

        if specific_id != None:
            if len(pledges[-1][0].get()) == 0:
                print("No alignment Done")
            else:
                for alignment in pledges[-1][0].get():
                    AlignmentPrinter().execute(alignment, pledges[0][0].get(), ref_pack)
                    plot = figure(width=800)
                    x = []
                    y = []
                    max_y = 0
                    for seed in pledges[2][0].get().vSoCs[alignment.stats.index_of_strip]:
                        plot.line([seed.start_ref, seed.start_ref + seed.size], [seed.start, seed.start + seed.size], line_width=6)
                        x.append( seed.start_ref + seed.size/2 )
                        y.append( seed.start + seed.size/2 )
                        if seed.start > max_y:
                            max_y = seed.start
                    for seed in pledges[2][0].get().vIngroup[alignment.stats.index_of_strip]:
                        plot.x(seed.start_ref, seed.start, color="green", size=10)
                    for seed in pledges[2][0].get().vHarmSoCs[alignment.stats.index_of_strip]:
                        plot.line([seed.start_ref, seed.start_ref + seed.size], [seed.start, seed.start + seed.size], color="red", line_width=5)

                    slope = pledges[2][0].get().vSlopes[alignment.stats.index_of_strip]
                    intercept = pledges[2][0].get().vIntercepts[alignment.stats.index_of_strip]
                    print(slope, intercept)

                    plot.line([intercept, intercept+max_y], [0, max_y*slope], color="green", line_width=4, line_alpha=0.5)

                    try:
                        X = np.array(x)[:, np.newaxis]
                        y = np.array(y)
                        ransac = linear_model.RANSACRegressor()
                        ransac.fit(X, y)
                        print(alignment.stats.index_of_strip, "py residual_threshold:", np.median(np.abs(y - np.median(y))))
                        plot.cross(X[ransac.inlier_mask_].flatten(), y[ransac.inlier_mask_], color="magenta", size=10)
                        line_x = np.arange(X.min(), X.max())
                        line_X = line_x[:, np.newaxis]
                        line_y = ransac.predict(line_X)
                        slope2, intercept2, r_value2, p_value2, error2 = stats.linregress(X[ransac.inlier_mask_].flatten(), y[ransac.inlier_mask_])
                        print("python regression ERROR:", error2)
                        plot.line(line_x, line_y, color="magenta", line_width=4, line_alpha=0.5)
                    except:
                        print("Python ransac failed")
                        pass

                    show(plot)


        # @note temporary debug code end

        if full_analysis:
            if not quitet:
                print("setting up optimal (", name, ") ...")
            for i, query_ in enumerate(queries):
                alignment = pledges[-1][i].get()[0]
                if alignment.begin_on_query == alignment.end_on_query or not full_analysis:
                    optimal_alignment_in[i][0].set(NucSeq())
                    optimal_alignment_in[i][1].set(NucSeq())
                    continue
                query = query_[0]
                optimal_alignment_in[i][0].set(
                    NucSeq(query[alignment.begin_on_query : alignment.end_on_query])
                )
                optimal_alignment_in[i][1].set(
                    ref_pack.extract_from_to(alignment.begin_on_ref, alignment.end_on_ref)
                )
            if not quitet:
                print("computing optimal (", name, ") ...")
            Pledge.simultaneous_get(optimal_alignment_out, 32)


        if not quitet:
            print("extracting results (", name, ") ...")
        result = []
        runtimes_result = []
        for i, alignments in enumerate(pledges[-1]):
            if len(alignments.get()) == 0:
                continue

            alignment = alignments.get()[0]
            #print(alignment.begin_on_ref)


            #check cigar for errors and complain if there are any
            #AlignmentPrinter(check_for_errors_only=True).execute(alignment, pledges[0][i].get(), ref_pack)

            if len(alignment) == 0:
                continue
            if cheat:
                ##
                # pretend backend is perfect...
                if not near(
                        alignment.begin_on_ref,
                        origin_pos,
                        alignment.end_on_ref,
                        origin_pos+orig_size):
                    for soc in pledges[2][i].get():
                        if len(soc) == 0:
                            continue
                        begin_on_ref = soc[0].start_ref
                        end_on_ref = begin_on_ref + soc[0].size
                        for seed in soc:
                            if begin_on_ref > seed.start_ref:
                                begin_on_ref = seed.start_ref
                            if end_on_ref < seed.start_ref + seed.size:
                                end_on_ref = seed.start_ref + seed.size
                        _, _, origin_pos, orig_size = queries[i]
                        x_hit = near(
                            begin_on_ref,
                            origin_pos,
                            end_on_ref,
                            origin_pos+orig_size)
                        if x_hit:
                            alignment.begin_on_ref = begin_on_ref
                            alignment.end_on_ref = end_on_ref
                ##
                # pretend backend is perfect end
            alignment2 = None
            optimal_alignment = None
            if full_analysis:
                if len(alignments.get()) > 1:
                    _, _, origin_pos, orig_size = queries[i]
                    align_1_hit = near(
                        alignment.begin_on_ref,
                        origin_pos,
                        alignment.end_on_ref,
                        origin_pos+orig_size)
                    if not align_1_hit:
                        print("Warning picked wrong alignment")
                        print("query id: ", i)
                    found_correct_alignment = align_1_hit
                    for index_, alignment2 in enumerate(alignments.get()[1:]):
                        index = index_ + 1
                        #
                        # check if the second alignment was accurate and the first was not...
                        #
                        align_2_hit = near(
                            alignment2.begin_on_ref,
                            origin_pos,
                            alignment2.end_on_ref,
                            origin_pos+orig_size)
                        if align_2_hit and not align_1_hit:
                            found_correct_alignment = True
                            print("correct was:", index)
                            def get_seed_coverage(alignment):
                                return (
                                    "Query: " + str(alignment.stats.initial_q_beg) + " - " +
                                    str(alignment.stats.initial_q_end) + " Reference: " +
                                    str(alignment.stats.initial_r_beg) + " - " +
                                    str(alignment.stats.initial_r_end)
                                    )
                            def print_seeds(alignment):
                                print("Num Seeds: ", alignment.stats.num_seeds_in_strip,
                                    " SoC index:", alignment.stats.index_of_strip)
                                for seed in pledges[2][i].get()[alignment.stats.index_of_strip]:
                                    print( (seed.start, seed.start_ref, seed.size) )
                            print("scores are:",alignment.get_score(),"?>=", alignment2.get_score())
                            print("Alignment 0 is hit: ", align_1_hit)
                            print("Alignment 0 seed coverage: ", get_seed_coverage(alignment))
                            print("Alignment 0 seeds: ")
                            print_seeds(alignment)
                            AlignmentPrinter().execute(alignment, pledges[0][i].get(), ref_pack)
                            print("Alignment", index, "is hit: ", align_2_hit)
                            print("Alignment", index, "seed coverage: ", get_seed_coverage(alignment2))
                            print("Alignment", index, "seeds: ")
                            print_seeds(alignment2)
                            AlignmentPrinter().execute(alignment2, pledges[0][i].get(), ref_pack)
                            picked_wrong_count += 1
                        assert(not sort_after_score or alignment.get_local_score() >= alignment2.get_local_score())
                        #
                        # end of check
                        #
                    #
                    # another check
                    #
                    if not found_correct_alignment:
                        correct_soc = -1
                        for index, soc in enumerate(pledges[2][i].get()):
                            if len(soc) == 0:
                                continue
                            begin_on_ref = soc[0].start_ref
                            end_on_ref = begin_on_ref + soc[0].size
                            for seed in soc:
                                if begin_on_ref > seed.start_ref:
                                    begin_on_ref = seed.start_ref
                                if end_on_ref < seed.start_ref + seed.size:
                                    end_on_ref = seed.start_ref + seed.size
                            _, _, origin_pos, orig_size = queries[i]
                            x_hit = near(
                                begin_on_ref,
                                origin_pos,
                                end_on_ref,
                                origin_pos+orig_size)
                            if x_hit:
                                correct_soc = index
                                break
                        if not correct_soc == -1:
                            print("found correct SOC (index:", correct_soc, ") but no correct alignment")
                            for seed in soc:
                                print( (seed.start, seed.start_ref, seed.size) )
                            # now look for correct harmonized SoC
                            found_ham = False
                            for index2, harm in enumerate(pledges[3][i].get()):
                                if len(harm) == 0:
                                    continue
                                for seed in harm:
                                    for seed2 in soc:
                                        if seed.start == seed2.start and seed.start_ref == seed2.start_ref and seed.size == seed2.size:
                                            print("Found correct harmonized SoC", index2)
                                            for seed in harm:
                                                print( (seed.start, seed.start_ref, seed.size) )
                                            found_ham = True
                                            break
                                    if found_ham:
                                        break
                                if found_ham:
                                    break
                            if not found_ham:
                                print("Could not find respective harmonized SoC")
                                print("found", len(pledges[2][i].get()), "SoCs; after harmonization", len(pledges[3][i].get()), "remain")

                            found_soc_index = False
                            for index2, alignment in enumerate(alignments.get()):
                                if alignment.stats.index_of_strip == index:
                                    print("Incorrect alignment for correct SoC. Index:", index2)
                                    AlignmentPrinter().execute(alignment, pledges[0][i].get(), ref_pack)
                                    found_soc_index = True
                                    break
                            if not found_correct_alignment:
                                print("aligner discared correct SoC.")




                    #
                    # end of check
                    #
                if len(optimal_alignment_out[i].get()) > 0:
                    optimal_alignment = optimal_alignment_out[i].get()[0]
                    if optimal_alignment.end_on_query == 0 or optimal_alignment.end_on_ref == 0:
                        optimal_alignment = None
                    else:
                        # adjust the scores to the slice that 
                        # was actually computed with the optimal alignment
                        optimal_alignment.begin_on_query += alignment.begin_on_query
                        optimal_alignment.end_on_query += alignment.begin_on_query
                        optimal_alignment.begin_on_ref += alignment.begin_on_ref
                        optimal_alignment.end_on_ref += alignment.begin_on_ref
            sample_id = int(alignment.stats.name)


            
            # save the mean amount of seeds
            #runtimes_result.append(
            #    (
            #        sample_id,
            #        pledges[1][i].get().computeAccSize() , 
            #        # / len(pledges[0][i].get()) = did not work well
            #        name
            #    )
            #)

            collect_ids.append(sample_id)

            total_time = 0
            runtimes["Seeding"] += pledges[1][i].exec_time
            total_time += pledges[1][i].exec_time
            pledges[1][i].exec_time = 0
            runtimes["SoC"] += pledges[2][i].exec_time
            total_time += pledges[2][i].exec_time
            pledges[2][i].exec_time = 0
            
            runtimes["SoC-Extraction"] += pledges[1][i].get().fExtraction / 1000.0
            runtimes["SoC-Sorting"] += pledges[1][i].get().fSorting / 1000.0
            runtimes["SoC-Linesweep"] += pledges[1][i].get().fLinesweep / 1000.0

            runtimes["Harm"] += pledges[3][i].exec_time
            total_time += pledges[3][i].exec_time
            pledges[3][i].exec_time = 0
            runtimes["DP"] += pledges[4][i].exec_time
            total_time += pledges[4][i].exec_time
            pledges[4][i].exec_time = 0
            runtimes["Map Qual"] += pledges[5][i].exec_time
            total_time += pledges[5][i].exec_time
            pledges[5][i].exec_time = 0

            max_nmw_area = 0
            nmw_area = 0
            max_diag_deviation_percent = 0.0
            max_diag_deviation = 0
            curr_diag_deviation = 0
            curr_max_diag_deviation = 0
            gap_size = 0
            cur_nmw_h = 0
            cur_nmw_w = 0
            seed_coverage = 0.0

            def print_alignments():
                print("MA alignment:")
                AlignmentPrinter().execute(
                    alignment,
                    pledges[0][i].get(),
                    ref_pack
                )
                if optimal_alignment != None:
                    print("SA alignment:")
                    AlignmentPrinter().execute(
                        optimal_alignment,
                        pledges[0][i].get(),
                        ref_pack
                    )

            if optimal_alignment != None and alignment.get_score() > optimal_alignment.get_score():
                print("WARNING: alignment computed better than optimal score",
                      alignment.get_score(), optimal_alignment.get_score()
                     )
                print_alignments()
                sw = SMW(True)
                sw.print = True
                sw.execute(
                    pledges[0][i].get(),
                    ref_pack.extract_from_to(alignment.begin_on_ref, alignment.end_on_ref)
                )
            if (warn_once and local and optimal_alignment != None and
                    alignment.get_score() < optimal_alignment.get_score()):
                warn_once = False
                print("got worse than optimal score", alignment.get_score(),
                      optimal_alignment.get_score()
                     )
                print("this warning is just printed once")
                print_alignments()

            seed_coverage_soc = 0.0
            if full_analysis:
                for match_type in alignment.extract():
                    if match_type == MatchType.seed:
                        seed_coverage += 1.0
                        if (
                                gap_size > 100 and
                                nmw.penalty_missmatch * gap_size < nmw.penalty_gap_open +
                                nmw.penalty_gap_extend * gap_size
                            ):
                            if (
                                    float(abs(curr_max_diag_deviation)) / gap_size >
                                    max_diag_deviation_percent
                                ):
                                max_diag_deviation_percent = float(abs(curr_max_diag_deviation))
                                max_diag_deviation_percent /= gap_size
                        if abs(curr_max_diag_deviation) > max_diag_deviation:
                            max_diag_deviation = abs(curr_max_diag_deviation)
                        curr_diag_deviation = 0
                        curr_max_diag_deviation = 0
                        gap_size = 0
                        nmw_area += cur_nmw_h * cur_nmw_w
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

                seed_coverage /= len(pledges[0][i].get())

                #check for how many irrelevant overlapped seeds where there...
                #first get all relevant seeds:
                if not use_chaining:
                    discovered_seeds = pledges[1][i].get().extract_seeds(fm_index, max_hits, True)
                    soc_seeds = pledges[2][i].get()[alignment.stats.index_of_strip]
                    num_soc_seeds += len(soc_seeds)
                    #num_coupled_seeds += len(pledges[3][i].get()[alignment.stats.index_of_strip])
                    num_seeds_total += pledges[1][i].get().num_seeds(max_hits)
                    #compute the area covered by relevant seeds
                    covered_area = [-1]*len(pledges[0][i].get())
                    for seed in discovered_seeds:
                        for pos in range(seed.start, seed.start + seed.size):
                            if covered_area[pos] < seed.size:
                                covered_area[pos] = seed.size
                #seed coverage after the soc
                    covered_area_soc = [False]*len(pledges[0][i].get())
                    for seed in soc_seeds:
                        for pos in range(seed.start, seed.start + seed.size):
                            covered_area_soc[pos] = True
                    for cov in covered_area_soc:
                        if cov:
                            seed_coverage_soc += 1.0
                    seed_coverage_soc /= len(pledges[0][i].get())
                #run over all discovered seeds and count the covered irelevant ones
                for seed in discovered_seeds:
                    if not (seed.start_ref + seed.size >= queries[i][2] and
                            seed.start_ref <= queries[i][2] + queries[i][3]):
                        num_irelevant_seeds += 1
                        covered = True
                        for pos in range(seed.start, seed.start + seed.size):
                            if covered_area[pos] <= seed.size:
                                covered = False
                        if covered:
                            num_covered_irelevant_seeds += 1

            sc2 = None
            if not optimal_alignment is None:
                sc2 = optimal_alignment.get_score()
            score2 = 0
            if not alignment2 is None:
                score2 = alignment2.get_score()

            num_seeds = None
            if do_minimizers:
                num_seeds = len(pledges[1][i].get())
            else:
                num_seeds = pledges[1][i].get().num_seeds(max_hits)


            for alignment in alignments.get():
                forw = ref_pack.unpacked_size_single_strand
                if alignment.begin_on_ref > forw:
                    b = alignment.begin_on_ref
                    alignment.begin_on_ref = forw*2 - alignment.end_on_ref
                    alignment.end_on_ref = forw*2 - b
                result.append(
                    (
                        sample_id,
                        alignment.begin_on_ref,
                        alignment.end_on_ref,
                        alignment.mapping_quality,
                        name,
                        1 if alignment.secondary else 0,
                        pledges[1][i].get().num_seeds_larger(15)
                    )
                )
                runtimes_result.append(
                    (
                        sample_id,
                        total_time,
                        name
                    )
                )
        if not quitet:
            print("submitting results (", name, ") ...")
        if len(result) > 0 and db_name != None:
            submitResults(db_name, result)
        if len(runtimes_result) > 0 and db_name != None:
            submitRuntimesAsList(db_name, runtimes_result, name)

        # this should not be necessary but for some reason it is...
        # @todo: fix memory managemet in python code
        for pledge in pledges[-1]:
            pledge.clear_graph()

        print(missed_list)
        if not missed_alignments_db is None and len(missed_list) > 0:
            insertQueries(conn, missed_list)

    if not quitet:
        if num_seeds_total > 0 :
            print("collected", num_irelevant_seeds,
                "irrelevant seeds. Thats", num_irelevant_seeds/len(all_queries), "per alignment and",
                100*num_irelevant_seeds/num_seeds_total, "percent")
            print("collected", num_covered_irelevant_seeds,
                "covered irrelevant seeds. Thats", num_covered_irelevant_seeds/len(all_queries),
                "per alignment and", 100*num_covered_irelevant_seeds/num_seeds_total, "percent")
        print("collected", num_seeds_total-num_irelevant_seeds, "relevant seeds")
        if num_soc_seeds > 0:
            print("having", num_soc_seeds, "seeds in the strip of consideration, having",
                num_coupled_seeds, "seeds after coupling, thats",
                100*(1-num_coupled_seeds/num_soc_seeds), "percent seeds discarded")
        last = 0
        num_missed = 0
        for ele in sorted(collect_ids):
            if ele != last + 1:
                num_missed += 1
            last  = ele
        print("Missed", num_missed, "samples")
        print("Picked wrong SoC", picked_wrong_count, "times")
        print("total runtimes:")
        for key, value in sorted(runtimes.items()):
            print(value, "\t(", key, ")")
        print("done")


def relevance(db_name, working_genome):
    all_queries = getQueriesAsASDMatrix("/MAdata/db/" + db_name)

    fm_index = FMIndex()
    fm_index.load(working_genome)

    def analyse(seeding):
        result = []
        for row in all_queries:
            for cell in row:
                relevant = 0
                total = 0
                for sequence, sample_id in cell:
                    segments = seeding.execute(fm_index, NucSeq(sequence))
                    seeds = segments.extract_seeds(fm_index, 0, 0, True)
                    total += len(seeds)
                    optimas = []
                    optimas_ = getOptima("/MAdata/db/" + db_name, sample_id)
                    for start, end in optimas_:
                        optimas.append( (start, end - start) )
                    s, seq = getOrigin("/MAdata/db/" + db_name, sample_id)
                    optimas.append( (s, len(seq)) )
                    for seed in seeds:
                        for origin, original_size in optimas:
                            if seed.start_ref + seed.size >= origin and seed.start_ref <= origin + original_size:
                                relevant += 1
                                break
                result.append(relevant / total)
        return mean(result), median(result), np.std(result)

    print("max. spanning")
    print(analyse(BinarySeeding(True)))
    print("SMEMs")
    print(analyse(BinarySeeding(False)))
    print("16-mer")
    print(analyse(OtherSeeding(True)))

#relevance("human_30000.db", human_genome)
#exit()
"""
class MA_Parameter(CommandLine):
    def __init__(self, index_str, fast, db_name, num_soc, max_hits, min_ambiguity, match, give_up):
        super().__init__()
        self.ma_home = "/usr/home/markus/workspace/aligner/"
        self.index_str = index_str
        self.threads = 32
        self.num_results = 3
        self.fast = "accurate"
        self.num_soc = num_soc
        self.max_hits = max_hits
        self.min_ambiguity = min_ambiguity
        self.match = match
        self.give_up = give_up
        if fast:
            self.fast = "fast"
        self.in_filename = ".tempMA" + self.fast + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.ma_home + "ma -a -t " + str(self.threads) + " -p " + self.fast
        cmd_str += " -S " + str(self.num_soc)
        cmd_str += " -A " + str(self.max_hits)
        cmd_str += " -B " + str(self.min_ambiguity)
        cmd_str += " --Match " + str(self.match)
        cmd_str += " -v " + str(self.give_up)
        return cmd_str + " -g " + self.index_str + " -i " + in_filename + " -n " + self.num_results

    def do_checks(self):
        return True
"""

def try_out_parameters(db_name, working_genome):
    for query_size, indel_size in [ (1000, 100), (200, 20), (30000, 100) ]:
        print("query_size", query_size, "indel_size", indel_size)
        createSampleQueries(working_genome, db_name, query_size, indel_size, 32, high_qual=False, quiet=True)
        for match in [1, 2, 3, 4, 5, 6]: # 3 seems to be best
            for num_soc in [3, 5, 30]:
                for min_ambiguity in [0, 1, 2, 3]:
                    for give_up in [.01, .02, .03, .05, .07, .08]: # seems to not matter
                        suffix = str(num_soc) + "_" + str(min_ambiguity) + "_" + str(give_up) + "_" + str(match)
                        test_my_approach(db_name, working_genome, "Fast_" + suffix, num_strips=num_soc, complete_seeds=False, full_analysis=False, local=True, max_nmw=3, min_ambiguity=min_ambiguity, give_up=give_up, match=match, quitet=True)

                        test_my_approach(db_name, working_genome, "Acc_" + suffix, num_strips=num_soc, complete_seeds=True, full_analysis=False, local=False, max_nmw=10, min_ambiguity=min_ambiguity, give_up=give_up, match=match, quitet=True)

        approaches = getApproachesWithData(db_name)
        final_result = []
        for approach_ in approaches:
            approach = approach_[0]
            results = getResults(db_name, approach, query_size, indel_size)

            ##
            # picks the correct alignment if multiple alignments are reported for one query...
            #
            def amunt_hits(results, approach):
                ret = []
                by_sample_id = {}
                runtime = 0
                for result in results:
                    runtime += result[10]
                    sample_id = result[-2]
                    if not sample_id in by_sample_id:
                        by_sample_id[sample_id] = []
                    by_sample_id[sample_id].append(result)

                hits = 0
                hits_first_c = 0
                for result_list in by_sample_id.values():
                    for result in result_list:
                        # get the interesting parts of the tuple...
                        _, _, _, result_start, result_end, original_start, _, num_indels, _, _, _, _, _, _, _, _, _, secondary_ = result
                        if near(result_start, original_start, result_end, original_start+query_size):
                            if num_indels == 0:
                                hits_first_c += 1
                            hits += 1
                            break
                return (hits, approach, runtime, hits_first_c)

            final_result.append( amunt_hits(results, approach) )

        #sort by second item and print in order
        print("F/A_numSoc_minAmbiguity_giveUp_matchScore")
        print("approach", "hits", "hits first column", "runtime", sep="\t")
        for hits, approach, runtime, hits_first_c in sorted(final_result, key=lambda x: x[0]):
            print(approach, hits, hits_first_c, runtime, sep="\t")

def test_my_approaches(db_name, genome, missed_alignments_db=None, specific_id=None, specific_query=None, be_mean=False, analyze_heuristics=False):
    full_analysis = False

    test_my_approach("/MAdata/db/"+db_name, genome, "MA Accurate PY", complete_seeds=True, full_analysis=full_analysis, local=False, specific_id=specific_id, specific_query=specific_query, be_mean=be_mean, max_hits=100, analyze_heuristics=analyze_heuristics)

    test_my_approach("/MAdata/db/"+db_name, genome, "MA Fast PY", complete_seeds=False, full_analysis=full_analysis, local=False, specific_id=specific_id, specific_query=specific_query, be_mean=be_mean, analyze_heuristics=analyze_heuristics)

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
        max_diag_deviations_percent = ([], [])
        max_nmw_areas = ([], [])
        nmw_areas = ([], [])
        mapping_qual = ([], [])
        score_dif = ([], [])

        approach = approach_[0]
        for result in getResults(db_name, approach):
            score, score2, optimal_score_this_region, result_start, result_end, original_start, num_mutation, num_indels, num_seeds, mapping_quality, run_time, index_of_chosen_strip, seed_coverage_chosen_strip, seed_coverage_alignment, num_seeds_chosen_strip, anchor_size, anchor_ambiguity, max_diag_deviation, max_diag_deviation_percent, max_nmw_area, nmw_area, original_size = result

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
            max_diag_deviations_percent[index].append(max_diag_deviation_percent)
            max_nmw_areas[index].append(max_nmw_area)
            nmw_areas[index].append(nmw_area)

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
            bar_plot(max_diag_deviations, "required nmw strip size to reach maximum possible score", 50),
            bar_plot(max_diag_deviations_percent, "required nmw strip size to reach maximum possible score", 50, False),
            bar_plot(max_nmw_areas, "maximum area needleman wunsch", 50),
            bar_plot(nmw_areas, "maximum area needleman wunsch", 50),
            ]))

##
# checks weather the database entries of two approaches are equal
#
def expecting_same_results(a, b, db_name, query_size = 100, indel_size = 10):
    results1 = getResults(db_name, a, query_size, indel_size)
    results2 = getResults(db_name, b, query_size, indel_size)

    for index in range(len(results1)):
        tup1 = results1[index]
        tup2 = results2[index]
        if tup1 is None:
            print("tup1 is none")
            print(tup2)
            return
        if tup2 is None:
            print("tup2 is none")
            print(tup1)
            return
        if tup1[3] != tup2[3]:#result_start
            print("Having unequal db entries")
            print(tup1)
            print(tup2)
            return
    print("Having equal db entries")

def analyse_all_approaches_depre(out, db_name, num_tries=1, print_relevance=False, allow_sw_hits=True, show_coverage=True):
    db_name = "/MAdata/db/" + db_name
    output_file(out)
    plots = [ [], [], [] ]
    json_save = [ ]

    approaches = getApproachesWithData(db_name)

    def makePicFromDict(d, title, print_out=False,desc2=None, set_max=None, inner=None):
        x_hover = []
        y_hover = []
        c_inner = []
        desc1_hover = []
        desc2_hover = []
        desc3_hover = []
        desc4_hover = []
        c_hover = []
        min_ = 0.0
        max_ = 0.0
        for x, value in d.items():
            for y, value in value.items():
                max_ = max(max_, value)
                min_ = min(min_, value)
        if set_max != None:
            max_ = set_max

        colors = heatmap_palette(light_spec_approximation, 127)
        w = 0
        h = 0
        w_keys = sorted(d.keys())
        width = 0
        height = 0
        if len(w_keys) > 1:
            width = w_keys[1] - w_keys[0]
            h_keys = sorted(d[w_keys[0]].keys())
            if len(h_keys) == 1:
                height = 1
            else:
                height = h_keys[1] - h_keys[0]
            for x, value in d.items():
                for y, value in value.items():
                    x_hover.append(x + width/2)
                    y_hover.append(y + height/2)
                    desc1_hover.append(str(x))
                    desc3_hover.append(str(y))
                    if inner != None:
                        if inner[x][y] == None or max_ == 0:
                            c_inner.append("#AAAAAA")
                            desc2_hover.append("")
                        else:
                            c_inner.append( colors[min(int(126*(inner[x][y]-min_) /max_), 126)] )
                            desc2_hover.append(str(inner[x][y]))
                    else:
                        desc2_hover.append(str(d[x][y]))
                    if desc2 != None:
                        desc4_hover.append(desc2[x][y])
                    if d[x][y] == None or max_ == 0:
                        c_hover.append("#AAAAAA")
                    else:
                        c_hover.append( colors[min(int(126*(d[x][y]-min_) /max_), 126)] )

        color_mapper = LinearColorMapper(
                palette=heatmap_palette(light_spec_approximation, 127),
                low=min_,
                high=max_
            )

        tick_formater = FuncTickFormatter(code="""
            return Math.max(Math.floor( (tick+1)/2),0) + '; ' +
                    Math.max(Math.floor( (tick)/2),0)"""
            )
        #tick_formater = FuncTickFormatter(code="return 'a')

        source = ColumnDataSource(data=dict(
            x=x_hover,
            y=y_hover,
            c_outer=c_hover,
            c=c_hover,
            desc1=desc1_hover,
            desc2=desc2_hover,
            desc3=desc3_hover,
        ))
        hover = HoverTool(tooltips=[
            ("#indel", "@desc1"),
            ("#mut", "@desc3"),
            ("val", "@desc2"),
        ])
        if desc2 != None:
            source = ColumnDataSource(data=dict(
                x=x_hover,
                y=y_hover,
                c_outer=c_hover,
                c=c_inner,
                desc1=desc1_hover,
                desc2=desc2_hover,
                desc3=desc3_hover,
                desc4=desc4_hover,
            ))
            hover = HoverTool(tooltips=[
                ("#indel", "@desc1"),
                ("#mut", "@desc3"),
                ("val", "@desc2"),
                ("fails", "@desc4"),
            ])

        plot = figure(title=title, plot_width=200, plot_height=250, tools=[hover])
        plot.axis.visible = False
        plot.toolbar.logo = None
        plot.toolbar_location = None

        plot.xaxis.formatter = tick_formater
        plot.xaxis.ticker = FixedTicker(ticks=[0, 4, 8, 12, 16, 20])
        if show_coverage:
            plot.rect('x', 'y', width, height, color='c', line_alpha=0, source=source)
        else:
            plot.rect('x', 'y', width, height, color='c', fill_alpha =0, line_alpha=0, source=source)
        plot.circle('x', 'y', radius=width*2/7, color='c_outer', line_alpha=0, source=source)

        #color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))

        #plot.add_layout(color_bar, 'left')
        plot.min_border_bottom = 50

        return plot

    for approach_ in approaches:
        approach = approach_[0]
        accuracy, coverage, runtime, alignments, fails = getAccuracyAndRuntimeOfAligner(db_name, approach, num_tries, allow_sw_hits)


        tot_runtime = ""
        # fix for old databases
        ## runtime_tup = [[getTotalRuntime(db_name, approach)[0][0]/ (32*11*86)]] 
        runtime_tup = getTotalRuntime(db_name, approach)
        if len(runtime_tup) > 0:
            tot_runtime = " [" + str(runtime_tup[0][0] * 1000 )[:5] + "ms]"

        print("plotting", approach)
        json_save.append( (approach, accuracy, coverage, runtime, alignments, fails, runtime_tup) )

        plots[0].append(makePicFromDict(accuracy, approach + tot_runtime, desc2=fails, inner=coverage, set_max=1))
        plots[1].append(makePicFromDict(runtime, None, set_max=4)) #, set_max=50
        plots[2].append(makePicFromDict(alignments, None, set_max=500))

    sw_accuracy, sw_coverage = getAccuracyAndRuntimeOfSW(db_name)
    plots.append([makePicFromDict(sw_accuracy, "sw accuracy", inner=sw_coverage)])

    with open(out + ".json", "w") as f:
        json.dump(json_save, f)

    save(layout(plots))
    print("wrote plot")

def analyse_all_approaches(out, db_name, query_size = 100, indel_size = 10):
    output_file(out)
    approaches = getApproachesWithData(db_name)
    yax = True

    # list of tuples:
    # (
    #   w,
    #   h,
    #   approach name,
    #   accuracy picture,    <- matrix with accuracy as values
    #   runtime picture,     <- matrix runtimes as values
    #   mapping quality list <- tuple (accurate, inaccurate) of vectors with map. qual. as values
    # )
    out_list = []

    for approach_ in approaches:
        approach = approach_[0]
        results = getResults(db_name, approach, query_size, indel_size)
        hits = {}
        tries = {}
        run_times = {}
        mapping_qual = ([], [])

        max_indel = 0
        if indel_size > 0:
            max_indel = 2*query_size/indel_size
        max_mut = int(query_size * 4 / 10)

        def init(d, x, y):
            if x not in d:
                d[x] = {}
            if y not in d[x]:
                d[x][y] = 0
            return d

        for score, score2, optimal_score, result_start, result_end, original_start, num_mutation, num_indels, num_seeds, num_seeds_chosen_strip, seed_coverage_chosen_strip, seed_coverage_alignment, mapping_quality, nmw_area, run_time, sequence, sample_id in results:
            hits = init(hits, num_mutation, num_indels)
            tries = init(tries, num_mutation, num_indels)
            run_times = init(run_times, num_mutation, num_indels)

            #print(result_start, original_start)
            if near(result_start, original_start, result_end, original_start+query_size):
                hits[num_mutation][num_indels] += 1
                mapping_qual[0].append(mapping_quality)
            else:
                mapping_qual[1].append(mapping_quality)
            tries[num_mutation][num_indels] += 1
            #in the current test scenario the total runtime is saved in every sample of the cell
            run_times[num_mutation][num_indels] = run_time

        def format(d, w, h, divideBy):
            pic = []
            if len(d.keys()) == 0:
                return [ [] ]
            else:
                w_keys = sorted(d.keys())
                h_keys = sorted(d[0].keys())

                for x in w_keys:
                    pic.append( [] )
                    for y in h_keys:
                        if x not in d or y not in d[x] or (not divideBy is None and divideBy[x][y] == 0):
                            pic[-1].append( float("nan") )
                        else:
                            if divideBy is None:
                                pic[-1].append( d[x][y] )
                            else:
                                pic[-1].append( d[x][y] / divideBy[x][y] )
            return pic

        out_list.append(
            (
                max_mut, 
                max_indel, 
                approach, 
                format(hits, max_mut, max_indel, tries),
                format(run_times, max_mut, max_indel, None),
                mapping_qual
            )
        )
    print("[COPY HERE]")
    print(out_list)
    print("[END COPY HERE]")

def compare_approaches(out, approaches, db_name, genome, query_size = 100, indel_size = 10):
    output_file(out)
    plots = []

    max_indel = 2*query_size/indel_size
    max_mut = query_size

    def get_data(db_name, approach, query_size, indel_size):
        results = getResults(db_name, approach, query_size, indel_size, genome)
        hits = {}
        tries = {}
        scores = {}

        def init(d, x, y):
            if x not in d:
                d[x] = {}
            if y not in d[x]:
                d[x][y] = 0
            return d

        for score, result_start, result_end, original_start, num_mutation, num_indels in results:
            hits = init(hits, num_mutation, num_indels)
            tries = init(tries, num_mutation, num_indels)
            if not score is None:
                scores = init(scores, num_mutation, num_indels)
            if near(result_start, original_start, result_end, original_start+query_size):
                hits[num_mutation][num_indels] += 1
            tries[num_mutation][num_indels] += 1
            if not score is None:
                scores[num_mutation][num_indels] += score

        return (hits, tries, scores)

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

    scores = sub(convertToList(l_1[2], max_indel, max_mut, l_1[1]),
            convertToList(l_2[2], max_indel, max_mut, l_2[1]))


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
            #Color palette for visualization.
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

        tick_formater = FuncTickFormatter(code="return Math.max(Math.floor( (tick+1)/2),0) + '; ' + Math.max(Math.floor( (tick)/2),0)")
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
    avg_score = makePic(scores, max_mut, max_indel, "score " + desc)

    if not avg_hits is None:
        plots.append(avg_hits)
    if not avg_score is None:
        plots.append(avg_score)

    save(row(plots))


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

    def get_random_pos_queries(l, amount, pack):
        if 4**l <= amount: # in this case there are not enough samples...
            print("computing all samples")
            for q in get_all_queries(l):
                yield q
        else: # compute the random unique samples
            print("computing samples from random location")
            already_seen = set()
            already_seen.add("")
            max_pos = pack.unpacked_size() - l - 1
            tries = 0
            for index in range(amount):
                sequence = ""
                while sequence in already_seen:
                    tries += 1
                    pos = random.randint(0, max_pos)
                    while pack.is_bridging(pos, l):
                        pos = random.randint(0, max_pos)
                    sequence = str(pack.extract_from_to(pos, pos+l))
                assert(sequence != "")
                already_seen.add(sequence)
                if (index+1) % 10000 == 0:
                    print(100*(index+1)/amount, "% num_tries:", tries/10000)
                    tries = 0
                yield sequence

    fm_index = FMIndex()
    fm_index.load(reference)

    pack = Pack()
    pack.load(reference)

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
        print(l, " of [", min_len, ", ", max_len, ")")
        num_queries_actual = 0
        data1.append( [] )
        for _ in range(r1max):
            data1[-1].append(0.0)
        data2.append( [] )
        for _ in range(r2size):
            data2[-1].append(0.0)
        for q in get_random_pos_queries(l, num_queries, pack):
            ambiguity = fm_index.get_ambiguity(NucSeq(q))
            #in this case an ambiguity of 0 might occur
            if l <= 4**num_queries and ambiguity == 0:
                continue
            #if not fm_index.test_sa_interval(NucSeq(q), pack):
            #    print(q)
            #    print("found error")
            #    exit()
            if ambiguity == 0:
                print(ambiguity)
                print(q)
            if ambiguity < r1max:
                data1[-1][ambiguity] += 1.0
            elif int(math.log2(ambiguity - r1max+1)) < r2size:
                data2[-1][int(math.log2(ambiguity - r1max+1))] += 1.0
            num_queries_actual += 1
        for index, number in enumerate(data1[-1]):
            data1[-1][index] = number / num_queries_actual
        for index, number in enumerate(data2[-1]):
            data2[-1][index] = number / num_queries_actual


    print(indices1)
    print(indices2)
    print("[copy here]")
    print( (min_len, max_len, data1, data2) )#the tuple to be copied
    print("[copy here end]")

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

def compute_bam_bai_for_reads(name):
    genome = "/MAdata/genome/GRCh38.p12_14_07"
    fasta_file_name = name + ".fastq"

    sam_tools_pref = "~/workspace/samtools/samtools "
    temp_prefix = "/mnt/ssd0/sra_reads/temp/" + fasta_file_name
    res_prefix = "/mnt/ssd0/sra_reads/res/" + fasta_file_name

    download_cmd = "~/workspace/sra_toolkit/bin/fastq-dump -O /mnt/ssd0/sra_reads/temp/ " + name
    print("(1/5) executing:", download_cmd)
    os.system(download_cmd)

    ma_cmd = "~/workspace/aligner/ma -x " + genome + " -i " + temp_prefix + " -o " + temp_prefix + ".sam"
    print("(2/5) executing:", ma_cmd)
    os.system(ma_cmd)

    to_bam_cmd = sam_tools_pref + "view -Sb " + temp_prefix + ".sam > " + temp_prefix + ".bam"
    print("(3/5) executing:", to_bam_cmd)
    os.system(to_bam_cmd)

    sort_cmd = sam_tools_pref + "sort -m 96G " + temp_prefix + ".bam > " + res_prefix + "_sorted.bam"
    print("(4/5) executing:", sort_cmd)
    os.system(sort_cmd)

    index_cmd = sam_tools_pref + "index " + res_prefix + "_sorted.bam > " + res_prefix + "_sorted.bam.bai"
    print("(5/5) executing:", index_cmd)
    os.system(index_cmd)

    print("done")

def filter_duplicates(in_f, out_f):
    with open(in_f, "r") as in_file:
        with open(out_f, "w") as out_file:
            read_names = set()
            bWrite = True
            for line in in_file:
                if line[0] == ">":
                    bWrite = not line in read_names
                    read_names.add(line)
                    print(line, bWrite)
                if bWrite:
                    out_file.write(line)

def compute_bam_bai_for_files(path_name_gen, task_name, pack, reference):
    sam_tools_pref = "~/workspace/samtools/samtools "
    merged_prefix = "/mnt/ssd0/sra_reads/temp/" + task_name
    res_prefix = "/mnt/ssd0/sra_reads/res/" + task_name

    merge_list = ""

    for path, file_name in path_name_gen:
        print("working on file", file_name)
        temp_prefix = "/mnt/ssd0/sra_reads/temp/" + file_name

        ma_cmd = "~/workspace/aligner/ma -x " + pack + " -i " + path + "/" + file_name + " -o " + temp_prefix + ".sam"

        mm_cmd = "~/workspace/minimap2/minimap2 -a " + reference + " " + path + "/" + file_name + " > " + temp_prefix + ".sam"

        ngmlr_cmd = "~/workspace/ngmlr/bin/ngmlr-0.2.8/ngmlr -t 32 -r " + reference + " -q " + path + "/" + file_name + " -o " + temp_prefix + ".sam"
        #os.system(ma_cmd)
        os.system(mm_cmd)
        #os.system(ngmlr_cmd)

        print("converting sam > bam                ")
        to_bam_cmd = sam_tools_pref + "view -Sb " + temp_prefix + ".sam > " + temp_prefix + ".bam"
        os.system(to_bam_cmd)

        merge_list += temp_prefix + ".bam "

    merge = sam_tools_pref + "merge -f " + merged_prefix + ".bam " + merge_list
    print("merging")
    os.system(merge)

    sort_cmd = sam_tools_pref + "sort -m 96G " + merged_prefix + ".bam > " + res_prefix + "_sorted.bam"
    print("sorting")
    os.system(sort_cmd)

    index_cmd = sam_tools_pref + "index " + res_prefix + "_sorted.bam > " + res_prefix + "_sorted.bam.bai"
    print("indexing")
    os.system(index_cmd)

    print("done")

def all_files_with_extension(path, ext):
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith(ext):
                yield root, f



#createPacBioReadsSimLord("/MAdata/db/testPacBio_human.db")
#createPacBioReadsSimLord("/MAdata/db/testPacBio_zebrafish.db")
#exit()


# extract fasta from bam:
# samtools view <in.bam> <chr>:<from>-<to> | awk '{OFS="\t"; print ">"$1"\n"$10}' - > <out.fasta>

#make_split_read_db_close(human_genome, 1000, 0.05, "split_read.db")
#make_split_read_with_del(human_genome, 5000, 1400, 0.01, "del_read.db")
#test_my_approaches("del_read.db", human_genome, specific_id=1 )
#exit()

#### compute_bam_bai_for_files(
####         all_files_with_extension("/mnt/ssd0/GIAB_HG002/", ".subreads.fasta"),
####         "test"
####     )
#### exit()

## filter_duplicates(
##     "/mnt/ssd0/PacBio/ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/" + 
##             "HG004_NA24143_mother/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/section.fasta",
##     "/mnt/ssd0/PacBio/ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/" + 
##             "HG004_NA24143_mother/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/section_filtered.fasta"
##             )

## compute_bam_bai_for_files(
##         [(
##         "/mnt/ssd0/PacBio/ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/" + 
##             "HG004_NA24143_mother/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37",
##         "section_filtered.fasta")
##         ],
##         "hg002_hg37_mm",
##         human_hg37_genome_full,
##         "/MAdata/chrom/human_hg37/full_sequence.fasta"
##     )
## exit()

# create_split_reads(human_genome, 1000, 10, .05, "split_reads_experiment.fasta")
# compute_bam_bai_for_files([("~/workspace/aligner","split_reads_experiment.fasta")], # "split_reads_experiment")
# exit()

## minIon reads

### #compute_bam_bai_for_reads("ERR2407705") # done
### compute_bam_bai_for_reads("ERR2407662")
### compute_bam_bai_for_reads("ERR2407706")
### 
### ## pacbio reads
### 
### compute_bam_bai_for_reads("SRR7515657")
### compute_bam_bai_for_reads("SRR7515658")
### compute_bam_bai_for_reads("SRR7515660")
### compute_bam_bai_for_reads("SRR7515664")
### compute_bam_bai_for_reads("SRR7515668")
### compute_bam_bai_for_reads("SRR7515669")
### 
### ## illumina reads
### 
### compute_bam_bai_for_reads("ERR1157312")
### compute_bam_bai_for_reads("ERR1157308")
### compute_bam_bai_for_reads("ERR1157301")
### compute_bam_bai_for_reads("ERR1157297")
### 
### exit()

##
# RUN SW for one sample
#
def run_sw_for_sample(db_name, genome, sample_id, gpu_id=0):
    sequence, _ = getQueries("/MAdata/db/" + db_name, sample_id)[0]
    
    print("loading pack...")
    ref_pack = Pack()
    ref_pack.load(genome)    
    ref = ref_pack.extract_complete()
    print("done")

    print("gpu computation...")
    alignment = libMA.testGPUSW(ContainerVector([NucSeq(sequence)]), ref, gpu_id)[0]
    print("done")

    print("SW-score:", alignment.iMaxScore)
    for maxpos in alignment.vMaxPos:
        print("SW-pos:", maxpos)

    origin, _ = getOrigin("/MAdata/db/" + db_name, sample_id)
    print("Orig-Pos:", origin)

    print("Recorded Optima:")
    for optima_start, optima_end in getOptima("/MAdata/db/" + db_name, sample_id):
        print(optima_start, "-", optima_end)


def check_sw_score(db_name, genome, sample_id, ref_start, gpu_id=0):
    sequence, _ = getQueries("/MAdata/db/" + db_name, sample_id)[0]
    
    print("loading pack...")
    ref_pack = Pack()
    ref_pack.load(genome)
    ref = ref_pack.extract_from_to(ref_start - len(sequence), ref_start + len(sequence)*2)
    print("done")

    print("gpu computation...")
    print(len(ref), "x", len(sequence))
    alignment = libMA.testGPUSW(ContainerVector([NucSeq(sequence)]), ref, gpu_id)[0]
    print("done")

    print("SW-score:", alignment.iMaxScore)
    for maxpos in alignment.vMaxPos:
        print("SW-pos:", maxpos + ref_start - len(sequence))

### origin, _ = getOrigin("/MAdata/db/human_30000_10.db", 648)
### print("original score:")
### check_sw_score("human_30000_10.db", human_genome, 648, origin)
### ref_pos = getRefPos("/MAdata/db/human_30000_10.db", 648, "MA Fast")
### for pos in ref_pos:
###     print("aligner score for position:", pos)
###     check_sw_score("human_30000_10.db", human_genome, 648, pos[0])
### exit()

#get_ambiguity_distribution(plasmodium_genome, 1, 100)
#get_ambiguity_distribution("/MAdata/genome/human_bugged", 1, 100)
#get_ambiguity_distribution(plasmodium_genome)
#print("random")
#get_ambiguity_distribution(random_genome)
#print("plasmodium")
#get_ambiguity_distribution(plasmodium_genome, 7, 17)
#print("mouse")
#get_ambiguity_distribution(mouse_genome)
#print("human")
#get_ambiguity_distribution(human_genome)
#exit()


#FileReader.testBufReader()
#exit()

#libMA.testKsw()
#exit()

#libMA.test_ransac()
#exit()


# ================================================================================================ #
# making the sw verified samples                                                                   #
# ================================================================================================ #

# plasmodium

## task 1: (done)
#createSampleQueries(plasmodium_genome, "sw_plasmodium_200.db", 200, 20, 2048, reset=True, gpu_id=0)

## task 2: (done)
#createSampleQueries(plasmodium_genome, "sw_plasmodium_1000.db", 1000, 100, 2048, reset=True, gpu_id=1)

## task 3: (aborted)
#createSampleQueries(plasmodium_genome, "sw_plasmodium_30000.db", 30000, 100, 2048, reset=True, gpu_id=0)

# human

## task 1:
#createSampleQueries(human_genome, "sw_human_200.db", 200, 20, 100, reset=True, gpu_id=0)

## task 2:
#createSampleQueries(human_genome, "sw_human_1000.db", 1000, 100, 50, reset=True, gpu_id=0)
#createSampleQueries(human_genome, "sw_human_3000.db", 30000, 100, 32, reset=True)



# zebrafish

## task 1:
#createSampleQueries(zebrafish_genome, "sw_zebrafish_200.db", 200, 20, 32, gpu_id=0)

## task 2:
#createSampleQueries(zebrafish_genome, "sw_zebrafish_1000.db", 1000, 100, 32, gpu_id=1)

#exit()

# ================================================================================================ #
# end making the sw verified samples                                                               #
# ================================================================================================ #

#test_my_approaches("zebrafish_n_200.db", zebrafish_genome_n, be_mean=True)
#analyse_all_approaches_depre("zebrafish_n_200.html","zebrafish_n_200.db", num_tries=1)

# ================================================================================================ #
#   DATABASE NAME                                                                                  #
# ================================================================================================ #


#db_name = "plasmodium_30000.db"
#db_name = "sw_plasmodium_200.db"

#createSampleQueries(zebrafish_genome, db_name, 200, 20, 32, gpu_id=0, smaller_box=True)

#resetResults(db_name)
#run_sw_for_sample(db_name, zebrafish_genome, 1, gpu_id=0)

#test_my_approaches("None", "/MAdata/genome/eColi_full", specific_query="AAAATGCCTTATGCTTAATGCGGGATATCCATCATTTCATCATGATGCTTATGGGAGCGGCTGGAACCAATATCACCATGCGCTGTCTCTTCGGATGAGTCCCCCTGGTCTGGTGCCTGGGCCTTGATGTTGCAGCAGACTTCAGGGTCACTCTAGACGACTATATCCTTTATCCATGGCATGGAATCGAATAGGTTCGTATTCAACCCGGCGAGTATATATTCGGCTGAACTGCGTGCAACGCTGAACGTTAAACATGCGAATCTGGGGGAAGGGATGGAAACTCGAACCCCCTGGCTGAGGGCGGCTGTGGCGTCCAAGAAATTTGTCGATGTAACCGGGTGAAGGTGAATATGACGGGTAATTATCAGTCAAAATGATTTAGTCCGGGGCAGACGTGGAATATACCAGGCAGGGCTATTCTAAAGCCTCAATTCAAGCTAGGTACCGTTTAGCGGACAATACTCGTGTGGCGGTCTAGCCATGGTTCGGTGTGGAATTGCCCCGTGGTACACGGGCGTGGTGGTGTGAACTGTGTCGTCTCGACATCAACTGACTTAAACTGGCTTCGGCGGTTTTCGTTTACAGGAATGTGATAGGGAAAATCTGGAGCGCATACAAAAAGGCCGCCACGTTAACGGTCGGGAACCTCTCCTATGGGCAGAACAAATGAATTAATGAAATTGCGGGGTTATCATCTCCCAGGTATTAATCCATACTAATCACAATAAGGTATATTTACTTCAACCAGGCATAAACCATTTTGGATTTTGTGCGTGGGAACAGCCTTAACTGTGTTAAAGGGGGAGGTGGAAAATAGAATGAGGAGTATCAGCAAGTAATACTCGTCCGCTTATACCACAACGTGGGATGAGAGGGATATGAAAAACATCAAGGCATGAGATAAGCTCTGCCTTGAAGAATAAATTCGCTGTTTACAGCGGGCTACTTGATTTCACTGTCTATCTGCAATTTCCGGTAGACGGTCAACGTTCAGCACTTCCATAGACATAGAAATATCCCCAATCGATTTCGCTCAAGCTGTACTGTTAGACCACTCAGGATCAATTTGCACGTTTATTTACAAATTTGACCTCAAGAAAATATCTTTACGCAACTGCGGCAGAAAGTGCGGTTCTGCAATCTGACACTGGGAGGCGTTCAGCCAACATCTAATCTGTGTGCCATGCTTTTTGCAATGTTTGGCTGTGTTTTTTTTCGCGAGGAGATAAAGAAAATCGAGTAATGTTAACATACACTATATCCTCCGAACAAGGCGTTTGCAGCAACGCCACTTTGCTTTCTTCTTCTAAGAAGCGGCAAAGGACGCTTCTTCTCCCAACAGACAGTGCTACGGTATCTCGCGTAGGCTTTACCACGCATTCTGGACGTTTATCGGTCGAAGAATGGACGCTGGTTCACCCTTGGAATTAGAGGCGCGCAAACTGAATTGCTCCTTCTGGGTTCACGCCGAACGAGTTTGATGCGCAGGCTCTGCCAGCAACATCTTCCGTGCTCAGCAAGTCACCTCTGCTTACCGTACAACCTGGTGTTATAGCCGTGTAACCTAGCACGGTTGGCTCTGTTAATAAGGGTTCTTCGCACATTTTCGTCGCGGTCGATGATTTCGACAGCAGATATGCCTAAGAATACGGCAGGTCGCGTACTGAGGAGTACCTTCCGGTTGGTGGCCAATTATGGCTGTCGTCGTGCAACGAATGATGAATTGCCATTAACGTACCTGGTTTCAATCCTCTTGCAGGGGAGTCACAAACGAGTACAATTCAAGAGATCCATCGCTTTCAGACATCAAGAACTTTGGCGACCCCTTCAGGGTGAGGGCACTTTATCGCGTGTTTGCGATGCGGCAGAATACGTAGAGATATTCAGTACGCTTTATCTAATTAACTGCCTGATTTAGCGTTGCATCGCCCTGAATGACGTTGACGAACATCGTAAACCGACCCGGCGTTCATCACTCCGATATCAGGTCGAGATTACGTCAGGCCGACTATCAAAACCTATCACGCACATTTTCTTTCCCACCTCATGGGCTCCAAAACCAGTGGCGATGCGCCGCGCAGGAGGTAACGTTCGTTACCAACACCCCCTTTCGCCGAAGTGATACAACATATGCGTGCCCATAAAATTCCTTGTTAATAAAGTGGATCAATTTAAACAGGTTGAACGGGTCAAGAGCTGTTTTCGACTAGACTGCCGTTCTGCGCCCGGCTTTGCCCATAAAAAGTTCTCCGCTGGGATTTGATCACTCAAGCACGTGATTCACTCCGATGGACAGCCAAGTTCACGCACATCAGGTTCGAGTAACAAAAGTATTTGAGTTTACCCGGTCCATCCACTTGCCCCTGCCAGCGCACGACCGCATGCTGCCATAGACCAATGAATGTTCCGATCGGAAATCAATCGGCCCCAGCCGCTACACTGTGGCTTGCAACAATCACGATCACATTGTGAGGCATCTAAATAGCTGACCGAACGCCGGGGATCATTAAACGCGTTTTTTGTGACGGCGTGTATTTGGCGCTCGGAACCTGCCCGGTGTGGGCAGCTCGGCTGGCGTCCTTTTATCTTTCCTTCGCAGGATAGCAGCCCCATCTTTCAATTTTTCGGGCTTTAAGTTGCGCCATCTTTGCAGCCACTTACTGCCAAGTAACCCGCATACACCCGGTTGCCGAAACCGCCTTATGCATCGCTGACCAAGTTTTTCCGGGCTTCCCAGTGCACTGACGTTGAAGTACACGGGGGGGCGGATGTGTTTAAAAATGCGGGGCCTGACGATGTTTGTCTTCCATGCGCCTCGACTGTGATGAACTTAAGGTCTGCCTCATGCAGATGAAACGCACAGATTATAGGAAGCTTACTAGCCTTGAAGCTCGATTGGCGTGTTTAGTACATCCTGGCCTTACTCATTTGCTTTTAATCATACGCTCAGGCGCGCGATTGATGGTTCCGAAGACTATAACGGGCTGTATTACTAGTCTGGGATCTATATTGGAGGCCAGTCGACCTCTCGCCATTTATTAGAGTAAAAGTCTATTCTGTGATAAAATGGCGCTTGGATTCATAGCTTCAAAAATACTCCCTGCAATTCAACCCATTTCGCACTCGTTCGAATACTTTGATTGTTGCTTACTTTACTCTTCCTATTTTCCGCACCCTAATCTGACAGCAGCTATCAGACCTTTTTAAATATCCAGTATGAGGAAAGTGTACTTGTTTTCACCAGACGCCAGGCTTTCTTAACCATTCCATTTCTTTAGCGGCATCATTAAAAACGCGATACCCGCTTATCTAAATCCACCAAGAAATTGTCTATCAGTAGAAACCGCCAGAAGACGTCTCGGAAACTTCAGTCAGGTCTCGCGAAACCTCGTCCAGGTCAGCTGGTTTTGCATTACGGATTGCATAGATAAATGAGCATAAGAGTTGTGACGGTTGGATCTATCAGGATGTCGCTATGAGATATGATCTGTATTTTTTATCGATTGTCAATGGTGGGGGGGATTGCTCGGTTTCTCATCAATGGTCGTCTGAAATGGACGCTGGCCGGCACTTTCCGAATCCTCCGGCCGTTGCGGCTTGTTGCGACCCGTAATTGTGTGGAACATCCAGATAGGTCCATCTACCCAGAAGTTGGACAATGGAAACCATCCGGCGACTAAAGGCCGCTGAAGAATTACCAATCATTACAAGGCTTACGCTCGGCCATTCGAAGTCTTAACATTCGCTGACGGTGTATTTTGCGATGTTAAAGAAACTATATCGCTGATGTCGTTATTTTGGTGGTGGATGGTAAAGCGGAGTGGAAAACATCAATCAAACCCCCAAAAGGATGTTATTTACGTGTCCGTGAAGTGATTTCTTGTTTACATCATTGGTCTGATCGATCATCTAACCCTTGAGCGACAGAAGAGAGTGGGGCAATTAACGTATTATGACGGGAGCAGTCCCTGCCGTCGGTTGAACATGTTACTTAATTTTGAAATTATCCGGACGCTATATATTCTTTCAGTTTTCTGTTTTTACGGGCGTGCGTTTCGGGTTTTGCACACACTCGCAGCCTCATATTCTCAGAAGCAGGGAGTCAATATTGTCGTGCCTATCGGGTTCTGGGGAAATAGATCTGCGATAAACACTTTGAATCATCCGATTTCTGCTGGCCCGAAGAGGGTAAGTTATCAAGGCGCCTTCATGAGATATCTCAGCTATGCACGTTACTCTGGACCCTTTTTTGGCCTAGAAGCGTATATGTATCGTAATCGTATCCCTTTATTTCGCCGGAGTACTTGAGCGCGCTGCTATGGCCTTTTTGAACTCATAGCATTGACATTCCTTCACGCAGCAAAGCTGTGGGATGGTTTAGCTAGTGGCAAAGGCAGCAATACAGCCTTTCCGACGATGCTCTTGATTGTTCATGGTAGGTCTCTCTGTTCAGAATAGTACGATGAACTGTTATAATATAACAATCCCTAACGGGAGAAGTCTGTGCAGTGGCTTTCCTGCGTTCGGTAAAATCACGAGTCTGCGGCCGATTTTATGCGGCGGGGGTTTACTTATAAGGGGGATGGATCAACGAACGCCGTCAGAAAATCTACACGCCTAGTTGGCTTCGGTGGGCGCAGGTTGTTTCATGATGTATAGCTTATGCCAGTTAAACCTGGCGATCGTTAAAGAGCAAAAGGTATGTTTCTTGTGGTGATTTATCGAAGCCAGTCAAGGTGACAGAGACCTATTAAATTGTGAAAAAAAGACGATCTTTTCGCGTCGTCTGCCAGAGAACTGATGGAGGAGGTGTTTGGTGCCCTCAGTTAGTCGATGATTCTGCCCTGGACTGGGCGGAAAGAAGACTGGTCAAATGCCGATATTGAACAAGGGTTAAATAGGAAACATTAATCGAGGCAAGGTTACTATTTTGCATTCCCGCACACCACCTGAAGATTTGCTGAGAGCAACATCTTTCCGTAATGGGGCAGAAAACAGCGACACTAACAACATAACTGGCATTATCGCCGGTGCGGTGGCATTATCTTTGTCGGCGCGGTTTCATATCCACGATGCAACGTGACGGGGAACGTTACTGTATCACACACACAACTTGGCAAGGCTGCGGCTGCTGGATCCTAGCTCCGGTGGTAAGTAGTGTTTGTGGATTGGCTCAGTAACTTGGCCAGAAGCGATATTAAGAAGATGGGCAGCGCAAGTGCCCGAAGAGCCGAGTGCTGTTTACTTAAACCAGAAACGGGCACTATGCAATCTGCGGCAGACCATTTGGCGATCCCATCCGATGCTCGGATTCAGTACATGCTATGAAGTCGAACTGGCGTGTGGATTGGCGACGACAGTGCGTCAGGCTACGAGAAGAAGCATGTCCGCAAAGCCATTGCCGTTAGGCAGTGGCGCTCGATCTGAATTGCGTGTTGTAGGGATAAGAAATGTAGAGAAAGCCCGGGCGAGCCCTGGGAAAAGGTAAAGCCGTTCTGATAACTCTATGATGCCCGCTTTCCGGGCTACATTCCCGCGGTCGGACCATTACCGGCTGTCACCGTCAAATACAAGCGCTGAGCATGAAGACGTAAACGGGGAACAACGCCAGCAAGGTACGACTGCGGACATGATGCCTAAAATCGTTCCGTAATCGCTCATATGGCAAGTTTTTTACCCTCAAAGGCGGATGACGTTGTGTGACAGCCCACGACCATGATGGCTCGTCGCGGACCCGTTGCAAATGACGGTGAATGGCAGCTGAGCAAGCAACATTTCGTAGGGCATTCGTTTTGCACAACTCGCGTTTTGTAATGACAGTCTGTTTTGCCGCTTGAAAGCGGCAGGCAAAACTTGCATCGCCTGTGCCGACTGGTTATAAAGGTGCGTTTTAACAGTCAATTTGGCGGAACACCTATGAAGCGATGTATACTTTCTGGGCAATAGTAAAAGTCACTTGGAGAATGAGCGATGCGGAATGGAGTCTTGTGTGATGGTTGCGGTCAGGATTTCGACCCGTCGCCATAAACTGATGCATAGAAGACAGCTACGAAATATATCGTACTTCGACTCTAACGTGCCTGTCGCCAGGCTCGATATAAAACCTGTCAATGTCGGAACACGAACGTCGTTTTGAGGTTGGAAACCCGACTGCATTAATTTAAACCCGCTGAAAATCTGCCAACATTCAGAAATGGTCTGCCAATGCACCTTGCGTTATCGTTTGCATGTGCGAAGGGTAAAGTATTTAACGCTTGCGTGGACATCCGCAACTTACTGTGGTTCGAAGCACGGGACAGATGCATTTGTGAACTATCTGTCTTGCGTCATTCGCCAGTGACAAGGAATCAGAAGGTGATTGTACGTGGACAGGAATCATATCTTAAATAATAAACCTGAACTGCGGCGCAGTGAAAAATTTGAATAGCTGCTGATGGTGCGCTTTCTTAAGTATACTGATCTGATTGCTCGCCATTTCGTAGCCCTGTCTTTTGAAGCGCTCTGGTCATAATGAGAGTTACCCAACGCTTATTATCTACATGAGTCAACTTCTCAATCTTTAAAAAACAAGTTACTGTGTTTTTAAATGTGATAGAACCTTATGACCAGGAATACCAGGATGAGAATAATCGCTTTATCAGACTTCAGGTATCACTGCATAAGAGTGTCTTGTTTACCGTAGTCTTTTCTCGATATCATTTACAGATGCTAAATCGATTTTTTTGGGCCGCTTCATATAGCTGTAAGAGAAAGACATTAAATCATCATAATCAACCGTTAGAATCTGGTAGGTTCTCAGTTTCCGTTTTTATCTCACGATGCGGTAATTACGGTGCGTTAATTATACAATTGGCCAATCGAATATTCTTTTATTCGCTTGTTTAACCGTGTTAGAACGGGGTTGGTAAAAACAATCTGCAAGATTTTAACTTTGTTTTCAATCTTGGAATCATGCTTTCCTTCAACTCGCCCACACAAGTTATAACTAGGACTCGATAATCCCGAAGTGGACCGGCACGCACACCTGGCTGCGGCACCGTGCGACTTTAGCTTATTTCCTGATTTTATCTACTTCTGTGACTGTGGAAATAGCTTGCTTTTTTCTGAGAAAATCATTGTTATACCTGGCTTTATCTAAACCCAGCCGTTTTCCCGGAGCGTTGTTGAACTCTCTTGTGAGCTTCACCAGTCTAGGGATTTTTGCGCTTCATTCAGCTTCGACTGGATGCCGTCAGTCGCAGTACTCGTTAATGAGAATGTCTTTCTGGGCGGATGCTTTCTCTCAGTTGTATCATGATTAGCAAAATCATACGCTGCAGACAAATTTGCGTCGCAACACCACACCCATTCATAACACCTCTTGTTTTGCTTCAAAATAACTTATCCTGCTATCCATAGAGTAAGGTTTTAAACTTATGCGCCGACTAAAACGGAGGCTGCCTGTAGGAATACTCGCTGTATTTAACGCCACTTGACTCTTTTATGGTTTCTCAAAGGTCTAGCCAGGGGATGACCTGATCGACGAGTATTTTTGTTGATAAAGACTAAATGCTACATCTGACGGTTTCATTGTCGTTTTAACCTACTTACGTTTTTCTGCAAACGATTTCATGTCCATAGATCATTCGCCTCTTGTTAAATATAGTTAATAGTATTGTAATGAAACCTGCCGGATGTTTTACAACTATTAAATAAATTTACCTACACTAATTCATAGTTAGCCAGAGGCGTGGAAATCGTCAATGTCCTTTATTTCTATTTAAACCATGATAGAATATCTAACAATGTTTAATGTTCCTTATGGCGATTTGCTTCTTATTCTATTTTTTAGCCGGGTGATATTCTGTTTCCATTCTCCTGCTGAGAAGCGCTCGTTCGCCAGCAGAGGCACGTTAAGCACACAGATAAGAGACCGAAAATCAAGCAAGCAATGCGCGAGCACATATTGATTAGCCTGAGATCAGTATTGATCTGCTGGCAAGAGACAGACTATGGTATATAAAAACAGCATAAACTTCAGTCAGATATATTATGTTGTGTTCATCACAGCCTGCGGATCTTCCGACGCGAAAATTGATGACTATTTCCGCCATTTAGCAGATCTGTTCACAGTGTGGCTTCCTCTCACTCGACAGCAGAATATACGTTGAACAGCTGCAATCGATTCTGAATCAACTGTTGTCCAGCACCCAGCTGCGGACTTACTCTTCAAGCAAGTCTGGTGAGTTTCTATGATGTTGAGTGGTGGAATGTTAGTGACGTGATCCATACTGACTTGAACGAGTGGTATGAGCATTTAGCCGCTGGAGATCGACAGCCAGGCTGGTGCACAAGCTCTGGATTACGGCCATTTTATTGTATGTCTTCCATCAGGATTACATCAGTGAAAAAATGGCGTGGATGTGCATGCGTTAAAATGAACAGATGCTGAGAACTGCTACAGAGCGCGGCGCGCAGTACCTGCCGAGCATAACGTCCGGCCATCTTGTATCCAAAGCACCGGACGCACGTTGCAGAAGTTCTATCGCGAGAACGACTCCGACCAAACAGCATGAAATCACGGGGATCGGTAAAACCAGTGTAAACGCGAAAAACTGGGCAGGAAGTGGAGTAAAAATTACGGATGGCAGAGTATCGCCATCGCGAATTCACTTCAATTCTGTTCTGTGCCCGTCTGCCCCCGCCGCGCCATTTGGGCTGCTTTTTGTTTTTTATAGTACACCGCTTGCTGCCGGCACAAGCATCACTTTACCGGTATCAATCCAGGTACGCAGGCCAGGCTGCGCATCGGCAAATGGCGTATATTTGCCAACGCGTCCATCACTAACCAGCGCCACCGGTTTAGTTATTGATAACCGTCGCATACACCAGACAATGGCCCGCTCGCATTGGTAAAGCCGGTTTTGGTTAACTGAATATTCCAAGTTATGCGATACACCAGAGTGATTAGTATTGCGGAACGGCCAGCGTATACGTCGGATTAGAGAAGTTGCCATATCTTCCCGGGTCGGTACTTAACTGCCCGATCAACGGATATTGTTTGGCTGGCAATGAAGCAGTTTTGGTGTAAGTCCCGGGCAGTTGAAACGTTTAGCACCGACCAATCCGGTAGGTTCAACAAAGCGCGTGTTGTTCCATTCCGAGCGATCGCTTTCGCATTCATTGCCTTAATAAAGGCTTTGTAACCACCGGGGATAATGGTGCGCACAGCATTGCCGCCGGCGCGGTTTTCTGAAGACATCAGCGCCAGCAACAGCATATCTTTCACGGGCTGATTTTCCGCCTATTCACGTCGGGCTACGGCGCGAGATAGACCCCCTTTCATTCTTCCGGGTCTGGGCTATGATGATCCACTTTTCAAGGTTTTTCATCCAGGCAGTACGTGCATCCAGCACAACCCATCGCGGTCATTAATTTGCTGATAGACGCATTATCGGAGCGGCACCAGATCCGGGTGTTCAATAGATCACTTTTTGGTATTCAGACAACAATTCATCAGCGCTACCGGAGGCAATTTCCGGTTGTGAAGCGGTGTAGCGGCTGCACGTTTTGTCAACGGCCCTGCGTGCAAACAGGCACAGCCAGCCATCACGGGCCAGGCTAGTAATAAAGAAACTCGAAATGTTCGGCATGATGAGCATCAGATAGTGGTTCACCGGCGCACGGGTTGCGCACCGCCGGAGTAAGGATATTACTGAGGTTAGCGACGCCGATCATAACGAGACCAAAAAGTGCGATCGTCATAAGGAGAAATCGTGAGGAAATGCTGCATTGCTGACATTTACGCCAGCAATGCAACGTCAAACGAACTTTCTAGAACAAACGATGAACCGTTAGCCCCATGTATAACGGCTTAGGCGAGCAGCACTTCCAGTACCAGCACGCCAATCGCCAGCGCTCGAACTAGAGAAGCTAAGGCCTTCCTCTTTGTGATATTCAGGAAGCTTCGGAATACCAATGTAAAGCAGGGTAGCCGGTGTAAAACAGCGCCACCGTTGCGACCAGCGCACACACCAGACCATGGATATAAGCGCCACCAGACCGCATTAAAAACAGCGGAGTTGTTCAACGTAGCCCAGCGAAAAGACGCAGGCAGTGCGCAAGTGACGGACCGCTGCGGATAATTACGCCGCCATCCACCAGATGACCCGCCCCAATCATCCGTCGACCCCAGCCCAGCGATAACGCCGTAAAACAAGATCAGCCAGCGGCCAAGTCCGGTAAACCAGGATAACCTTCAGGATAGTGCCATCGCAACATTCCACGCCAATCTGTGTAGTGCCAATGAAGGCGCAAAATCAACCGGAATCGCGCCATCAGCAAAACGTGTTGAGGATGTAATGGTGAGAAATACGTTTCGTTTTCGCGATTAATCACCTGCATTTCACGATCGGGATGGGAAAAACCAGTCCCCAGACATGGCTCATACCCCGCCCTCCTTGTGTGTGAGCTTTCCAATGAACCTGACAGTTCAAGTTATAATGCAAGGCTTGTGCATTATTTTTGTGTTCGCCCATGAATTTTCACTGTCTGATGAACGTCCTTTTATCAGCGTGAATGATTCACAGGGTCGTAGTGCTTACTGGCAACCAAAGGGAGACAGACTGGCCTATGGATCTCAATTACACTTATCTCACAATATGGTTAATCCGCTGCTGGTGATCGGTAGCCTGGCGGAAGGTGAAACCGTGACTTACTGGGAGGCGTTGACGGCGCATCAGCGGGCTATTAAAGTTCCCACGCTGGTGGTACTTCTGTGGCGCTTGGCGGCATGATTGGCGACCAGTGTGCTTCTATCTGGCGGGCGGCGGTGTGGCGGCAAGCTGTTACGCCGTTTCTCGGAAACATCAGGATAAAATTTGAGCGGGCGCAGAAAACCTTATCAACGCCATCCGTATTGTTTGTCATTGGTGACGCGCTTTATGTATGCTTTCGGGTGATTGGCCGCGAACGCATGAGTTGGGTGCCAGCCGAGCTGCCGCCGAAATCTTTCTGCCGCTTGGAATATTCTCGGCGCATTTGCCTGGGACGTTGATTTTTAACCACTATTGCGTTTACGCTGGTGGTCAGGTGATTGCGCCGTAGGTTGCACAAACTCTACGACCATGCATTTGAAGCACTGGGTCTGGTTGATTTCTGGTGGCGTGGTTCTGGTGGGTGGGCGTGCGCTGGTGGCCTGAAACGACAGCTGGCGAAGAAAAGCCGGATCATCAGGCGTAAAACCATTGCCCCTGGATAAGGCTTTCCGCCGCACTCCGTACATTCCGGGGTACATCGCCTGATACGACGCTTGAGCCGTCTTTTCATGCCGACCAAAAATATTACCCAATTGAAAATATACGCACCTCACCGCGCTATTCTGGTTGAATGTGGATTCGCCAGACATAAAACGCCGCCATCCACTATCAACGACTGCCCGGTGGTGTAATTTGCGCCCTCCGAACAAAGCCACACCACCAGGCTGGCAATCTCATGCGTTGCGGCCAAAAACGCCGCATAGGGAATCGAAGGCTCCGCGTCGGGCTTCACGTCGCTGCATCCATGCCAGATTCATGGCGTTGGGCGATCGCCCCAGGCGCGGATCTGCGTTCTGACCAAAATCTTATGCCTGACCAGCTCCCAGCCCATCGCTTTTGGTTAACCCACCGAACGCATGTTTAGCGGCTGGTGTAGGCGCTGGCATCCGGCAGCGCTGTTATGTTCATAGTACACAGACGTTAATGTTGAGATGCGACCGCCCTGCCCTTGTTTCAGCACGATCTGACGAGCCGCAATTTGCGAGCATAAGAGATGCACCATTCGACATCAACGGTAAAAATCTTTGCGCCACTCATCCTAAAGAGCCATATCAAGAAAACGGCGCTTTGGTCATTGCACCCGCATTATTCACCAGCACATAATGGCCCCAGCCTGTTGAGATGAGGTTTCTCCAGCGCCAGGGCCCCTTAGGTAGATTGCCGAGATCCAGCATGGCACCATCTCCGCACGTACGGACCGTGGGCTAAGCTACCTCACGCGCGGTATCTTTTGCCCCTTCTTCATCTGAGTTGCCAGGTAATACCAACTACAGAACCCCAATGCTGCGCCAGTAATAACGCGCACTCTTTGCCGATCCCCGAATCGGAGGGCGGTAATAATCGCAACCTGTGCCATCGAGTTCTCCACTTAACGCTGAAATAAACGTTAAAGTATAGAAGGCGCATATCATCAGCGTTGTACCCCCCGCCGCAACGCACCAGTGAGTTGAATGGGAGGCATCCAGCCACTGCCCTTGCAATAACAGGCATTGGCTCCGGCTCACGCAGCGCGGGCGATTCTGGCTTCGCTGACGCGGGATACCAGAATGATGTCCCGCGTTACCAAAGCGGCGCCTGCGCAGAGACTTACCACCACGCAAGGCATCGCGCTCAATCTTGCGTCCTGATGCTGGTTTTTTCTTCCGCGCAGTGTC")
#exit()
#test(db_name, plasmodium_genome, only_overall_time=False, long_read_aligners=False, processor=31)

#analyse_all_approaches_depre(db_name + ".html", db_name, num_tries=1)


# ================================================================================================ #
# running through all sample sets                                                                  #
# ================================================================================================ #




# [7, 8, 6]:
# [0, 1, 2]:
# [5, 4, 3]:
# [9, 10]:
"""
for task_id in [0]:

    processor=2 #task_id*2

    data_set = [
        ("sw_plasmodium_200.db",  plasmodium_genome, False, True, 100), #
        ("sw_plasmodium_1000.db", plasmodium_genome, False, True, 10), #
        ("plasmodium_30000.db",   plasmodium_genome, True, False, 1), #

        ("sw_human_200.db",  human_genome, False, True, 0), # # 3
        ("sw_human_1000.db", human_genome, False, True, 10), #
        ("human_30000.db",   human_genome, True, False, 0), #

        ("sw_zebrafish_200.db",  zebrafish_genome, False, True, 100), # # 6
        ("sw_zebrafish_1000.db", zebrafish_genome, False, True, 10), #
        ("zebrafish_30000.db",   zebrafish_genome, True, False, 1), #
        
        ("sw_human_1000_10.db", human_genome, False, True, 10), # # 9
        ("human_30000_10.db",   human_genome, True, False, 1), #

        ("zebrafish_30000_10.db",   zebrafish_genome, False, False, 0), # 11
    ]

    db_name, working_genome, long_read_aligners, short_read_aligners, runtime_sample_multiplier = data_set[task_id]

    #createSampleQueries(working_genome, db_name, 30000, 10, 32)
    #resetResults(db_name)

    #test(db_name, working_genome, only_overall_time=True, long_read_aligners=long_read_aligners, short_read_aligners=short_read_aligners, processor=task_id*2, runtime_sample_multiplier=10)
    test(db_name, working_genome, only_overall_time=True, long_read_aligners=True, short_read_aligners=False, processor=processor, runtime_sample_multiplier=runtime_sample_multiplier)

    analyse_all_approaches_depre(db_name + ".html", db_name, num_tries=1)


exit()
"""

#split_reads("/MAdata/genome/eColi_full", "/MAdata/chrom/eColi/GCA_000005845.2_ASM584v2_genomic.fna")

split_reads(
        "/MAdata/genome/GRCh38.p12", 
        "/MAdata/chrom/human/GCA_000001405.27_GRCh38.p12_genomic.fna"
    )
#genome_dup_reads()
exit()

#test("sw_human_1000.db", human_genome, long_read_aligners=False, short_read_aligners=False, specific_sample=420)
#exit()


# ================================================================================================ #
# running blasr and graphmap                                                                       #
# ================================================================================================ #

# [7, 8, 6]:
# [0, 1, 2]:
# [5, 4, 3]:
# [9, 10]:

#first_accurate_SoC("sw_human_200.db",  human_genome, "soc_test")
#exit()

#createSampleQueries(plasmodium_genome, "sw_plasmodium_200.db",     200, 20, 32, gpu_id=0)
#createSampleQueries(plasmodium_genome, "sw_plasmodium_1000.db",   1000, 100, 32, gpu_id=1)

#createSampleQueries(e_coli_genome, "sw_eColi_200.db",     200, 20, 32, gpu_id=0)
#createSampleQueries(e_coli_genome, "sw_eColi_1000.db",   1000, 100, 32, gpu_id=0)
#createSampleQueries(e_coli_genome, "sw_eColi_30000_10.db",  30000, 10, 32, gpu_id=1)

for task_id in [18]:

    processor= task_id*2

    data_set = [
        ("sw_plasmodium_200.db",  plasmodium_genome, True, True, 100), # 0
        ("sw_plasmodium_1000.db", plasmodium_genome, True, True, 10), #
        ("plasmodium_30000.db",   plasmodium_genome, True, False, 0), #

        ("sw_human_200.db",  human_genome, False, True, 0), # # 3
        ("sw_human_1000.db", human_genome, False, True, 0), #
        ("human_30000.db",   human_genome, True, False, 1), #

        ("sw_zebrafish_200.db",  zebrafish_genome, False, True, 0), # # 6
        ("sw_zebrafish_1000.db", zebrafish_genome, False, True, 0), #
        ("zebrafish_30000.db",   zebrafish_genome, True, False, 0), #
        
        ("sw_human_1000_10.db", human_genome, False, True, 10), # # 9
        ("human_30000_10.db",   human_genome, True, False, 1), #

        ("zebrafish_30000_10.db",   zebrafish_genome, False, False, 0), # 11
        
        ("sw_eColi_200.db",   e_coli_genome, True, True, 100), # 12
        ("sw_eColi_1000.db",  e_coli_genome, True, True, 10), #
        ("eColi_30000.db", e_coli_genome, True, True, 1), #
        ("sw_eColi_30000_10.db", e_coli_genome, True, True, 0), #

        ("testPacBio.db", e_coli_genome, True, True, 0), # 16
        ("testPacBio_human.db", "/MAdata/genome/GRCh38.p12", True, True, 0), #
        ("testPacBio_zebrafish.db", "/MAdata/genome/zebrafish_full", True, True, 0), #
    ]

    db_name, working_genome, long_read_aligners, short_read_aligners, runtime_sample_multiplier = data_set[task_id]

    #createSampleQueries(working_genome, db_name, 30000, 10, 32)
    #resetResults(db_name)

    #test(db_name, working_genome, only_overall_time=True, long_read_aligners=long_read_aligners, short_read_aligners=short_read_aligners, processor=task_id*2, runtime_sample_multiplier=10)
    test(db_name, working_genome, only_overall_time=True, long_read_aligners=False, short_read_aligners=False, processor=1, runtime_sample_multiplier=runtime_sample_multiplier)

    analyse_all_approaches_depre(db_name + ".html", db_name, num_tries=1)