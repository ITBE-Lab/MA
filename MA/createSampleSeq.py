##
# @file createSampleSeq.py
# @brief Creates and stores an array of sample query sequences
# @author Markus Schmidt
from .aligner import *
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.palettes import d3
from bokeh.models import LinearAxis, Range1d
from bokeh.models.formatters import FuncTickFormatter
from bokeh.models import LinearColorMapper, LogTicker, ColorBar
import sqlite3
import numpy
import random
from random import shuffle
from math import floor
import operator
import os
import csv 

def mutate(char):
    num = 0
    if char == "c" or char == "C":
        num = 1
    elif char == "t" or char == "T":
        num = 2
    elif char == "g" or char == "G":
        num = 3
    elif char == "a" or char == "A":
        num = 0
    else:
        print("error while mutating")
        return 'n'

    num += random.randint(1,3)
    num %= 4
    if num == 0:
        return 'a'
    elif num == 1:
        return 'c'
    elif num == 2:
        return 't'
    elif num == 3:
        return 'g'
    else:
        print("error while mutating")
        return 'n'

def setUpDbTables(conn, reset = False):
    c = conn.cursor()

    if reset:
        c.execute("DROP TABLE IF EXISTS samples")
        c.execute("DROP TABLE IF EXISTS results")

    c.execute("""
                CREATE TABLE IF NOT EXISTS samples
                (
                    sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    num_mutation INTEGER,
                    num_indels INTEGER,
                    original_size INTEGER,
                    origin INTEGER,
                    paired_origin INTEGER,
                    indel_size REAL,
                    sequence TEXT,
                    paired_sequence TEXT,
                    reference TINYTEXT
                )
                """)

    c.execute("""
                CREATE TABLE IF NOT EXISTS results
                (
                    sample_id INTEGER,
                    score REAL,
                    score2 REAL,
                    optimal_score_this_region REAL,
                    result_start INTEGER,
                    result_end INTEGER,
                    num_seeds INTEGER,
                    index_of_chosen_strip INTEGER,
                    seed_coverage_chosen_strip INTGER,
                    seed_coverage_alignment INTGER,
                    num_seeds_chosen_strip INTEGER,
                    anchor_size INTEGER,
                    anchor_ambiguity INTEGER,
                    max_diag_deviation INTEGER,
                    max_diag_deviation_percent INTEGER,
                    max_nmw_area INTEGER,
                    nmw_area INTEGER,
                    run_time REAL,
                    mapping_quality REAL,
                    approach TINYTEXT
                )
                """)

    c.execute("""
                CREATE INDEX IF NOT EXISTS results_index_1 ON results
                (
                    approach
                )
                """)

    c.execute("""
                CREATE INDEX IF NOT EXISTS results_index_2 ON results
                (
                    sample_id
                )
                """)

    c.execute("""
                CREATE INDEX IF NOT EXISTS samples_index_1 ON samples
                (
                    num_indels
                )
                """)

    c.execute("""
                CREATE INDEX IF NOT EXISTS samples_index_2 ON samples
                (
                    original_size
                )
                """)

    c.execute("""
                CREATE INDEX IF NOT EXISTS samples_index_3 ON samples
                (
                    num_mutation
                )
                """)

    c.execute("""
                CREATE INDEX IF NOT EXISTS samples_index_4 ON samples
                (
                    sample_id
                )
                """)
    conn.commit()

def insertQueries(conn, queries_list):
    c = conn.cursor()
    insert_str = """
        INSERT INTO samples 
        (
            num_mutation,
            num_indels,
            original_size,
            origin,
            paired_origin,
            indel_size,
            sequence,
            paired_sequence,
            reference
        ) 
        """

    if len(queries_list[0]) == 7:
        insert_str += "VALUES (?,?,?,?,-1,?,?,\"\",?)"
    else:
        insert_str += "VALUES (?,?,?,?,?,?,?,?,?)"
    c.executemany(insert_str, queries_list)
    conn.commit()

def getNewQueries(db_name, approach, reference, res_mut = 1, res_indel = 1, size = None, give_orig_pos = False, give_orig_size = False):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    select = "sequence, sample_id"
    if give_orig_pos:
        select += ", origin"
    if give_orig_size:
        select += ", original_size"
    if size is None:
        return c.execute("""
                        SELECT """ + select + """
                        FROM samples
                        WHERE sample_id NOT IN 
                            (
                                SELECT DISTINCT sample_id
                                FROM results
                                WHERE approach == ?
                            )
                        AND reference == ?
                        AND num_mutation % ? == 0
                        AND num_indels % ? == 0
                        """, (approach, reference, res_mut, res_indel)).fetchall()
    return c.execute("""
                        SELECT """ + select + """
                        FROM samples
                        WHERE sample_id NOT IN 
                            (
                                SELECT DISTINCT sample_id
                                FROM results
                                WHERE approach == ?
                            )
                        AND reference == ?
                        AND num_mutation % ? == 0
                        AND num_indels % ? == 0
                        AND original_size == ?
                        """, (approach, reference, res_mut, res_indel, size)).fetchall()

def getQueriesFor(db_name, reference, mut=1, indel=1, size=100):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT sequence, sample_id
                        FROM samples
                        WHERE reference == ?
                        AND num_mutation == ?
                        AND num_indels == ?
                        AND original_size == ?
                        """, (reference, mut, indel, size)).fetchall()

def getQueriesAsASDMatrix(db_name, approach, reference, give_seed_length = False):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    indel_amounts = c.execute("""
                        SELECT DISTINCT num_indels
                        FROM samples
                        WHERE reference == ?
                        ORDER BY num_indels
                        """, (reference,)).fetchall()
    mut_amounts = c.execute("""
                        SELECT DISTINCT num_mutation
                        FROM samples
                        WHERE reference == ?
                        ORDER BY num_mutation
                        """, (reference,)).fetchall()
    result = []
    for indel_amount in indel_amounts:
        result.append( [] )
        for mut_amount in mut_amounts:
            select = "sequence, sample_id, origin"
            if give_seed_length:
                select += ", original_size"
            elements = c.execute("""
                        SELECT """ + select + """
                        FROM samples
                        WHERE sample_id NOT IN 
                            (
                                SELECT DISTINCT sample_id
                                FROM results
                                WHERE approach == ?
                            )
                        AND reference == ?
                        AND num_mutation == ?
                        AND num_indels == ?
                        """, (approach, reference, mut_amount[0], indel_amount[0])).fetchall()
            if len(elements) == 0:
                continue
            result[-1].append(elements)
    return result



def getOriginOf(db_name, sample_id):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT origin
                        FROM samples
                        WHERE sample_id == ?
                        """, (sample_id,)).fetchone()[0]

def submitResults(db_name, results_list):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    if len(results_list[0]) == 8:
        c.executemany("""
                        INSERT INTO results 
                        (
                            sample_id,
                            score,
                            score2,
                            result_start,
                            result_end,
                            num_seeds,
                            index_of_chosen_strip,
                            seed_coverage_chosen_strip,
                            seed_coverage_alignment,
                            num_seeds_chosen_strip,
                            anchor_size,
                            anchor_ambiguity,
                            max_diag_deviation,
                            max_diag_deviation_percent,
                            max_nmw_area,
                            nmw_area,
                            run_time,
                            mapping_quality,
                            approach
                        )
                        VALUES (?,?,0,?,?,?,0,0,0,0,0,0,0,0,0,0,?,?,?)
                        """, results_list)
    else: # len == 20
        c.executemany("""
                        INSERT INTO results 
                        (
                            sample_id,
                            score,
                            score2,
                            optimal_score_this_region,
                            result_start,
                            result_end,
                            num_seeds,
                            index_of_chosen_strip,
                            seed_coverage_chosen_strip,
                            seed_coverage_alignment,
                            num_seeds_chosen_strip,
                            anchor_size,
                            anchor_ambiguity,
                            max_diag_deviation,
                            max_diag_deviation_percent,
                            max_nmw_area,
                            nmw_area,
                            mapping_quality,
                            run_time,
                            approach
                        )
                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                        """, results_list)
    conn.commit()

def clearResults(db_name, ref, approach):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("""
                DELETE FROM results 
                WHERE approach = ?
                """, (approach,))
    conn.commit()

def getApproachesWithData(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT DISTINCT approach
                        FROM results
                    """).fetchall()

##
# returns a list of tuples.
# the tuple consists of:
#
# ( score, result_start, original_start, num_mutation, num_indels, num_seeds, run_time )
#
# returns the results for the specified reference and approach only
# if no reference is specified the tuple contains the reference
def getResults(db_name, approach, size=None, indel_size=None, reference=None):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    if size is None and indel_size is None:
        return c.execute("""
                            SELECT 
                                results.score,
                                results.score2,
                                results.optimal_score_this_region,
                                results.result_start,
                                results.result_end,
                                samples.origin,
                                samples.num_mutation,
                                samples.num_indels,
                                results.num_seeds,
                                results.mapping_quality,
                                results.run_time,
                                results.index_of_chosen_strip,
                                results.seed_coverage_chosen_strip,
                                results.seed_coverage_alignment,
                                results.num_seeds_chosen_strip,
                                results.anchor_size,
                                results.anchor_ambiguity,
                                results.max_diag_deviation,
                                results.max_diag_deviation_percent,
                                results.max_nmw_area,
                                results.nmw_area,
                                samples.original_size
                            FROM samples
                            JOIN results ON results.sample_id = samples.sample_id
                            AND results.approach == ?
                            ORDER BY results.sample_id
                            """, (approach, )).fetchall()
    elif reference is None:
        return c.execute("""
                            SELECT 
                                results.score,
                                results.score2,
                                results.optimal_score_this_region,
                                results.result_start,
                                results.result_end,
                                samples.origin,
                                samples.num_mutation,
                                samples.num_indels,
                                results.num_seeds,
                                results.num_seeds_chosen_strip,
                                results.seed_coverage_chosen_strip,
                                results.seed_coverage_alignment,
                                results.mapping_quality,
                                results.nmw_area,
                                results.run_time,
                                samples.sequence,
                                results.sample_id
                            FROM samples
                            JOIN results ON results.sample_id = samples.sample_id
                            WHERE results.approach == ?
                            AND samples.original_size == ?
                            AND samples.indel_size == ?
                            ORDER BY results.sample_id
                            """, (approach, size, indel_size)).fetchall()
    return c.execute("""
                        SELECT 
                            results.score,
                            results.result_start,
                            results.result_end,
                            samples.origin,
                            samples.num_mutation,
                            samples.num_indels
                        FROM samples
                        JOIN results ON results.sample_id = samples.sample_id
                        AND results.approach == ?
                        WHERE samples.original_size == ?
                        AND samples.indel_size == ?
                        AND samples.reference == ?
                        """, (approach, size, indel_size, reference)).fetchall()

def near(start_align, start_orig, end_align, end_orig):
    return end_align >= start_orig and start_align <= end_orig

def analyzeAccuracy(db_name, out_file_name, approaches, res_mut, res_indel, size, 
        indel_size, reference=None):
    plots = []
    output_file(out_file_name)

    # analyze results
    max_indels = int(size/indel_size)*4
    for approach in approaches:
        # create empty matrix
        accurateMatrix = []
        samplesMatrix = []
        for _ in floor(size/res_mut):
            accurateMatrix.append( [] )
            samplesMatrix.append( [] )
            for _ in floor(max_indels/res_indel):
                accurateMatrix[-1].append( [] )
                samplesMatrix[-1].append( [] )

        # collect results
        results = getResults(db_name, approach, size, indel_size, reference)

        # fill in the matrix
        for score, result_start, result_end, original_start, num_mutation, num_indels in results:
            samplesMatrix[num_mutation][num_indels] += 1
            if near(original_start, result_start, original_start+size, result_end):
                accurateMatrix[num_mutation][num_indels] += 1

        # divide the matrices
        for x, row in enumerate(accurateMatrix):
            for y, ele in enumerate(row):
                accurateMatrix[x][y] = ele / samplesMatrix[x][y]

        # plot results
        color_mapper = LinearColorMapper(palette="Viridis256", low=0, high=1)
        plot = figure(title="accuracy " + approach,
                x_range=(0,max_indels), y_range=(0,size),
                x_axis_label='num indels', y_axis_label='num mutations'
            )
        plot.image(image=[accurateMatrix], color_mapper=color_mapper,
                dh=[size], dw=[max_indels], x=[0], y=[0])
        color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
        plot.add_layout(color_bar, 'left')
        plots.append(plot)

    show(row(plots))

def get_query(ref_seq, q_len, mutation_amount, indel_amount, indel_size, in_to_del_ratio=0.5):
    q = ""
    original_nuc_dist = [0, 0, 0, 0, 0]
    modified_nuc_dist = [0, 0, 0, 0, 0]

    #
    # get a non- bridging sequence - just randomly try until one is okay
    #

    q_from = 0
    #do - while loop
    while True:#do
        q_from = random.randint(0, ref_seq.unpacked_size()/2 - q_len)
        q_to = q_from + q_len
    #while
        if ref_seq.is_bridging(q_from, q_len):
            continue
        q = str(ref_seq.extract_from_to(q_from, q_to))
        break

    for nuc in q:
        if nuc == 'A' or nuc == 'a':
            original_nuc_dist[0] += 1
        elif nuc == 'C' or nuc == 'c':
            original_nuc_dist[1] += 1
        elif nuc == 'G' or nuc == 'g':
            original_nuc_dist[2] += 1
        elif nuc == 'T' or nuc == 't':
            original_nuc_dist[3] += 1
        else:
            original_nuc_dist[4] += 1

    #
    # apply modifications
    # deletions, mutations, insertions
    # the order is important.
    # code makes sure that the same nucleotide is never modified twice
    # also indels must be at least one nuc apart from each other...
    #

    ##
    # helper function
    # returns a list (amount elements) of random indices that are part of [0, interval_length]
    # and at least min_distance apart from each other
    #
    # if it is not possible to fit enough indices into [0, interval_length] the last indices 
    # are behind interval_length
    def get_random_spots(amount, interval_length, min_distance):
        spots = []
        for index in range(interval_length - min_distance):
            spots.append(index)
        shuffle(spots)
        spots = spots[:amount]
        spots.sort()
        for index in range(1,len(spots)):
            if spots[index-1] + min_distance + 1 > spots[index]:
                spots[index] = spots[index-1] + min_distance + 1
        if len(spots) >=1 and spots[-1] > interval_length:
            start = spots[0]
            for index in range(len(spots)):
                spots[index] -= start
        return spots

    deletion_amount = floor((1-in_to_del_ratio) * indel_amount)
    insertion_amount = floor( in_to_del_ratio*(indel_amount+1))

    # deletion
    deletion_spots = get_random_spots(deletion_amount, len(q), indel_size + 1)
    for pos in reversed(deletion_spots):
        q = q[:pos] + q[pos + indel_size:]

    # mutations
    mutation_spots = []
    for index in range(len(q)):
        mutation_spots.append(index)
    shuffle(mutation_spots)
    for pos in mutation_spots[:mutation_amount]:
        l = len(q)
        q = q[:pos] + mutate(q[pos]) + q[pos+1:]
        if not len(q) == l:
            print("ERROR: mutation changed query length:" + str(l) + " != " + str(len(q)))
            print(q)

    # insertion
    insertion_spots = get_random_spots(insertion_amount, len(q), 1)
    # pos are sorted in order to we need to reverse them in order 
    # to not insert twice at the same location
    for pos in reversed(insertion_spots):
        for _ in range(indel_size):
            char = random.randint(1,4)
            if char == 1:
                q = q[:pos] + "a" + q[pos:]
            elif char == 2:
                q = q[:pos] + "c" + q[pos:]
            elif char == 3:
                q = q[:pos] + "t" + q[pos:]
            elif char == 4:
                q = q[:pos] + "g" + q[pos:]
            else:
                print("ERROR: indel produced impossible nucleotide")

    for nuc in q:
        if nuc == 'A' or nuc == 'a':
            modified_nuc_dist[0] += 1
        elif nuc == 'C' or nuc == 'c':
            modified_nuc_dist[1] += 1
        elif nuc == 'G' or nuc == 'g':
            modified_nuc_dist[2] += 1
        elif nuc == 'T' or nuc == 't':
            modified_nuc_dist[3] += 1
        else:
            modified_nuc_dist[4] += 1

    return (q_from, q, original_nuc_dist, modified_nuc_dist)

def createSampleQueries(ref, db_name, size, indel_size, amount, reset = True, high_qual=False, in_to_del_ratio=0.5, smaller_box = False, only_first_row = False, validate_using_sw = True):
    conn = sqlite3.connect(db_name)

    setUpDbTables(conn, reset)

    ref_seq = Pack()
    ref_seq.load(ref)

    max_indels = 1
    if indel_size != 0:
        max_indels = int(size/indel_size)*2
    queries_list = []

    nuc_distrib_count_orig = [0,0,0,0,0]
    nuc_distrib_count_mod = [0,0,0,0,0]

    skip_x = max(1, int(size/50))
    skip_y = max(2, int(max_indels/50))

    if high_qual:
        skip_x = 1
        skip_y = 2
    
    if skip_y % 1 == 1:
        skip_y += 1

    max_x = int(size * 4 / 10)
    max_y = max_indels

    if smaller_box:
        max_x = int(max_x / 5)
        max_y = int(max_y / 5)

    if only_first_row:
        max_y = 1

    smw = SMW(False)
    genome = Pledge(NucSeq())
    if validate_using_sw:
        genome_seq = ref_seq.extract_complete()
        genome.set(genome_seq)

    #
    # iterate over the given range of mutations indels and number of sequences
    #
    for mutation_amount in range(0, max_x, skip_x):
        for indel_amount in range(0, max_y, skip_y):

            #dont put sequences that can't be found
            if mutation_amount + indel_amount/2 * indel_size >= size - 1:
                continue

            validated_queries = []
            while len(validated_queries) < amount:
                #
                # extract random query sequences and modify them
                #
                proposed_queries = []
                optimal_alignment_in = []
                while len(proposed_queries) < amount:
                    q_from, query, original_nuc_dist, modified_nuc_dist = get_query(
                        ref_seq, size, mutation_amount, indel_amount, indel_size, in_to_del_ratio)
                    if validate_using_sw:
                        p1 = Pledge(NucSeq())
                        p1.set(NucSeq(query))
                        optimal_alignment_in.append( (p1, genome) )
                    proposed_queries.append( (q_from, query, original_nuc_dist, modified_nuc_dist) )
                #
                # validate the query sequences using SW
                #
                if validate_using_sw:
                    optimal_alignment_out = []
                    for a, b in optimal_alignment_in:
                        optimal_alignment_out.append(smw.promise_me(a, b))

                    Pledge.simultaneous_get(optimal_alignment_out, 32)

                    num_discarded = 0
                    for index, alignments in enumerate(optimal_alignment_out):
                        if len(alignments.get()) > 1:
                            num_discarded += 1
                            continue
                        alignment = alignments.get()[0]
                        q_from = proposed_queries[index][0]
                        q_to = q_from + len(proposed_queries[index][1])
                        if not(q_from <= alignment.begin_on_ref and q_to >= alignment.begin_on_ref):
                            num_discarded += 1
                            continue
                        # if we arrive here the query is okay
                        validated_queries.append(proposed_queries[index])

                    if num_discarded > 0:
                        print("discarded", num_discarded, "queries due to SW.")
                        print("mutation amount:", mutation_amount, "indel amount:", indel_amount)

            #
            # append the validated queries
            #
            for q_from, query, original_nuc_dist, modified_nuc_dist in range(validated_queries):
                nuc_distrib_count_orig = list(map(operator.add, nuc_distrib_count_orig, original_nuc_dist))
                nuc_distrib_count_mod = list(map(operator.add, modified_nuc_dist, nuc_distrib_count_mod))

                #
                # construct the query tuple
                #
                queries_list.append((
                        mutation_amount,
                        indel_amount,
                        size,
                        q_from,
                        indel_size,
                        query,
                        ref
                    ))

                #
                # save the queries to the database
                #
                if len(queries_list) > 10000:
                    print("saving...")
                    insertQueries(conn, queries_list)
                    queries_list = []
                    print("done saving")
            #for _ in range(amount) end

        #for indel_amount in range(max_indels) end
        print(mutation_amount, "/", max_x)
    #for mutation_amount in range(size) end
    print("saving...")
    if len(queries_list) > 0:
        insertQueries(conn, queries_list)
    print("done saving")

    print("nuc distrib (A, C, G, T, N) originally: ", nuc_distrib_count_orig)
    print("nuc distrib (A, C, G, T, N) changed: ", nuc_distrib_count_mod)
    print("nuc distrib (A, C, G, T, N) changed by: ", list(map(operator.mul, map(
        operator.sub, nuc_distrib_count_mod, nuc_distrib_count_orig), [1./float(sum(nuc_distrib_count_orig))]*5)))
    print("total amount: ", sum(nuc_distrib_count_orig))
#function

def create_as_sequencer_reads(db_name, amount, technology="HS25", paired=False):
    print("setting up db...")
    conn = sqlite3.connect(db_name)
    setUpDbTables(conn, True)
    size = 0
    chrom = []
    for i in range(1,23):
        chrom.append("chr" + str(i) + ".fna")
    chrom.append("chrX.fna")
    chrom.append("chrY.fna")
    os.system("mkdir .temp_art")
    print("done")
    print("running simulator...")
    for in_file in chrom:
        print(in_file)
        command = "~/art_bin_MountRainier/"
        if technology == "HS25":
            size = 1000
            command += "art_illumina -ss HS25 -sam -i /mnt/ssd0/chrom/human/" + in_file + " -l " + str(size)
            if paired:
                command += " -p -m 200 -s 10"
            command += " -q -c " + str(amount) + " -o .temp_art/" + in_file

        os.system(command)
    print("done")
    print("extracting sequences...")
    genome = "/mnt/ssd0/genome/human"
    pack = Pack()
    pack.load(genome)

    origin = 0
    first = True
    sequence = ""
    sequences = []
    for in_file in chrom:
        print(in_file)
        with open(".temp_art/" + in_file + ".sam", "r") as csv_file:
            reader = csv.reader(csv_file, delimiter='\t')
            skip = 3
            for row in reader:
                if skip > 0:
                    skip -= 1
                    continue
                if paired:
                    if first:
                        origin = pack.start_of_sequence(row[2]) + int(row[3])
                        sequence = row[9]
                    else:
                        sequences.append( [0,0,size,origin,pack.start_of_sequence(row[2]) + int(row[3]) - 1,0,sequence,row[9],genome] )
                    first = not first
                else:
                    sequences.append( [0,0,size,pack.start_of_sequence(row[2]) + int(row[3]) - 1,0,row[9],genome] )

    os.system("rm -r .temp_art/")
    print("done")
    print("commiting...")
    insertQueries(conn, sequences)
    print("done")



db_name = "/mnt/ssd1/alignmentSamples.db"

#createSampleQueries("/mnt/ssd0/chrom/human/all", db_name, 1000, 100, 50, 64)