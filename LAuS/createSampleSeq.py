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

def mutate(char):
    num = 0
    if char == "c":
        num = 1
    elif char == "t":
        num = 2
    elif char == "g":
        num = 3
    num += random.randint(1,3)
    num %= 4
    if num == 0:
        return 'a'
    elif char == 1:
        return 'c'
    elif char == 2:
        return 't'
    else:
        return 'g'

def setUpDbTables(conn):
    c = conn.cursor()

    c.execute("""
                CREATE TABLE IF NOT EXISTS samples
                (
                    sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    num_mutation INTEGER,
                    num_insertions INTEGER,
                    num_indels INTEGER,
                    original_size INTEGER,
                    origin INTEGER,
                    indel_size REAL,
                    sequence TEXT,
                    reference TINYTEXT
                )
                """)

    c.execute("""
                CREATE TABLE IF NOT EXISTS results
                (
                    sample_id INTEGER,
                    score REAL,
                    result_start INTEGER,
                    approach TINYTEXT
                )
                """)

    c.execute("""
                CREATE INDEX IF NOT EXISTS results_index_1 ON results
                (
                    approach
                )
                """)
    conn.commit()

def insertQueries(conn, queries_list):
    c = conn.cursor()
    c.executemany("""
                    INSERT INTO samples 
                    (
                        num_mutation,
                        num_indels,
                        original_size,
                        origin,
                        indel_size,
                        sequence,
                        reference
                    )
                    VALUES (?,?,?,?,?,?,?)
                    """, queries_list)
    conn.commit()

def getNewQueries(db_name, approach, reference, res_mut, res_indel):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT sequence, sample_id
                        FROM samples
                        WHERE sample_id NOT IN 
                            (
                                SELECT DISTINCT sample_id
                                FROM results
                                WHERE approach = ?
                                AND reference = ?
                            )
                        AND num_mutation % ? == 0
                        AND num_indels % ? == 0
                        """, (approach, reference, res_mut, res_indel)).fetchall()

def getNewQueries(db_name, approach, reference):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT sequence, sample_id
                        FROM samples
                        WHERE sample_id NOT IN
                            (
                                SELECT DISTINCT sample_id
                                FROM results
                                WHERE approach = ?
                                AND reference = ?
                            )
                        """, (approach, reference)).fetchall()

def submitResults(db_name, results_list):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.executemany("""
                    INSERT INTO results 
                    (
                        sample_id,
                        score,
                        result_start,
                        num_seeds,
                        run_time,
                        approach
                    )
                    VALUES (?,?,?,?,?,?)
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

##
# returns a list of tuples.
# the tuple consists of:
#
# ( score, result_start, original_start, num_mutation, num_indels, num_seeds, run_time )
#
# returns the results for the specified reference and approach only
# if no reference is specified the tuple contains the reference
def getResults(db_name, approach, size, indel_size, reference=None):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    if reference is None:
        return c.execute("""
                            SELECT 
                                (
                                    results.score,
                                    results.result_start,
                                    samples.origin,
                                    samples.num_mutation,
                                    samples.num_indels,
                                    samples.num_seeds,
                                    samples.run_time
                                )
                            FROM samples
                            JOIN results ON results.sample_id=samples.sample_id
                            WHERE results.approach == ?
                            AND samples.original_size == ?
                            AND samples.indel_size == ?
                            """, (approach, size, indel_size)).fetchall()
    return c.execute("""
                        SELECT 
                            (
                                results.score,
                                results.result_start,
                                samples.origin,
                                samples.num_mutation,
                                samples.num_indels
                            )
                        FROM samples
                        JOIN results ON results.sample_id=samples.sample_id
                        WHERE results.approach == ?
                        AND samples.original_size == ?
                        AND samples.indel_size == ?
                        AND samples.reference == ?
                        """, (approach, size, indel_size, reference)).fetchall()

min_accuracy = 100
def near(index, index_2):
    return index + min_accuracy > index_2 and index - min_accuracy < index_2

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
        for score, result_start, original_start, num_mutation, num_indels in results:
            samplesMatrix[num_mutation][num_indels] += 1
            if near(original_start, result_start):
                accurateMatrix[num_mutation][num_indels] += 1

        # divide the matrices
        for x, row in enumerate(accurateMatrix):
            for y, ele in enumerate(row):
                accurateMatrix[x][y] = ele / samplesMatrix[x][y]

        # plot results
        color_mapper = LinearColorMapper(palette="Viridis256", low=0, high=1)
        plot = figure(title="accuracy " + aligner_desc[index],
                x_range=(0,max_indels), y_range=(0,size),
                x_axis_label='num indels', y_axis_label='num mutations'
            )
        plot.image(image=[accurateMatrix], color_mapper=color_mapper,
                dh=[size], dw=[max_indels], x=[0], y=[0])
        color_bar = ColorBar(color_mapper=color_mapper, border_line_color=None, location=(0,0))
        plot.add_layout(color_bar, 'left')
        plots.append(plot)

    show(row(plots))

def get_query(ref_seq, q_len, mutation_amount, indel_amount, indel_size):
    q = ""

    #
    # get a non- bridging sequence - just randomly try until one is okay
    #

    q_from = 0
    #do - while loop
    while True:#do
        q_from = random.randint(0, ref_seq.unpacked_size() - q_len)
        q_to = q_from + q_len
    #while
        if ref_seq.is_bridging(q_from, q_len):
            continue
        q = str(ref_seq.extract_from_to(q_from, q_to))
        break

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
        spots = spots[:deletion_amount]
        spots.sort()
        for index in range(1,len(spots)):
            if spots[index-1] + min_distance + 1 > spots[index]:
                spots[index] = spots[index-1] + min_distance + 1
        if len(spots) >=1 and spots[-1] > interval_length:
            start = spots[0]
            for index in range(len(spots)):
                spots[index] -= start
        return spots

    deletion_amount = floor(indel_amount/2)
    insertion_amount = floor( (indel_amount+1) /2)

    #deletion
    deletion_spots = get_random_spots(deletion_amount, len(q), indel_size + 1)
    for pos in reversed(deletion_spots):
        q = q[:pos] + q[pos + indel_size:]

    #mutations
    mutation_spots = []
    for index in range(len(q)):
        mutation_spots.append(index)
    shuffle(mutation_spots)
    for pos in mutation_spots[:mutation_amount]:
        q = q[:pos-1] + mutate(q[pos]) + q[pos:]

    #insertion
    insertion_spots = get_random_spots(insertion_amount, len(q), 1)
    for pos in insertion_spots:
        for _ in range(indel_size):
            char = random.randint(1,4)
            if char == 1:
                q = q[:pos] + "a" + q[pos:]
            elif char == 2:
                q = q[:pos] + "c" + q[pos:]
            elif char == 3:
                q = q[:pos] + "t" + q[pos:]
            else:
                q = q[:pos] + "g" + q[pos:]

    return (q_from, q)

def createSampleQueries(ref, db_name, size, indel_size, amount):
    conn = sqlite3.connect(db_name)

    setUpDbTables(conn)

    ref_seq = Pack()
    ref_seq.load(ref)

    max_indels = int(size/indel_size)
    queries_list = []

    #
    # iterate over the given range of mutations indels and number of sequences
    #
    for mutation_amount in range(0, size, 10):
        for indel_amount in range(0, max_indels, 1):

            for _ in range(amount):
                #
                # extract the query sequence
                #
                q_from, query = get_query(ref_seq, size, mutation_amount, indel_amount, indel_size)

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
        print(str(mutation_amount) + "/" + str(size))
    #for mutation_amount in range(size) end
    print("saving...")
    if len(queries_list) > 0:
        insertQueries(conn, queries_list)
    print("done saving")
#function

db_name = "/mnt/ssd1/alignmentSamples.db"

#createSampleQueries("/mnt/ssd0/chrom/human/all", db_name, 1000, 100, 50, 64)