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
        c.execute("DROP TABLE IF EXISTS samples_optima")
        c.execute("DROP TABLE IF EXISTS results")
        c.execute("DROP TABLE IF EXISTS runtimes")

    c.execute("""
                CREATE TABLE IF NOT EXISTS samples
                (
                    sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    num_mutation INTEGER,
                    num_indels INTEGER,
                    original_size INTEGER,
                    origin INTEGER,
                    sequence TEXT
                )
                """)

    c.execute("""
                CREATE TABLE IF NOT EXISTS samples_optima
                (
                    sample_id INTEGER,
                    optima INTEGER
                )
                """)

    c.execute("""
                CREATE TABLE IF NOT EXISTS results
                (
                    sample_id INTEGER,
                    start INTEGER,
                    end INTEGER,
                    run_time REAL,
                    mapping_quality REAL,
                    approach TINYTEXT,
                    secondary TINYINT
                )
                """)

    c.execute("""
                CREATE TABLE IF NOT EXISTS runtimes
                (
                    num_mutation INTEGER,
                    num_indels INTEGER,
                    run_time REAL,
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
                CREATE INDEX IF NOT EXISTS samples_index_1 ON runtimes
                (
                    num_indels,
                    num_mutation
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

    c.execute("""
                CREATE INDEX IF NOT EXISTS samples_optima_index_1 ON samples_optima
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
            sequence
        ) 
        VALUES (?,?,?,?,?)
        """
    c.executemany(insert_str, queries_list)
    conn.commit()

def getQueries(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT sequence, sample_id
                        FROM samples
                        """).fetchall()


def getIndelsAndMuts(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    indel_amounts = c.execute("""
                        SELECT DISTINCT num_indels
                        FROM samples
                        ORDER BY num_indels
                        """).fetchall()
    mut_amounts = c.execute("""
                        SELECT DISTINCT num_mutation
                        FROM samples
                        ORDER BY num_mutation
                        """).fetchall()

    return indel_amounts, mut_amounts

def getQueriesAsASDMatrix(db_name):
    indel_amounts, mut_amounts = getIndelsAndMuts(db_name)
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    result = []
    for indel_amount in indel_amounts:
        result.append( [] )
        for mut_amount in mut_amounts:
            elements = c.execute("""
                        SELECT sequence, sample_id
                        FROM samples
                        AND num_mutation == ?
                        AND num_indels == ?
                        """, (mut_amount[0], indel_amount[0])).fetchall()
            if len(elements) == 0:
                continue
            result[-1].append(elements)
    return result

def getNumQueriesDict(db_name, reference):
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
    result = {}
    for mut_amount in mut_amounts:
        result[mut_amount[0]] = {}
        for indel_amount in indel_amounts:
            elements = c.execute("""
                        SELECT DISTINCT sample_id
                        FROM samples
                        WHERE reference == ?
                        AND num_mutation == ?
                        AND num_indels == ?
                        """, (reference, mut_amount[0], indel_amount[0])).fetchall()
            result[mut_amount[0]][indel_amount[0]] = len(elements)
    return result

def submitOptima(db_name, optima_list):
    results_list = []
    # convert the list
    for sample_id, score, positions in optima_list:
        for position in positions:
            results_list.append( (sample_id, position, score) )

    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("""DELETE FROM samples_optima""")
    c.executemany("""
                    INSERT INTO samples_optima 
                    (
                        sample_id,
                        optima
                    )
                    VALUES (?,?)
                    """, results_list)
    conn.commit()

def submitResults(db_name, results_list):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("""
                DELETE FROM results 
                WHERE approach = ?
                """, (results_list[-2],))
    c.executemany("""
                    INSERT INTO results 
                    (
                        sample_id,
                        start,
                        end,
                        mapping_quality,
                        approach,
                        secondary
                    )
                    VALUES (?,?,?,?,?,?,?)
                    """, results_list)
    conn.commit()

def submitRuntimes(db_name, results_list):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("""
                DELETE FROM runtimes 
                WHERE approach = ?
                """, (approach,))
    c.executemany("""
                    INSERT INTO runtimes 
                    (
                        num_mutation,
                        num_indels,
                        run_time,
                        approach
                    )
                    VALUES (?,?,?,?)
                    """, results_list)
    conn.commit()

def getApproachesWithData(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT DISTINCT approach
                        FROM results
                    """).fetchall()

def getOptimalPositions(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT 
                            samples.sample_id,
                            samples_optima.optima,
                            samples.original_size,
                            samples.num_mutation,
                            samples.num_indels
                        FROM samples
                        JOIN samples_optima ON samples.sample_id = samples_optima.sample_id
                    """).fetchall()

def getResults(db_name, approach):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT 
                            samples.sample_id,
                            results.secondary,
                            results.start,
                            results.end,
                            results.mapping_quality
                        FROM samples
                        JOIN results ON samples.sample_id = results.sample_id
                    """).fetchall()

def getRuntimes(db_name, approach):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    return c.execute("""
                        SELECT 
                            num_mutation,
                            num_indels,
                            run_time
                        FROM runtimes
                        WHERE approach = ?
                    """, (approach,) ).fetchall()

def near(start_align, start_orig, end_align, end_orig):
    return end_align >= start_orig and start_align <= end_orig

def getAccuracyAndRuntimeOfAligner(db_name, approach, max_tries):
    # get the indel and mut amounts so that we know the size of the resulting matrix
    indel_amounts, mut_amounts = getIndelsAndMuts(db_name)
    hits = [[0]*len(indel_amounts)]*len(mut_amounts)
    samples = [[0]*len(indel_amounts)]*len(mut_amounts)

    # store the aligner results in a dict of sample id
    aligner_results = {}
    for sample_id, secondary, start, end, mapping_quality in getResults(db_name, approach):
        if not sample_id in aligner_results:
            aligner_results[sample_id] = []
        if secondary: # append secondary ones
            aligner_results[sample_id].append( (start, end, mapping_quality) )
        else: # prepend the primary alignment so that we always use that one first
            aligner_results[sample_id] = [ (start, end, mapping_quality) ] + aligner_results[sample_id]

    # store the optimal positions in a dict of sample id
    optimal_pos = {}
    for sample_id, optima_end, original_size, num_mutation, num_indels in getOptimalPositions(db_name):
        if not sample_id in optimal_pos:
            optimal_pos[sample_id] = ([], num_mutation, num_indels)
        optimal_pos[sample_id].append( (optima_end-original_size, optima_end) )

    # iterate over all samples
    for sample_id, value in optimal_pos.items():
        positions, num_mutation, num_indels = value
        # record the sample amount
        samples[num_indels][num_mutation] += 1

        # figure out if within the first max_tries tries of the aligner 
        # there was a hit to one optimal position
        hit = False
        for start, end, mapping_quality in aligner_results[sample_id][:max_tries]:
            for start_orig, end_orig, in positions:
                if near(start, start_orig, end, end_orig):
                    hit = True
                    break
            if hit:
                break
        if hit:
            hits[num_indels][num_mutation] += 1

    # compute the accuracy
    accuracy = []
    for hit_row, sample_row in zip(hits, samples):
        accuracy.append( [] )
        for hit, sample in zip(hit_row, sample_row):
            if sample == 0: # if there were no samples we want to return NAN
                accuracy[-1].append( float('nan') )
            else: # compute the accuracy for that cell
                accuracy[-1].append( hit/sample )

    runtime = [[float('nan')]*len(indel_amounts)]*len(mut_amounts)
    for num_mutation, num_indels, run_time in getRuntimes(db_name, approach):
        runtime[num_indels][num_mutation] = run_time

    return accuracy, runtime

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
    # revcomp, deletions, mutations, insertions
    # the order is important.
    # code makes sure that the same nucleotide is never modified twice
    # also indels must be at least one nuc apart from each other...
    #

    #reverse complement
    if random.randint(0,1) == 1:
        comp = {
            'A' : 'T',
            'T' : 'A',

            'G' : 'C',
            'C' : 'G'
        }# dict
        q_ = ""
        for nuc in reversed(q):
            q_ += comp[nuc.upper()]
        q = q_



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
            command += "art_illumina -ss HS25 -sam -i /MAdata/chrom/human/" + in_file + " -l " + str(size)
            if paired:
                command += " -p -m 200 -s 10"
            command += " -q -c " + str(amount) + " -o .temp_art/" + in_file

        os.system(command)
    print("done")
    print("extracting sequences...")
    genome = "/MAdata/genome/human"
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



db_name = "/MAdata/alignmentSamples.db"

#createSampleQueries("/MAdata/chrom/human/all", db_name, 1000, 100, 50, 64)