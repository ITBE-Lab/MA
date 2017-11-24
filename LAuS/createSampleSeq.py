##
# @file createSampleSeq.py
# @brief Creates and stores an array of sample query sequences
# @author Markus Schmidt
from .aligner import *
import sqlite3
import numpy
import random
from math import floor

def mutate(char):
    num = 0
    if char == "C":
        num = 1
    elif char == "T":
        num = 2
    elif char == "G":
        num = 3
    num += random.randint(1,3)
    num %= 4
    if num == 0:
        return 'A'
    elif char == 1:
        return 'C'
    elif char == 2:
        return 'T'
    else:
        return 'G'

def setUpDbTables(conn):
    c = conn.cursor()

    c.execute("""
                CREATE TABLE IF NOT EXISTS samples
                (
                    sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    num_mutation INTEGER,
                    num_indels INTEGER,
                    original_size INTEGER,
                    origin INTEGER,
                    indel_size_expectation REAL,
                    indel_size_variance REAL,
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
                    result_size INTEGER,
                    approach TINYTEXT
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
                        indel_size_expectation,
                        indel_size_variance,
                        sequence,
                        reference
                    )
                    VALUES (?,?,?,?,?,?,?,?)
                    """, queries_list)
    conn.commit()

def getNewQueries(ref, approach, reference):
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
                        """, (approach,reference)).fetchall()

def submitResults(ref, results_list):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.executemany("""
                    INSERT INTO results 
                    (
                        sample_id,
                        score,
                        result_start,
                        result_size,
                        approach
                    )
                    VALUES (?,?,?,?,?)
                    """, results_list)
    conn.commit()

def clearResults(ref, approach):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("""
                DELETE FROM results 
                WHERE approach = ?
                """, (approach,))
    conn.commit()

def get_del_ins_size(indel_size_expectation, indel_size_variance):
    return int(numpy.random.normal(indel_size_expectation, indel_size_variance))

def createSampleQueries(ref, db_name, size, indel_size_expectation, indel_size_variance, amount):
    conn = sqlite3.connect(db_name)

    setUpDbTables(conn)

    ref_seq = Pack()
    ref_seq.load(ref)

    max_indels = int(size/indel_size_expectation)*4
    queries_list = []

    #
    # iterate over the given range of mutations indels and number of sequences
    #
    for mutation_amount in range(size):
        for indel_amount in range(max_indels):


            for _ in range(amount):
                #
                # extract the query sequence
                #

                #do - while loop
                while True:#do
                    q_from = random.randint(0, ref_seq.unpacked_size()/2 - size)
                    q_to = q_from + size
                #while
                    if not ref_seq.is_bridging(q_from, size):
                        break
                query = str(ref_seq.extract_from_to(q_from, q_to))
                #
                # modify the sequence
                #

                #mutations
                for _ in range( floor(mutation_amount) ):
                    pos = random.randint(0,len(query)-1)
                    query = query[:pos-1] + mutate(query[pos]) + query[pos:]
                #deletion
                for _ in range( floor( indel_amount/2 ) ):
                    del_ins_size = get_del_ins_size(indel_size_expectation, indel_size_variance)
                    if len(query) > del_ins_size + 2:
                        pos = random.randint(0,len(query)-del_ins_size)
                        query = query[:pos-1] + query[pos + del_ins_size:]
                    else:
                        query = ""
                        break
                #insertion
                for _ in range( floor( (indel_amount+1)/2 ) ):
                    pos = 0
                    if len(query) > 2:
                        pos = random.randint(1,len(query)-1)
                    for _ in range(get_del_ins_size(indel_size_expectation, indel_size_variance)):
                        char = random.randint(1,4)
                        if char == 1:
                            query = query[:pos] + "a" + query[pos:]
                        elif char == 2:
                            query = query[:pos] + "c" + query[pos:]
                        elif char == 3:
                            query = query[:pos] + "t" + query[pos:]
                        else:
                            query = query[:pos] + "g" + query[pos:]
                if query == "":
                    query = "A"
                #
                # construct the query tuple
                #

                queries_list.append((
                        mutation_amount,
                        indel_amount,
                        size,
                        q_from,
                        indel_size_expectation,
                        indel_size_variance,
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