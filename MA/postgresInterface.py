import psycopg2
from .aligner import *

class CustomCursor():
    def __init__(self, con_str):
        self.conn = psycopg2.connect(con_str)
        self.cur = self.conn.cursor()
    def __enter__(self):
        return self.cur
    def __exit__(self, type, value, traceback):
        if type is None:
            self.conn.commit()
        else:
            print("Warning - did not commit due to exception")
        self.cur.close()
        self.conn.close()


class DbWriter(Module):
    def __init__(self, pack):
        self.conn = psycopg2.connect("host=192.168.1.10 port=5432 dbname=MA user=postgres")
        self.cur = self.conn.cursor()

        self.cur.execute(
                "INSERT INTO sam_header (genome_name) VALUES (%s) RETURNING id",
                ("testGenome", )
            )
        genome_id = self.cur.fetchone()
        self.cur.execute(
                "INSERT INTO run (aligner_name, header_id) VALUES (%s, %s) RETURNING id",
                ("testAligner", genome_id)
            )
        self.run_id = self.cur.fetchone()

    def get_input_type(self):
        return [Alignment(), NucSeq()]

    def get_output_type(self):
        return Nil()

    def execute(self, *input):
        align = input[0]
        #query = input[1][align.begin_on_query:align.end_on_query]
        #ref = input[2]

        length = align.end_on_ref - align.begin_on_ref

        self.cur.execute(
                "INSERT INTO alignment (cigar, position, mapping_quality, query_id, run_id, length) VALUES (%s, %s, %s, %s, %s, %s)",
                (align.cigarString(), align.begin_on_ref, align.mapping_quality, align.stats.name, self.run_id, length)
            )

    def finalize(self):
        if not (self.conn is None and self.curr is None):
            self.conn.commit()
            self.cur.close()
            self.conn.close()
            self.conn = None
            self.cur = None

    def __del__(self):
        self.finalize()





## with CustomCursor("host=192.168.1.10 port=5432 dbname=MA user=postgres") as cur:
##     cur.execute(
##             "INSERT INTO sam_header (genome_name) VALUES (%s) RETURNING id",
##             ("testGenome", )
##         )
##     genome_id = cur.fetchone()
##     cur.execute(
##             "INSERT INTO run (aligner_name, header_id) VALUES (%s, %s) RETURNING id",
##             ("testAligner", genome_id)
##         )
##     run_id = cur.fetchone()
##     cur.execute(
##             "INSERT INTO alignment (position, run_id) VALUES (%s, %s)",
##             (100, run_id)
##         )