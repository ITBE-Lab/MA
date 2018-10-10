import psycopg2
from .aligner import *
import pysam
import os

## class CustomCursor():
##     def __init__(self, con_str):
##         self.conn = psycopg2.connect(con_str)
##         self.cur = self.conn.cursor()
##     def __enter__(self):
##         return self.cur
##     def __exit__(self, type, value, traceback):
##         if type is None:
##             self.conn.commit()
##         else:
##             print("Warning - did not commit due to exception")
##         self.cur.close()
##         self.conn.close()
##
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


class PyDbWriter(Module):
    def __init__(self):
        #192.168.1.10
        self.conn = psycopg2.connect(
            "host=192.168.1.10 port=5432 dbname=MA user=postgres")
        self.cur = self.conn.cursor()

        genome_id = 0  # used dummy
        ## @todo
        ## self.cur.execute(
        ##         "SELECT id FROM sam_header WHERE genome_name = %s",
        ##         (pack_name, )
        ##     )
        ## genome_id = self.cur.fetchone()
        ## if genome_id is None:
        ##     self.cur.execute(
        ##             "INSERT INTO sam_header (genome_name) VALUES (%s) RETURNING id",
        ##             ("testGenome", )
        ##         )
        ##     genome_id = self.cur.fetchone()
        self.cur.execute(
            "INSERT INTO run (aligner_name, header_id) VALUES (%s, %s) RETURNING id",
            ("MA-py", genome_id))
        self.run_id = self.cur.fetchone()

    def execute(self, *input):
        align = input[0]
        query = input[1]
        pack = input[2]

        length = align.end_on_ref - align.begin_on_ref

        self.cur.execute(
            "INSERT INTO alignment (cigar, position, mapping_quality, query_id, run_id, sam_flags, contig, query_sequence, length) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)",
            (align.cigarString(pack), align.getSamPosition(pack),
             align.mapping_quality, align.stats.name, self.run_id,
             align.getSamFlag(pack), align.getContig(pack),
             align.getQuerySequence(query, pack), length))

    def finalize(self, runtime=None):
        if not (self.conn is None and self.cur is None):
            if not runtime is None:
                self.cur.execute("UPDATE run SET runtime=%s WHERE id=%s",
                                 (runtime, self.run_id))
            self.conn.commit()
            self.cur.close()
            self.conn.close()
            self.conn = None
            self.cur = None


def get_alignments_from_db(run_id, contig, start, end):
    conn = psycopg2.connect(
        "host=192.168.1.10 port=5432 dbname=MA user=postgres")
    cur = conn.cursor()

    cur.execute(
        """
            SELECT cigar, position, mapping_quality, query_id, sam_flags, query_sequence, length,
            contig_other, position_other
            FROM alignment
            WHERE run_id = %s
            AND contig = %s
            AND position + length > %s
            AND position < %s
            ORDER BY position
            """, (run_id, contig, start, end))

    ret = cur.fetchall()

    cur.close()
    conn.close()

    return ret


def read_cigar(cigar):
    symbol_to_operation = {
        "I": 1,
        "D": 2,
        "=": 7,
        "X": 8,
    }
    if cigar[0] == '*':
        return
    number_start = 0
    q_len = 0
    while number_start < len(cigar):
        symbol_start = number_start
        #increase symbol_start while cigar[symbol_start] is a number...
        while cigar[symbol_start] in [str(x) for x in range(10)]:
            symbol_start += 1
        num = int(cigar[number_start:symbol_start])
        sym = cigar[symbol_start]
        if sym in ["=", "X", "I"]:
            q_len += num
        yield num, symbol_to_operation[sym], q_len
        number_start = symbol_start + 1


def reads_to_bam(run_id, contig, start, end, pack, bam_file_name):
    alignment_list = get_alignments_from_db(run_id, contig, start, end)
    print("extracted", len(alignment_list), "alignments")
    ref_name_list = []
    contig_id = None
    for index, name in enumerate(pack.contigNames()):
        ref_name_list.append(name)
        if name == contig:
            contig_id = index
    ref_len_list = []
    for l in pack.contigLengths():
        ref_len_list.append(l)
    with pysam.AlignmentFile(
            bam_file_name,
            'wb',
            reference_names=ref_name_list,
            reference_lengths=ref_len_list) as bam_file:
        for cigar, position, mapping_quality, query_id, sam_flags, query_sequence, length_, contig_other, position_other in alignment_list:
            pys_alignment = pysam.AlignedSegment()
            pys_alignment.query_name = query_id
            pys_alignment.query_sequence = query_sequence
            pys_alignment.flag = sam_flags
            #if sam_flags | 0x080 == 0:
            #    continue
            pys_alignment.reference_id = contig_id
            pys_alignment.reference_start = position
            pys_alignment.next_reference_id = contig_id
            if contig_other != "*":
                for index, name in enumerate(pack.contigNames()):
                    if name == contig_other:
                        pys_alignment.next_reference_id = index
                pys_alignment.next_reference_start = position_other
            #pys_alignment.cigarstring = cigar
            cigar_list = []
            length = 0
            for amount, operation, length in read_cigar(cigar):
                cigar_list.append((operation, amount))
            pys_alignment.cigar = cigar_list
            pys_alignment.mapping_quality = max(
                min(int(mapping_quality * 254), 254), 0)
            pys_alignment.template_length = length  #len(query_sequence)
            bam_file.write(pys_alignment)


def reads_to_bam_w_bai(run_id, contig, start, end, pack, prefix):
    # allow abbreviations for the chromosomes
    name_trans_dict = {
        "chr1": "CM000663.2",
        "chr2": "CM000664.2",
        "chr3": "CM000665.2",
        "chr4": "CM000666.2",
        "chr5": "CM000667.2",
        "chr6": "CM000668.2",
        "chr7": "CM000669.2",
        "chr8": "CM000670.2",
        "chr9": "CM000671.2",
        "chr10": "CM000672.2",
        "chr11": "CM000673.2",
        "chr12": "CM000674.2",
        "chr13": "CM000675.2",
        "chr14": "CM000676.2",
        "chr15": "CM000677.2",
        "chr16": "CM000678.2",
        "chr17": "CM000679.2",
        "chr18": "CM000680.2",
        "chr19": "CM000681.2",
        "chr20": "CM000682.2",
        "chr21": "CM000683.2",
        "chr22": "CM000684.2",
        "chrX": "CM000685.2",
        "chrY": "CM000686.2"
    }
    if contig in name_trans_dict:
        contig = name_trans_dict[contig]

    reads_to_bam(run_id, contig, start, end, pack, prefix + ".bam")
    sam_tools_pref = "~/workspace/samtools/samtools "
    index_cmd = sam_tools_pref + "index " + prefix + ".bam > " + prefix + ".bam.bai"
    os.system(index_cmd)
