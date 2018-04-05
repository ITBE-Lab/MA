import os
import subprocess
from MA import *
import sys
import time


class CommandLine(Module):

    def __init__(self):
        self.elapsed_time = 0
        self.check_existence = []
        self.skip_count = 0
        self.warn_once = True

    def check(self):
        okay = True
        for index, existence in enumerate(self.check_existence):
            if existence > 1:
                print("Error: index", index, "had", existence, "alignments associated")
                okay = False
            if existence == 0:
                self.skip_count += 1
        return okay

    def final_checks(self):
        if not self.do_checks():
            return
        if self.skip_count > 0:
            print("Aligner skipped", self.skip_count, "queries")

    def __get_sam(self, index_str, queries):
        f = open(self.in_filename, "w")
        for index, query in enumerate(queries):
            check = ""
            f.write(">" + str(index) + " sequence" + str(index) + "\n" )
            query_string = str(query)
            assert(len(query) > 0)
            assert(len(query_string) > 0)
            while len(query_string) >= 60:
                line = query_string[:60]
                query_string = query_string[60:]
                f.write(line + "\n")
                check += line
            f.write(query_string + "\n")
            check += query_string
            if not check == str(query):
                print(check, "!=", str(query))
                assert(False)
        f.close()

        #assemble the shell command
        cmd_str = self.create_command(self.in_filename)
        #print(cmd_str)
        #exit()

        start_time = time.time()
        result = subprocess.run(cmd_str, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
        self.elapsed_time = time.time() - start_time

        if result.returncode != 0:
            print("subprocess returned with ERROR:")
            print("Error:")
            print(result.stderr.decode('utf-8')[:-1])
            exit()

        os.remove(self.in_filename)

        sam_file = result.stdout.decode('utf-8')

        return sam_file

    def __align(self, index_str, queries, pack):
        sam = self.__get_sam(index_str, queries)
        #print(sam)

        lines = sam.split("\n")
        while len(lines) > 0 and ( len(lines[0]) == 0 or lines[0][0] is '@' ):
            lines = lines[1:]

        #transform sam file into list data structure
        alignments = []
        #check how often each query was aligned
        del self.check_existence[:]

        for _ in range(len(queries)):
            self.check_existence.append(0)

        for line in lines:
            if len(line) == 0:
                continue
            columns = line.split("\t")
            #print(columns[0], columns[2], columns[3], columns[4])
            try:
                #print(str(int(columns[3])) + " " + str(pack.start_of_sequence(columns[2])))
                align_length = 0
                ##
                # brief helper function to read cigars...
                def read_cigar(cigar):
                    if cigar[0] == '*':
                        return
                    number_start = 0
                    while number_start < len(cigar):
                        symbol_start = number_start
                        #increase symbol_start while cigar[symbol_start] is a number...
                        while cigar[symbol_start] in [str(x) for x in range(10)]:
                            symbol_start += 1
                        yield int(cigar[number_start:symbol_start]), cigar[symbol_start]
                        number_start = symbol_start + 1

                assert(len(columns) >= 5)
                for amount, char in read_cigar(columns[5]):
                    if char in ['M', 'X', '=', 'D', 'N']:
                        align_length += amount
                    elif not char in ['I', 'S', 'H', 'P']:#sanity check...
                        print("Error: got wierd cigar symbol", char, "in cigar", columns[5])
                        exit()
                start = pack.start_of_sequence(columns[2]) + int(columns[3])
                alignment = Alignment(start, start + align_length)
                alignment.mapping_quality = int(columns[4])/255
                alignment.stats.name = columns[0]
                # set the secondary flag for secondary alignments
                alignment.secondary = True if int(columns[1]) | 0x100 > 0 else False
                if not alignment.secondary:
                    self.check_existence[int(columns[0])] += 1
                alignments.append(alignment)
            except Exception as e:
                print("Error:", e)
                print(line)
                traceback.print_exc()
                pass

        #transform list into alignment data structure
        ret = ContainerVector(Alignment())
        del ret[:]

        for alignment in alignments:
            ret.append(alignment)

        if self.do_checks():
            if not self.check():
                return None

        return ret

    def get_input_type(self):
        return [ContainerVector(NucSeq()), Pack()]

    def get_output_type(self):
        return ContainerVector(Alignment())

    def execute(self, *input):
        #print(input)
        queries = input[0]
        pack = input[1]

        return self.__align(self.index_str, queries, pack)

class Bowtie2(CommandLine):
    def __init__(self, index_str, threads, num_results, db_name):
        super().__init__()
        self.bowtie2_home = "/usr/home/markus/workspace/bowtie2/bowtie2-2.3.3.1/"
        self.index_str = index_str + "bowtie2"
        self.threads = threads
        self.num_results = num_results
        self.in_filename = ".tempBowtie2" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.bowtie2_home + "bowtie2 -p " + str(self.threads)
        index_str = "-x " + self.index_str
        input_str = "-f -U " + in_filename
        return cmd_str + " " + index_str + " " + input_str + " -k " + self.num_results

    def do_checks(self):
        return False

class Minimap2(CommandLine):
    def __init__(self, index_str, threads, num_results, db_name):
        super().__init__()
        self.minimap2_home = "/usr/home/markus/workspace/minimap2/"
        self.index_str = index_str + ".mmi"
        self.threads = threads
        self.num_results = num_results
        self.in_filename = ".tempMinimap2" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.minimap2_home + "minimap2 -t " + str(self.threads) + " -a "
        return cmd_str + " " + self.index_str + " " + in_filename + " --secondary=yes -N " + self.num_results

    def do_checks(self):
        return False

class Blasr(CommandLine):
    def __init__(self, index_str, threads, num_results, genome_str, db_name):
        super().__init__()
        self.blasr_home = "/usr/home/markus/workspace/blasr/build/bin/"
        self.index_str = index_str + "blasr"
        self.genome_str = genome_str
        self.threads = threads
        self.num_results = num_results
        self.in_filename = ".tempBlasr" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.blasr_home + "blasr " + in_filename
        return cmd_str + " " + self.genome_str + " --sam --bestn " + self.num_results + " --nproc " + str(self.threads) + " --hitPolicy leftmost --sa " + self.index_str

    def do_checks(self):
        return False

class BWA_MEM(CommandLine):
    def __init__(self, index_str, threads, num_results, db_name):
        super().__init__()
        self.bwa_home = "/usr/home/markus/workspace/bwa/"
        self.index_str = index_str + "bwa"
        self.threads = threads
        self.num_results = num_results
        self.in_filename = ".tempBwaMem" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.bwa_home + "bwa mem -t " + str(self.threads)
        if int(self.num_results) > 1:
            cmd_str += " -a"
        return cmd_str + " " + self.index_str + " " + in_filename

    def do_checks(self):
        return False

class BWA_SW(CommandLine):
    def __init__(self, index_str, threads, num_results, db_name):
        super().__init__()
        self.bwa_home = "/usr/home/markus/workspace/bwa/"
        self.index_str = index_str + "bwa"
        self.threads = threads
        self.num_results = num_results
        self.in_filename = ".tempBwaSw" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.bwa_home + "bwa bwasw -t " + str(self.threads)
        return cmd_str + " " + self.index_str + " " + in_filename

    def do_checks(self):
        return False


class G_MAP(CommandLine):
    def __init__(self, index_str, threads, num_results, genome_str, db_name):
        super().__init__()
        self.g_home = "/usr/home/markus/workspace/graphmap/bin/Linux-x64/"
        self.genome_str = genome_str
        self.threads = threads
        self.num_results = num_results
        self.index_str = index_str
        self.in_filename = ".tempBlasr" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.g_home + "graphmap align -r " + self.genome_str
        return cmd_str + " -d " + in_filename + " -t " + str(self.threads) + " -Z"

    def do_checks(self):
        return False

class MA(CommandLine):
    def __init__(self, index_str, threads, num_results, fast, db_name, num_soc=None):
        super().__init__()
        self.ma_home = "/usr/home/markus/workspace/aligner/"
        self.index_str = index_str
        self.threads = threads
        self.num_results = num_results
        self.fast = "accurate"
        self.num_soc = num_soc
        if fast:
            self.fast = "fast"
        self.in_filename = ".tempMA" + self.fast + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.ma_home + "ma -a -t " + str(self.threads) + " -p " + self.fast
        if not self.num_soc is None:
            cmd_str += " -G -S " + str(self.num_soc)
        return cmd_str + " -g " + self.index_str + " -i " + in_filename + " -n " + self.num_results

    def do_checks(self):
        return True

human_genome = "/MAdata/genome/human"

def test(
            db_name,
            reference,
            no_time=True
        ):
    print("working on " + db_name)
    ref_pack = Pack()
    ref_pack.load(reference)
    reference_pledge = Pledge(Pack())
    reference_pledge.set(ref_pack)

    num_threads = 1
    num_results = "3"

    l = [
        #("MINIMAP 2", Minimap2(reference, num_threads, num_results, db_name)),
        #("BLASR", Blasr(reference, num_threads, num_results, "/MAdata/chrom/human/n_free.fasta", db_name)),
        #("MA Fast", MA(reference, num_threads, num_results, True, db_name)),
        #("MA Accurate", MA(reference, num_threads, num_results, False, db_name)),
        ("BWA MEM", BWA_MEM(reference, num_threads, num_results, db_name)),
        #("BWA SW", BWA_SW(reference, num_threads, num_results, db_name)),
        #("BOWTIE 2", Bowtie2(reference, num_threads, num_results, db_name)),
        #("GRAPH MAP", G_MAP(reference, num_threads, num_results, "/MAdata/chrom/human/n_free.fasta", db_name)),
    ]

    for name, aligner in l:
        print("evaluating " + name)
        clearResults("/MAdata/db/"+db_name, reference, name) # CAREFUL WITH THE CLEARING

        matrix = getQueriesAsASDMatrix("/MAdata/db/"+db_name, name, reference)

        if no_time:
            def red(mat):
                return [ j for i in mat for j in i ]
            matrix = [[ red(red(matrix)) ]]

        #c = 1
        for row in matrix:
            #if c <= 0:
            #    print("break")
            #    break
            #c -= 1
            #count = 3
            for queries in row:
                #if count <= 0:
                #    print("break")
                #    break
                #count -= 1
                print(".", end="", flush=True)#print a line of dots
                #print("extracting " + str(len(queries)) + " samples (" + name + ")...")
                #setup the query pledges
                query_list = ContainerVector(NucSeq())
                # @todo temp bugfix
                del query_list[:]
                origins = []
                assert(len(query_list) == 0)
                for sequence, _, origin in queries:
                    query_list.append(NucSeq(sequence))
                    origins.append(origin)


                #print("setting up (" + name + ") ...")

                query_vec_pledge = Pledge(ContainerVector(NucSeq()))

                #fullfill the made promises
                query_vec_pledge.set(query_list)

                result_pledge = aligner.promise_me(query_vec_pledge, reference_pledge)

                #print("computing (" + name + ") ...")
                result_pledge.get()

                result = []
                for alignment in result_pledge.get():
                    #print(alignment.begin_on_ref, queries[int(alignment.stats.name)][-1])
                    result.append(
                        (
                            queries[int(alignment.stats.name)][1],
                            float('nan'),
                            alignment.begin_on_ref,
                            alignment.end_on_ref,
                            float('nan'),
                            aligner.elapsed_time,
                            alignment.mapping_quality,
                            name,
                            1 if alignment.secondary else 0
                        )
                    )
                #print("submitting results (" + name + ") ...")
                if len(result) > 0:
                    submitResults("/MAdata/db/"+db_name, result)
        print("")#print a newline
        aligner.final_checks()#do a final consistency check
    print("done working on " + db_name)

def test_all():
    test("test.db", human_genome)
    #test("default.db", human_genome)
    #test("short.db", human_genome)
    #test("long.db", human_genome)
    #test("shortIndels.db", human_genome)
    #test("longIndels.db", human_genome)
    #test("insertionOnly.db", human_genome)
    #test("deletionOnly.db", human_genome)
    #test("zoomLine.db", human_genome)
    #test("zoomSquare.db", human_genome)

    #test("illumina.db", human_genome)



def createSampleQueries(ref, db_name, size, indel_size, amount, reset = True, high_qual=False, in_to_del_ratio=0.5, smaller_box = False, only_first_row = False, validate_using_bwa=False, quiet=False):
    conn = sqlite3.connect("/MAdata/db/" + db_name)

    setUpDbTables(conn, reset)

    ref_seq = Pack()
    ref_seq.load(ref)
    reference_pledge = Pledge(Pack())
    reference_pledge.set(ref_seq)

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

    bwa = BWA_MEM(ref, 32, 1, db_name)
    genome = Pledge(NucSeq())
    if validate_using_bwa:
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

            #
            # extract random query sequences and modify them
            #
            queries = []
            query_list = ContainerVector(NucSeq())
            # @todo temp bugfix
            del query_list[:]
            for _ in range(amount):
                q_from, query, original_nuc_dist, modified_nuc_dist = get_query(ref_seq, size, mutation_amount, indel_amount, indel_size, in_to_del_ratio)
                queries.append( [q_from, query, original_nuc_dist, modified_nuc_dist] )

                #set up the validation pledges if we need them
                if validate_using_bwa:
                    query_list.append( NucSeq(query) )
            #
            # validate the query sequences using BWA MEM
            #
            if validate_using_bwa:
                query_vec_pledge = Pledge(ContainerVector(NucSeq()))
                query_vec_pledge.set(query_list)
                result_pledge = bwa.promise_me(query_vec_pledge, reference_pledge)

                problem_queries = []
                for alignment in result_pledge.get():
                    #check if there was an ambiguous alignment
                    if alignment.mapping_quality <= 0:
                        problem_queries.append( int(alignment.stats.name) )
                #get unique sorted (descending) elements
                problem_queries = list(reversed(sorted(set(problem_queries))))
                if len(problem_queries) > 0:
                    print("deleting", len(problem_queries), "queries")
                for index in problem_queries:
                    del queries[index]

            #
            # append the validated queries
            #
            for q_from, query, original_nuc_dist, modified_nuc_dist in queries:
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
                    if not quiet:
                        print("saving...")
                    insertQueries(conn, queries_list)
                    queries_list = []
                    if not quiet:
                        print("done saving")
            #for _ in range(amount) end

        #for indel_amount in range(max_indels) end
        if not quiet:
            print(mutation_amount, "/", max_x)
    #for mutation_amount in range(size) end
    if not quiet:
        print("saving...")
    if len(queries_list) > 0:
        insertQueries(conn, queries_list)
    if not quiet:
        print("done saving")

        print("nuc distrib (A, C, G, T, N) originally: ", nuc_distrib_count_orig)
        print("nuc distrib (A, C, G, T, N) changed: ", nuc_distrib_count_mod)
        print("nuc distrib (A, C, G, T, N) changed by: ", list(map(operator.mul, map(
            operator.sub, nuc_distrib_count_mod, nuc_distrib_count_orig), [1./float(sum(nuc_distrib_count_orig))]*5)))
        print("total amount: ", sum(nuc_distrib_count_orig))
#function
