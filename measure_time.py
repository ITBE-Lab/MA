import os
import subprocess
from MA import *
import sys
import time

class CommandLine(Module):

    def __init__(self):
        self.elapsed_time = 0
        self.skip_count = 0
        self.warn_once = True

    def final_checks(self):
        if not self.do_checks():
            return
        if self.skip_count > 0:
            print("Aligner skipped", self.skip_count, "queries")

    def __get_sam(self, index_str, queries):
        f = open(self.in_filename, "w")
        for query in queries:
            check = ""
            f.write(">" + query.name + " sequence" + query.name + "\n" )
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

    ##
    # @origin https://stackoverflow.com/questions/10321978/integer-to-bitfield-as-a-list
    def bitfield(self, num):
        return list(reversed([True if digit=='1' else False for digit in bin(num)[2:]]))

    def secondary(self, string):
        bits = self.bitfield(int(string))
        if len(bits) <= 8:
            return False
        return bits[8]

    def __align(self, index_str, queries, pack):
        sam = self.__get_sam(index_str, queries)
        #print(sam)

        lines = sam.split("\n")
        while len(lines) > 0 and ( len(lines[0]) == 0 or lines[0][0] is '@' ):
            lines = lines[1:]

        #transform sam file into list data structure
        alignments = []

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
                alignment.secondary = True if self.secondary(columns[1]) else False
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
    num_results = "5"

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

        matrix = getQueriesAsASDMatrix("/MAdata/db/"+db_name)

        if no_time:
            def red(mat):
                return [ j for i in mat for j in i ]
            matrix = [[ red(red(matrix)) ]]

        #c = 1
        for mut_amount, row in enumerate(matrix):
            #if c <= 0:
            #    print("break")
            #    break
            #c -= 1
            #count = 3
            for indel_amount, queries in enumerate(row):
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
                assert(len(query_list) == 0)
                for sequence, sample_id in queries:
                    query_list.append(NucSeq(sequence))
                    query_list[-1].name = str(sample_id)


                #print("setting up (" + name + ") ...")

                query_vec_pledge = Pledge(ContainerVector(NucSeq()))

                #fullfill the made promises
                query_vec_pledge.set(query_list)

                result_pledge = aligner.promise_me(query_vec_pledge, reference_pledge)

                #print("computing (" + name + ") ...")
                result_pledge.get()

                result = []
                times = []
                for alignment in result_pledge.get():
                    #print(alignment.begin_on_ref, queries[int(alignment.stats.name)][-1])
                    result.append(
                        (
                            int(alignment.stats.name),
                            alignment.begin_on_ref,
                            alignment.end_on_ref,
                            alignment.mapping_quality,
                            name,
                            1 if alignment.secondary else 0
                        )
                    )
                    if not no_time:
                        times.append(
                            (
                                mut_amount,
                                indel_amount,
                                aligner.elapsed_time,
                                name
                            )
                        )
                #print("submitting results (" + name + ") ...")
                if len(result) > 0:
                    submitResults("/MAdata/db/"+db_name, result)
                if len(times) > 0 and not no_time:
                    submitRuntimes("/MAdata/db/"+db_name, times)
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

