import os
import subprocess
from MA import *
import sys
import time

class CommandLine(Module):

    def __init__(self):
        self.elapsed_time = 0

    def __get_sam(self, index_str, queries):
        f = open(self.in_filename, "w")
        for index, query in enumerate(queries):
            f.write(">" + str(index) + " sequence" + str(index) + "\n" )
            query_string = str(query)
            assert(len(query) > 0)
            assert(len(query_string) > 0)
            while len(query_string) > 60:
                line = query_string[:60]
                query_string = query_string[61:]
                f.write(line + "\n")
            f.write(query_string + "\n")


        f.close()
        #assemble the shell command
        cmd_str = self.create_command(self.in_filename)

        start_time = time.time()
        result = subprocess.run(cmd_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        self.elapsed_time = time.time() - start_time

        os.remove(self.in_filename)

        if result.returncode != 0:
            print("subprocess returned with ERROR:")
            print(result.stderr.decode('utf-8')[:-1])
            return None

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

        for _ in range(len(queries)):
            alignments.append(Alignment())

        for line in lines:
            if len(line) == 0:
                continue
            columns = line.split("\t")
            #print(columns[0], columns[2], columns[3], columns[4])
            try:
                #print(str(int(columns[3])) + " " + str(pack.start_of_sequence(columns[2])))
                start = pack.start_of_sequence(columns[2]) + int(columns[3])
                alignments[int(columns[0])] = Alignment(start)
            except:
                print("oh oh:" + line)
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
    def __init__(self, index_str, threads, db_name):
        self.bowtie2_home = "/usr/home/markus/workspace/bowtie2/bowtie2-2.3.3.1/"
        self.index_str = index_str + "bowtie2"
        self.threads = threads
        self.in_filename = ".temp" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.bowtie2_home + "bowtie2 -p " + str(self.threads)
        index_str = "-x " + self.index_str
        input_str = "-f -U " + in_filename
        return cmd_str + " " + index_str + " " + input_str

class Minimap2(CommandLine):
    def __init__(self, index_str, threads, db_name):
        self.minimap2_home = "/usr/home/markus/workspace/minimap2/"
        self.index_str = index_str + ".mmi"
        self.threads = threads
        self.in_filename = ".temp" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.minimap2_home + "minimap2 -t " + str(self.threads) + " -a "
        return cmd_str + " " + self.index_str + " " + in_filename

class Blasr(CommandLine):
    def __init__(self, index_str, threads, genome_str, db_name):
        self.blasr_home = "/usr/home/markus/workspace/blasr/build/bin/"
        self.index_str = index_str + "blasr"
        self.genome_str = genome_str
        self.threads = threads
        self.in_filename = ".temp" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.blasr_home + "blasr " + in_filename
        return cmd_str + " " + self.genome_str + " -m 1 --bestn 1 --nproc " + str(self.threads) + " --hitPolicy leftmost --sa " + self.index_str

class BWA_MEM(CommandLine):
    def __init__(self, index_str, threads, db_name):
        self.bwa_home = "/usr/home/markus/workspace/bwa/"
        self.index_str = index_str + "bwa"
        self.threads = threads
        self.in_filename = ".temp" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.bwa_home + "bwa mem -t " + str(self.threads)
        return cmd_str + " " + self.index_str + " " + in_filename

class BWA_SW(CommandLine):
    def __init__(self, index_str, threads, db_name):
        self.bwa_home = "/usr/home/markus/workspace/bwa/"
        self.index_str = index_str + "bwa"
        self.threads = threads
        self.in_filename = ".temp" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.bwa_home + "bwa bwasw -t " + str(self.threads)
        return cmd_str + " " + self.index_str + " " + in_filename

class MA(CommandLine):
    def __init__(self, index_str, threads, fast, db_name):
        self.ma_home = "/usr/home/markus/workspace/aligner/"
        self.index_str = index_str
        self.threads = threads
        self.fast = "accurate"
        if fast:
            self.fast = "fast"
        self.in_filename = ".temp" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.ma_home + "ma -a -t " + str(self.threads) + " -p " + self.fast
        return cmd_str + " -g " + self.index_str + " -i " + in_filename

human_genome = "/mnt/ssd0/genome/human"

def test(
            db_name,
            reference
        ):
    print("working on " + db_name + "...")
    ref_pack = Pack()
    ref_pack.load(reference)
    reference_pledge = Pledge(Pack())
    reference_pledge.set(ref_pack)

    num_threads = 1

    l = [
        ("BOWTIE 2", Bowtie2(reference, num_threads, db_name)),
        ("MINIMAP 2", Minimap2(reference, num_threads, db_name)),
        #("BLASR", Blasr(reference, num_threads, "/mnt/ssd0/genome/humanbwa")),
        ("BWA MEM", BWA_MEM(reference, num_threads, db_name)),
        ("BWA SW", BWA_SW(reference, num_threads, db_name)),
        #("MA Fast", MA(reference, num_threads, True, db_name)),
        #("MA Accurate", MA(reference, num_threads, False, db_name)),
    ]

    for name, aligner in l:
        print("evaluating " + name + "...")
        clearResults("/mnt/ssd1/"+db_name, reference, name) # CAREFUL WITH THE CLEARING
        result = []

        matrix = getQueriesAsASDMatrix("/mnt/ssd1/"+db_name, name, reference)
        for row in matrix:
            for queries in row:
                print(".", end="", flush=True)#print a line of dots
                #print("extracting " + str(len(queries)) + " samples (" + name + ")...")
                #setup the query pledges
                query_list = ContainerVector(NucSeq())
                # @todo temp bugfix
                del query_list[:]
                origins = []
                assert(len(query_list) == 0)
                for sequence, sample_id, origin in queries:
                    query_list.append(NucSeq(sequence))
                    origins.append(origin)


                #print("setting up (" + name + ") ...")

                query_vec_pledge = Pledge(ContainerVector(NucSeq()))

                #fullfill the made promises
                query_vec_pledge.set(query_list)

                result_pledge = aligner.promise_me(query_vec_pledge, reference_pledge)

                #print("computing (" + name + ") ...")
                result_pledge.get()

                for index in range(len(queries)):
                    alignment = result_pledge.get()[index]
                    #print(index, "->", alignment.begin_on_ref(), origins[index])
                    #total_time = result_pledge.exec_time / 1
                    sample_id = queries[index][1]
                    result.append(
                        (
                            sample_id,
                            float('nan'),
                            alignment.begin_on_ref,
                            alignment.end_on_ref,
                            float('nan'),
                            aligner.elapsed_time,
                            alignment.mapping_quality,
                            name
                        )
                    )
        #print("submitting results (" + name + ") ...")
        if len(result) > 0:
            submitResults("/mnt/ssd1/"+db_name, result)
    print("", end="")#print a newline
    print("done working on " + db_name)

def test_all():
    test("test.db", human_genome)
    test("default.db", human_genome)
    test("short.db", human_genome)
    test("long.db", human_genome)
    test("shortIndels.db", human_genome)
    test("longIndels.db", human_genome)

    #test("insertionOnly.db", human_genome)
    #test("deletionOnly.db", human_genome)
    #test("zoomLine.db", human_genome)
    #test("zoomSquare.db", human_genome)

    #test("illumina.db", human_genome)