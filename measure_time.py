import os
import subprocess
from MA import *
from MA.createSampleSeq import *
import sys
import time

class CommandLine(Module):

    def __init__(self):
        self.elapsed_time = 0.0
        self.skip_count = 0
        self.warn_once = True
        self.processor = 0

    def final_checks(self):
        if not self.do_checks():
            return
        if self.skip_count > 0:
            print("Aligner skipped", self.skip_count, "queries")

    def __get_sam(self, queries):
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

        bitmask = [1]
        for i in range(0, self.processor):
            bitmask.append( 0 )
        hex_num = ""
        while len(bitmask) > 0:
            num = 0
            for index, val in enumerate(reversed(bitmask[-4:])):
                num += val*2**index
            hex_num = str(num) + hex_num
            bitmask = bitmask[:-4]

        taskset = "taskset " + hex_num + " "

        # default command
        if True:
            #assemble the shell command
            cmd_str = taskset + self.create_command(self.in_filename)

            start_time = time.time()
            result = subprocess.run(cmd_str, stdout=subprocess.PIPE, shell=True)
            self.elapsed_time = time.time() - start_time

            if result.returncode != 0:
                print("call command:", cmd_str)
                print("subprocess returned with ERROR:")
                print("Error:")
                print(result.stderr.decode('utf-8'))
                exit()

            os.remove(self.in_filename)

            sam_file = result.stdout.decode('utf-8')

            return sam_file

        else:
            #assemble the shell command
            #cmd_str = taskset + self.create_command(self.in_filename)
            cmd_str = taskset + self.create_command(self.in_filename) + " -o .tempSamOut"
            #cmd_str = self.create_command(self.in_filename)

            start_time = time.time()
            #result = subprocess.run(cmd_str, stdout=subprocess.PIPE, shell=True)
            result = subprocess.run(cmd_str, shell=True)
            self.elapsed_time = time.time() - start_time

            if result.returncode != 0:
                print("call command:", cmd_str)
                print("subprocess returned with ERROR:")
                print("Error:")
                print(result.stderr.decode('utf-8'))
                exit()

            os.remove(self.in_filename)

            #sam_file = result.stdout.decode('utf-8')
            sam_file = ""
            with open('.tempSamOut', 'r') as f:
                sam_file=f.read()

            return sam_file

    ##
    # @origin https://stackoverflow.com/questions/10321978/integer-to-bitfield-as-a-list
    def bitfield(self, num):
        return list(reversed([True if digit=='1' else False for digit in bin(num)[2:]]))

    def check_flag(self, string, flag_bit):
        bits = self.bitfield(int(string))
        if len(bits) <= flag_bit:
            return False
        return bits[flag_bit]

    def secondary(self, string):
        return self.check_flag(string, 8)

    def supplementary(self, string):
        return self.check_flag(string, 11)

    def __align(self, queries, pack):
        sam = self.__get_sam(queries)
        #print(sam)

        lines = sam.split("\n")
        if self.output_type() == "SAM":
            while len(lines) > 0 and ( len(lines[0]) == 0 or lines[0][0] is '@' ):
                lines = lines[1:]

        #transform sam file into list data structure
        alignments = []

        # @todo add me as a separate item in the database
        determine_secondary_by_order = False
        secondary_ordered = {}
        secondary_list = []

        for line in lines:
            # ignore empty lines
            if len(line) == 0:
                continue
            # ignore dummy samples
            if line[:6] == "dummy_":
                continue
            if self.output_type() == "SAM":
                columns = line.split("\t")

                ## Verbose output
                # for column in columns:
                #     if len(column) <= 20:
                #         print(column, end=" ")
                #     else:
                #         print(column[:20], "...", sep="", end=" ")
                # print()

                #print(line)
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
                    qLen = 0
                    rLen = 0
                    for amount, char in read_cigar(columns[5]):
                        if char in ['M', 'I', '=', 'X', 'S']:
                            qLen += amount
                        if char in ['M', 'D', '=', 'X', 'N']:
                            rLen += amount
                        if char in ['M', 'X', '=', 'D']:
                            align_length += amount
                        # sanity check...
                        elif not char in ['M', 'I', 'D', '=', 'X', 'S', 'N', 'H', 'P']:
                            print("Error: got wierd cigar symbol", char, "in cigar", columns[5])
                            exit()
                    start = pack.start_of_sequence(columns[2]) + int(columns[3])
                    # 0 is not correct actually but we do not save the query pos anyway...
                    alignment = Alignment(start, 0, start + rLen, qLen)
                    alignment.mapping_quality = int(columns[4])/255
                    alignment.stats.name = str(int(columns[0]))
                    # set the secondary flag for secondary alignments
                    alignment.secondary = True if self.secondary(columns[1]) else False

                    if self.supplementary(columns[1]):
                        alignment.secondary = True
                        secondary_list.append(alignment.stats.name)

                    #overwrite secondary by order of the output
                    if determine_secondary_by_order:
                        alignment.secondary = alignment.stats.name in secondary_ordered
                        secondary_ordered[alignment.stats.name] = True

                    alignments.append(alignment)
                except Exception as e:
                    print("Error:", e)
                    print(line)
                    traceback.print_exc()
                    pass
            elif self.output_type() == "BLASR":
                columns = line.split(" ")
                #print(columns)
                start = int(columns[6])
                length = int(columns[7]) - int(columns[6])
                #
                # for blasr it seems like we have to invert to the correct position 
                # in case of output on the reverse complement
                # also the inversion is chromosome specific i.e. in my world:
                # s1 - s2 - s3 - ~s3 - ~s2 - ~s1
                # whereas in blasr's world:
                # s1 - ~s1 - s2 - ~s2 - s3 - ~s3
                #
                if columns[3] == '1':
                    start = pack.length_of_sequence(columns[1]) - (start + length)
                start += pack.start_of_sequence(columns[1])
                # 0 is not correct actually but we do not save the query pos anyway...
                alignment = Alignment(start, 0, start + length, int(columns[11]))
                alignment.mapping_quality = int(columns[4]) # score instead of map. qual
                alignment.stats.name = columns[0]
                alignment.secondary = False # does not give primary secondary information
                
                #overwrite secondary by order of the output
                if determine_secondary_by_order:
                    alignment.secondary = alignment.stats.name in secondary_ordered
                    secondary_ordered[alignment.stats.name] = True

                alignments.append(alignment)

            elif self.output_type() == "FINDER":
                columns = line.split("\t")
                #print(columns)
                r_start = int(columns[4])
                r_length = int(columns[5])
                r_start += pack.start_of_sequence(columns[3])
                q_start = int(columns[1])
                q_length = int(columns[2])

                r_pos = r_start + int(r_length / 2)
                q_pos = q_start + int(q_length / 2)

                alignment = Alignment(r_pos, q_pos, r_pos, q_pos)
                alignment.mapping_quality = int(columns[8]) # acc seed length instead of map. qual
                alignment.stats.name = columns[0]
                alignment.secondary = not columns[6] == "true"
                #print(columns[6], "?=", columns[6] == "true")
                
                #overwrite secondary by order of the output
                if determine_secondary_by_order:
                    alignment.secondary = alignment.stats.name in secondary_ordered
                    secondary_ordered[alignment.stats.name] = True

                alignments.append(alignment)


        #transform list into alignment data structure
        ret = ContainerVector(Alignment())
        del ret[:]

        if len(secondary_list) > 0 and len(secondary_list) < 5:
            print("WARNING: Aligner outputted chimeric alignment for:", secondary_list)
        elif len(secondary_list) > 0:
            print("WARNING: Aligner outputted chimeric alignment for:", secondary_list[:5], "...")

        for alignment in alignments:
            ret.append(alignment)

        return ret

    def get_input_type(self):
        return [ContainerVector(NucSeq()), Pack()]

    def get_output_type(self):
        return ContainerVector(Alignment())

    def output_type(self):
        return "SAM"

    def execute(self, *input):
        #print(input)
        queries = input[0]
        pack = input[1]

        try:
            return self.__align(queries, pack)
        except:
            print("aligner crashed...")
            # return an empty alignment vector...
            ret = ContainerVector(Alignment())
            del ret[:]
            return ret

class Bowtie2(CommandLine):
    def __init__(self, index_str, num_results, db_name):
        super().__init__()
        self.bowtie2_home = "/usr/home/markus/workspace/bowtie2/bowtie2-2.3.3.1/"
        self.index_str = index_str + "bowtie2"
        self.num_results = num_results
        self.in_filename = ".tempBowtie2" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.bowtie2_home + "bowtie2 "
        index_str = "-x " + self.index_str
        input_str = "-f --very-fast -U " + in_filename
        return cmd_str + " " + index_str + " " + input_str + " -k " + self.num_results

    def do_checks(self):
        return False

class Minimap2(CommandLine):
    def __init__(self, index_str, num_results, db_name, z_drop=None, presetting=None):
        super().__init__()
        self.minimap2_home = "/usr/home/markus/workspace/minimap2/"
        self.index_str = index_str + ".mmi"
        self.num_results = num_results
        self.in_filename = ".tempMinimap2" + db_name + ".fasta"
        self.presetting = presetting
        if not z_drop is None:
            self.z_drop = " -z " + str(z_drop)
        else:
            self.z_drop = ""

    def create_command(self, in_filename):
        cmd_str = self.minimap2_home + "minimap2 -c -a "
        if int(self.num_results) > 1:
            in_filename += " --secondary=yes -N " + self.num_results
        if not self.presetting is None:
            cmd_str += " -x " + self.presetting
        return cmd_str + " " + self.index_str + " " + in_filename + self.z_drop

    def do_checks(self):
        return False

class Blasr(CommandLine):
    def __init__(self, index_str, num_results, genome_str, db_name):
        super().__init__()
        self.blasr_home = "/usr/home/markus/workspace/blasr/build/bin/"
        self.index_str = index_str + "blasr"
        self.genome_str = genome_str
        self.num_results = num_results
        self.in_filename = ".tempBlasr" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.blasr_home + "blasr " + in_filename
        # --nproc 32
        return cmd_str + " " + self.genome_str + " -m 1 --bestn " + self.num_results + " --sa " + self.index_str

    def do_checks(self):
        return False

    def output_type(self):
        return "BLASR"

class Ngmlr(CommandLine):
    def __init__(self, genome_file, db_name):
        super().__init__()
        self.ngmlr_home = "/usr/home/markus/workspace/ngmlr/bin/ngmlr-0.2.8/"
        self.genome_file = genome_file
        self.in_filename = ".tempNgmlr" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.ngmlr_home + "ngmlr -r " + self.genome_file + " -q " + in_filename
        print(cmd_str)
        return cmd_str

    def do_checks(self):
        return False

class BWA_MEM(CommandLine):
    def __init__(self, index_str, num_results, db_name, z_drop=None, presetting=None):
        super().__init__()
        self.bwa_home = "/usr/home/markus/workspace/bwa/"
        self.index_str = index_str + "bwa"
        self.num_results = num_results
        self.in_filename = ".tempBwaMem" + db_name + ".fasta"
        self.presetting = presetting
        if not z_drop is None:
            self.z_drop = " -d " + str(z_drop)
        else:
            self.z_drop = ""

    def create_command(self, in_filename):
        cmd_str = self.bwa_home + "bwa mem "
        if int(self.num_results) > 1:
            cmd_str += " -a"
        if not self.presetting is None:
            cmd_str += " -x " + self.presetting
        return cmd_str + " " + self.index_str + " " + in_filename + self.z_drop

    def do_checks(self):
        return False

class BWA_SW(CommandLine):
    def __init__(self, index_str, num_results, db_name, s=None):
        super().__init__()
        self.bwa_home = "/usr/home/markus/workspace/bwa/"
        self.index_str = index_str + "bwa"
        self.num_results = num_results
        self.in_filename = ".tempBwaSw" + db_name + ".fasta"
        self.s = s

    def create_command(self, in_filename):
        cmd_str = self.bwa_home + "bwa bwasw "
        if not self.s is None:
            return cmd_str + " " + self.index_str + " " + in_filename + " -s " + self.s + " -T 15 -c 10"
        return cmd_str + " " + self.index_str + " " + in_filename

    def do_checks(self):
        return False
        

class GEM(CommandLine):
    def __init__(self, index_str, num_results, genome_str, db_name):
        super().__init__()
        self.gem_home = "/usr/home/markus/workspace/gemtools/GEMTools/bin/"
        self.index_str = index_str + ".gem"
        self.num_results = num_results
        self.genome_str = genome_str
        self.in_filename = ".tempGEM" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.gem_home + "gt.map2sam"
        return cmd_str + " -I " + self.index_str + " -r " + self.genome_str + " -i " + in_filename

    def do_checks(self):
        return False


class G_MAP(CommandLine):
    def __init__(self, index_str, num_results, genome_str, db_name):
        super().__init__()
        self.g_home = "/usr/home/markus/workspace/graphmap/bin/Linux-x64/"
        self.genome_str = genome_str
        self.num_results = num_results
        self.index_str = index_str
        self.in_filename = ".tempBlasr" + db_name + ".fasta"

    def create_command(self, in_filename):
        cmd_str = self.g_home + "graphmap align -v 0 -r " + self.genome_str
        minZ = ""
        if int(self.num_results) > 1:
            minZ = " -Z"
        return cmd_str + " -d " + in_filename + minZ

    def do_checks(self):
        return False

class MA(CommandLine):
    def __init__(self, index_str, num_results, fast, db_name, finder_mode=False, other_dp_scores=False, soc_width=None):
        super().__init__()
        self.ma_home = "/usr/home/markus/workspace/aligner/"
        self.index_str = index_str
        self.num_results = num_results
        self.fast = "acc"
        if fast:
            self.fast = "fast"
        self.in_filename = ".tempMA" + self.fast + db_name + ".fasta"
        self.finder_mode = finder_mode
        self.other_dp_scores = other_dp_scores
        self.soc_width = soc_width

    def create_command(self, in_filename):
        cmd_str = self.ma_home + "ma -t 1 -m " + self.fast
        if self.finder_mode:
            cmd_str += " -d"
        if self.other_dp_scores:
            cmd_str += " --Match 1 --MisMatch 1 --Gap 1 --Extend 1"
        if not self.soc_width is None:
            cmd_str += " --SoCWidth " + self.soc_width
        return cmd_str + " -x " + self.index_str + " -i " + in_filename + " -n " + self.num_results

    def do_checks(self):
        return True

    def output_type(self):
        if self.finder_mode:
            return "FINDER"
        else:
            return "SAM"

human_genome = "/MAdata/genome/human"

def split_reads(reference, ref_seq):
    ref_pack = Pack()
    ref_pack.load(reference)
    reference_pledge = Pledge(Pack())
    reference_pledge.set(ref_pack)

    num_results = "2"

    warned_for_n = False

    l = [
        ("MA Fast", MA(reference, num_results, True, "")),
        ("BWA MEM", BWA_MEM(reference, num_results, "")),
        ("MINIMAP2", Minimap2(reference, num_results, "")),
        ("NGMLR", Ngmlr(ref_seq, "")),
    ]

    result_list = []
    query_list = ContainerVector(NucSeq())
    reads = createPacBioReadsSimLord(ref_seq, reference)

    for sample_id, sample in enumerate(reads):
        origin_a, len_a, origin_b, len_b, sequence = sample
        query_list.append(NucSeq(sequence))
        query_list[-1].name = str(sample_id)

    query_vec_pledge = Pledge(ContainerVector(NucSeq()))
    query_vec_pledge.set(query_list)
    for name, aligner in l:
        result_pledge = aligner.promise_me(query_vec_pledge, reference_pledge)

        tries = 0
        found_list = []
        for _ in range(len(reads)):
            found_list.append([False, False])

        aligned_sample_ids = set()

        for alignment in result_pledge.get():
            sample_id = int(alignment.stats.name)
            aligned_sample_ids.add(sample_id)
            # alignment.begin_on_ref, alignment.end_on_ref

            origin_a, len_a, origin_b, len_b, sequence = reads[sample_id]

            if near(alignment.begin_on_ref, origin_a, alignment.end_on_ref, origin_a + len_a):
                found_list[sample_id][0] = True
            if near(alignment.begin_on_ref, origin_b, alignment.end_on_ref, origin_b + len_b):
                found_list[sample_id][1] = True
            tries += 1

        both = 0
        one = 0
        none = 0
        for sample_id, found in enumerate(found_list):
            if found[0] and found[1]:
                both += 1
            elif found[0] or found[1]:
                one += 1
            else:
                none += 1
                print("WARNING: aligner", name, "found none for", sample_id, reads[sample_id])

        result_list.append( (name, both, one, none, len(reads), tries, aligner.elapsed_time) )

    for name, both, one, none, num, tries, aligner.elapsed_time in result_list:
        print(name, "found \tboth:", both, "one:", one, "none:", none,
            "of", num, "split reads using", tries, "alignments and",
            aligner.elapsed_time, "seconds")

def genome_dup_reads():
    reference = "/MAdata/genome/GRCh38.p12"
    ref_seq = "/MAdata/chrom/human/GCA_000001405.27_GRCh38.p12_genomic.fna"

    ref_pack = Pack()
    ref_pack.load(reference)
    reference_pledge = Pledge(Pack())
    reference_pledge.set(ref_pack)

    num_results = "2"

    warned_for_n = False

    l = [
        ("MA Fast", MA(reference, num_results, True, "")),
        ("BWA MEM", BWA_MEM(reference, num_results, "")),
        ("MINIMAP2", Minimap2(reference, num_results, "")),
        ("NGMLR", Ngmlr(ref_seq, "")),
    ]

    result_list = []
    query_list = ContainerVector(NucSeq())
    reads = createPacBioReadsSimLordGenDup()

    for sample_id, sample in enumerate(reads):
        origin, length, sequence = sample
        query_list.append(NucSeq(sequence))
        query_list[-1].name = str(sample_id)

    query_vec_pledge = Pledge(ContainerVector(NucSeq()))
    query_vec_pledge.set(query_list)
    for name, aligner in l:
        result_pledge = aligner.promise_me(query_vec_pledge, reference_pledge)

        tries = 0
        found_list = []
        for _ in range(len(reads)):
            found_list.append(False)

        for alignment in result_pledge.get():
            sample_id = int(alignment.stats.name)
            # alignment.begin_on_ref, alignment.end_on_ref

            origin, length, sequence = reads[sample_id]

            if near(alignment.begin_on_ref, origin, alignment.end_on_ref, origin + length):
                found_list[sample_id] = True
            tries += 1

        one = 0
        none = 0
        for sample_id, found in enumerate(found_list):
            if found:
                one += 1
            else:
                none += 1
                #print("WARNING: aligner", name, "found none for", sample_id, reads[sample_id])

        result_list.append( (name, one, none, len(reads), tries, aligner.elapsed_time) )

    for name, one, none, num, tries, aligner.elapsed_time in result_list:
        print(name, "found \t:", one, "missed:", none,
            "of", num, "split reads using", tries, "alignments and",
            aligner.elapsed_time, "seconds")

def test(
            db_name,
            reference,
            only_overall_time=True,
            long_read_aligners=True,
            short_read_aligners=True,
            runtime_sample_multiplier=0,
            processor=None,
            specific_sample=None,
        ):
    print("working on " + db_name)
    ref_pack = Pack()
    ref_pack.load(reference)
    reference_pledge = Pledge(Pack())
    reference_pledge.set(ref_pack)

    num_results = "1"

    warned_for_n = False

    g_map_genome = "/MAdata/chrom/" + reference.split('/')[-1] + "/n_free.fasta"

    l = [
        ##("MA Accurate -w 10", MA(reference, num_results, False, db_name, soc_width="10")),
        ##("MA Accurate -w 100", MA(reference, num_results, False, db_name, soc_width="100")),
        ##("MA Accurate -w 300", MA(reference, num_results, False, db_name, soc_width="300")),
        ##("BWA SW s=30", BWA_SW(reference, num_results, db_name, s="300")),


        ("MA Fast", MA(reference, num_results, True, db_name)),
        #("MA Basic", MA(reference, num_results, True, db_name, finder_mode=True)),

        ("BWA MEM", BWA_MEM(reference, num_results, db_name)),
        ("MINIMAP2", Minimap2(reference, num_results, db_name)),
        #  #
        #  ("BWA MEM pacbio", BWA_MEM(reference, num_results, db_name, presetting="pacbio")),
        #  ("BWA MEM ont2d", BWA_MEM(reference, num_results, db_name, presetting="ont2d")),
        #  ("BWA MEM intractg", BWA_MEM(reference, num_results, db_name, presetting="intractg")),
        #  ("MINIMAP 2 map-pb", Minimap2(reference, num_results, db_name, presetting="map-pb")),
        #  ("MINIMAP 2 map-ont", Minimap2(reference, num_results, db_name, presetting="map-ont")),
        #  ("MINIMAP 2 asm10", Minimap2(reference, num_results, db_name, presetting="asm10")),
        #  #
        #  # ("BWA MEM 0 zDrop", BWA_MEM(reference, num_results, db_name, z_drop=0)),
        #  # ("MINIMAP 2 0 zDrop", Minimap2(reference, num_results, db_name, z_drop=0)),
        #  #

        #("MA Accurate", MA(reference, num_results, False, db_name)),

        #  ("BWA SW", BWA_SW(reference, num_results, db_name)),
    
        ("NGMLR", Ngmlr(g_map_genome, db_name)),
    ]

    ## if long_read_aligners:
    ##     l.extend([
    ##             ("GRAPH MAP", G_MAP(reference, num_results, g_map_genome, db_name)),
    ##             
    ##         ])
    ## 
    ## if short_read_aligners:
    ##     l.extend([
    ##             ("BOWTIE 2", Bowtie2(reference, num_results, db_name)),
    ##             ("BLASR", Blasr(reference, num_results, g_map_genome, db_name)),
    ##         ])

    for name, aligner in l:
        print("evaluating " + name)

        if not processor is None:
            aligner.processor = processor

        matrix = None

        if specific_sample is None:
            matrix = getQueriesAsASDMatrix("/MAdata/db/"+db_name)

            if only_overall_time:
                def red(mat):
                    return [ j for i in mat for j in i ]
                matrix = [[ red(red(matrix)) ]]
            clearApproach("/MAdata/db/"+db_name, name)
        else:
            matrix = [[ getSpecificQuery("/MAdata/db/"+db_name, specific_sample) ]]

        #c = 1
        total_time = 0.0
        total_queries = 0
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

                query_list_remove_load_time = ContainerVector(NucSeq())

                assert(len(query_list) == 0)
                assert(len(query_list_remove_load_time) == 0)

                for sequence, sample_id in queries:
                    # check for N's
                    if not warned_for_n:
                        for nuc in sequence:
                            if nuc not in ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']:
                                warned_for_n = True
                                print("Queries contain:", nuc)

                    for _ in range(0, runtime_sample_multiplier):
                        query_list.append(NucSeq(sequence))
                        query_list[-1].name = "dummy_" + str(sample_id)
                    query_list.append(NucSeq(sequence))
                    query_list[-1].name = str(sample_id)

                    if len(query_list_remove_load_time) == 0:
                        query_list_remove_load_time.append(NucSeq(sequence))
                        query_list_remove_load_time[-1].name = str(sample_id)
                assert(len(query_list_remove_load_time) <= 1)

                #print("setting up (" + name + ") ...")

                query_vec_pledge = Pledge(ContainerVector(NucSeq()))

                #fullfill the made promises
                query_vec_pledge.set(query_list)

                #print("num queries:", len(query_list))
                

                result_pledge = aligner.promise_me(query_vec_pledge, reference_pledge)

                #print("computing (" + name + ") ...")

                result = []
                times = []
                for alignment in result_pledge.get():
                    #print(alignment.begin_on_ref, queries[int(alignment.stats.name)][-1])
                    # @note the query position in the alignment is not set correctly

                    result.append(
                        (
                            int(alignment.stats.name),
                            alignment.begin_on_ref,
                            alignment.end_on_ref,
                            alignment.mapping_quality,
                            name,
                            1 if alignment.secondary else 0,
                            0 # dummy value
                        )
                    )
                if not only_overall_time:
                    times.append(
                        (
                            mut_amount,
                            indel_amount,
                            aligner.elapsed_time / len(query_list),
                            name
                        )
                    )

                #print("submitting results (" + name + ") ...")
                if len(result) > 0:
                    submitResults("/MAdata/db/"+db_name, result)
                if len(times) > 0 and not only_overall_time:
                    submitRuntimes("/MAdata/db/"+db_name, times)

                #update the overall total runtime
                if True: # substract the index load time
                    total_queries += len(query_list) - 1
                    total_time_this = aligner.elapsed_time
                    aligner.elapsed_time = 0
                    # let the aligner run on one single query in order to remove the index load time
                    query_vec_pledge.set(query_list_remove_load_time)
                    result_pledge = aligner.promise_me(query_vec_pledge, reference_pledge)
                    result_pledge.get()
                    assert(total_time_this > aligner.elapsed_time)
                    total_time += (total_time_this - aligner.elapsed_time) 
                else: # do not subtract index load times.
                    total_queries += len(query_list)
                    total_time += aligner.elapsed_time

                #just overwrite the value multiple times
                if total_queries > 0:
                    putTotalRuntime("/MAdata/db/" + db_name, name, total_time / total_queries)

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

