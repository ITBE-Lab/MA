from .aligner import *

class SeedPrinter(Module):
    """prints out a list of Seeds.
Execution:
   Expects seedsList, query, ref
       seedsList: the list of Seeds
       query: the query as NucSeq
       ref: the reference as Pack
   returns: nothing
    """
    def __init__(self):
        pass

    #override
    def get_input_type(self):
        return [ContainerType.seedsList, ContainerType.nucSeq, ContainerType.packedNucSeq]

    #override
    def get_output_type(self):
        return ContainerType.nothing

    #override
    def execute(self, input):
        seedList = input[0]
        query = str(input[1])
        ref_pack = input[2]

        for seed in seedList:
            printstr = "(" + str(seed.start()) + "," + str(seed.end()) + "): "
            if seed.start() > 0:
                printstr += query[seed.start()-1] + "|"
            printstr += query[seed.start():seed.end()-1]
            if seed.end() < len(query)-1:
                printstr += "|" + query[seed.end()+1]
            printstr += " (" + str(seed.start_ref()) + "," + str(seed.end_ref()) + "): "
            if seed.start_ref() > 0:
                printstr += str(ref_pack.extract_from_to(seed.start_ref()-1, seed.start_ref()))
                printstr += "|"
            printstr += str(ref_pack.extract_from_to(seed.start_ref(), seed.end_ref()-1))
            if seed.end_ref() < ref_pack.unpacked_size()-1:
                printstr += "|"
                printstr += str(ref_pack.extract_from_to(seed.end_ref()-1, seed.end_ref()))
            print(printstr)

        return None
