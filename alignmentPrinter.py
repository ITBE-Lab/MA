from LAuS import *


class AlignmentPrinter(Module):
    """prints a Alignment.
Execution:
   Expects align, query, ref
       align: the Alignment
       query: the query as NucSeq
       ref: the reference as Pack
   returns: nothing
    """

    def __init__(self, nuc_per_line = 100):
        #Module.__init__()
        self.nuc_per_line = nuc_per_line

    #overrride
    def get_input_type(self):
        return [ContainerType.alignment, ContainerType.nucSeq, ContainerType.packedNucSeq]

    #overrride
    def get_output_type(self):
        return [ContainerType.nothing]

    #overrride
    def execute(self, input):
        align = input[0]
        query = input[1]
        ref_pack = input[2]

        ref = ref_pack.extract_from_to(align.begin_on_ref(), len(align) + align.begin_on_ref())

        lines = [
            "some desc for aligment here...:"
        ]
        counter = 0
        ind_query = 0
        ind_ref = 0

        while counter < len(align):
            #append three more lines if the current lines are full
            if counter % self.nuc_per_line == 0:
                desc = str(counter) + "-" + str(counter + self.nuc_per_line)
                lines.extend(["", desc, "", "", ""])

            #check for match or missmatch
            print ind_ref
            print ind_query
            if align[counter] is MatchType.match:
                lines[-3] += ref[ind_ref]
                lines[-2] += '|'
                lines[-1] += query[ind_query]
                ind_ref += 1
                ind_query += 1
            elif align[counter] is MatchType.missmatch:
                lines[-3] += ref[ind_ref]
                lines[-2] += ' '
                lines[-1] += query[ind_query]
                ind_ref += 1
                ind_query += 1
            elif align[counter] is MatchType.insertion:
                lines[-3] += '-'
                lines[-2] += ' '
                lines[-1] += query[ind_query]
                ind_query += 1
            elif align[counter] is MatchType.deletion:
                lines[-3] += ref[ind_ref]
                lines[-2] += ' '
                lines[-1] += '-'
                ind_ref += 1
            counter += 1

        for line in lines:
            print line

        return None
