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

        ref = ref_pack.extract_from_to(align.begin_on_ref(), align.end_on_ref())

        lines = [
            "some desc for aligment here...:"
        ]
        counter = 0
        ind_query = 0
        ind_ref = 0

        desc = str(0) + "-" + str(0 + self.nuc_per_line)
        lines.extend(["", desc, "", "", ""])

        atLeastOneMistake = False

        while counter < len(align):
            #append three more lines if the current lines are full
            if counter % self.nuc_per_line == 0:
                if counter > 0:
                    lines[-3] += "\treference"
                    lines[-1] += "\tquery"
                    desc = str(counter) + "-" + str(counter + self.nuc_per_line)
                    lines.extend(["", desc, "", "", ""])

            #perform double check for messup:
            if ind_ref >= len(ref) and (align[counter] == MatchType.match or align[counter] == MatchType.deletion or align[counter] == MatchType.missmatch):
                print "This should not happen... (ref)"
                print ind_ref
                print len(ref)
                break
            if ind_query >= len(query) and (align[counter] == MatchType.match or align[counter] == MatchType.insertion or align[counter] == MatchType.missmatch):
                print "This should not happen... (query)"
                print ind_query
                print len(query)
                break

            #check for match or missmatch
            #print str(ind_ref) + " of " + str(align.end_on_ref() - align.begin_on_ref())
            if align[counter] is MatchType.match:
                lines[-3] += ref[ind_ref]
                if ref[ind_ref] != query[ind_query]:
                    lines[-2] += 'X'
                    atLeastOneMistake = True
                else:
                    lines[-2] += '|'
                lines[-1] += query[ind_query]
                ind_ref += 1
                ind_query += 1
            if align[counter] is MatchType.seed:
                lines[-3] += ref[ind_ref]
                lines[-2] += 'I'
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

        while len(lines[-1]) < 100:
            lines[-3] += '-'
            lines[-2] += ' '
            lines[-1] += '-'
        lines[-3] += "\treference"
        lines[-1] += "\tquery"

        for line in lines:
            print line
        if atLeastOneMistake:
            print "WARNING: the alignment contains errors!"

        return None
