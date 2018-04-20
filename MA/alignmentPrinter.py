##
# @package MA.alignmentPrinter
# @brief Implements @ref MA.alignmentPrinter.AlignmentPrinter "AlignmentPrinter".
# @file alignmentPrinter.py
# @brief Implements @ref MA.alignmentPrinter.AlignmentPrinter "AlignmentPrinter".
# @author Markus Schmidt

from .aligner import *


##
# @brief Print the given Alignment.
# @details 
# Prints the alignment to the console.
# @ingroup module
#
class AlignmentPrinter(Module):

    def __init__(self, nuc_per_line = 80, check_for_errors_only=False):
        self.nuc_per_line = nuc_per_line
        self.check_for_errors_only = check_for_errors_only

    ##
    # @brief returns the @ref ContainerType "container types" alignment, nucSeq, packedNucSeq.
    # @details
    # Reimplemented from MA.aligner.Module.get_input_type.
    def get_input_type(self):
        return [Alignment(), NucSeq(), Pack()]

    ##
    # @brief returns the @ref ContainerType "container type" nothing.
    # @details
    # Reimplemented from MA.aligner.Module.get_output_type.
    def get_output_type(self):
        return Nil()

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, *input):
        align = input[0]
        query = input[1][align.begin_on_query:align.end_on_query]
        ref = input[2].extract_from_to(align.begin_on_ref, align.end_on_ref)

        lines = [
            "score: " + str(align.get_score()) + " cigar length: " + str(len(align)) + 
            " soc index:" + str(align.stats.index_of_strip) + 
            " map qual:" + str(align.mapping_quality) + 
            " query name:" + str(align.stats.name),
            "reference: " + str(align.begin_on_ref) + " - " + str(align.end_on_ref),
            "query: " + str(align.begin_on_query) + " - " + str(align.end_on_query)
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
            if ind_ref >= len(ref) and (align[counter] == MatchType.match or align[counter] == MatchType.seed or align[counter] == MatchType.deletion or align[counter] == MatchType.missmatch):
                print("This should not happen... (ref)")
                print(align[counter])
                print(ind_ref)
                print(len(ref))
                #@todo it would be nice if the cpp code would simply support slicing the alignment
                s = []
                for index in range(counter, len(align)):
                    s.append(align[counter])
                print("remaining cigar:", s)
                atLeastOneMistake = True
                break
            if ind_query >= len(query) and (align[counter] == MatchType.match or align[counter] == MatchType.seed or align[counter] == MatchType.insertion or align[counter] == MatchType.missmatch):
                print("This should not happen... (query)")
                print(align[counter])
                print(ind_query)
                print(len(query))
                s = []
                for index in range(counter, len(align)):
                    s.append(align[counter])
                print("remaining cigar:", s)
                atLeastOneMistake = True
                break

            #check for match or missmatch
            #print str(ind_ref) + " of " + str(align.end_on_ref() - align.begin_on_ref())
            if align[counter] is MatchType.match:
                if ref[ind_ref] != query[ind_query]:
                    lines[-2] += 'x'
                    atLeastOneMistake = True
                else:
                    lines[-2] += '|'
                lines[-3] += ref[ind_ref]
                lines[-1] += query[ind_query]
                ind_ref += 1
                ind_query += 1
            if align[counter] is MatchType.seed:
                if ref[ind_ref] != query[ind_query]:
                    lines[-2] += 'X'
                    atLeastOneMistake = True
                else:
                    lines[-2] += 'I'
                lines[-3] += ref[ind_ref]
                lines[-1] += query[ind_query]
                ind_ref += 1
                ind_query += 1
            elif align[counter] is MatchType.missmatch:
                if ref[ind_ref] == query[ind_query]:
                    lines[-2] += ':'
                    atLeastOneMistake = True
                else:
                    lines[-2] += ' '
                lines[-3] += ref[ind_ref]
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

        while len(lines[-1]) < self.nuc_per_line:
            lines[-3] += '-'
            lines[-2] += ' '
            lines[-1] += '-'
        lines[-3] += "\treference"
        lines[-1] += "\tquery"

        if not self.check_for_errors_only and not atLeastOneMistake:
            for line in lines:
                print(line)
        if atLeastOneMistake:
            for line in lines:
                print(line)
            print("WARNING: the alignment contains errors!")

        return None
