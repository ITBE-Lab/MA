from .aligner import *
from .alignmentPrinter import *

def test_aligner():
    reference_pack = Pack()
    reference_pack.load("/MAdata/genome/eColi")
    reference_index = FMIndex()
    reference_index.load("/MAdata/genome/eColi")

    query = NucSeq("CGCGCGCGCGCGCGCG")

    # initialize modules
    seeding_module = BinarySeeding()
    soc_module = StripOfConsideration()
    harm_module = Harmonization()
    nw_module = NeedlemanWunsch()
    printer_module = AlignmentPrinter()

    #execute modules
    seeds = seeding_module.execute(reference_index, query)
