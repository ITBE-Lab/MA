from .aligner import *
from .alignmentPrinter import *

def test_aligner():
    configureFast()

    reference_pack = Pack()
    reference_pack.load("C:\\Users\\Markus\\Desktop\\MAdata\\eColi")
    reference_index = FMIndex()
    reference_index.load("C:\\Users\\Markus\\Desktop\\MAdata\\eColi")

    query = reference_pack.extract_from_to(10000, 10200)

    # initialize modules
    seeding_module = BinarySeeding()
    soc_module = StripOfConsideration()
    harm_module = Harmonization()
    nw_module = NeedlemanWunsch()
    printer_module = AlignmentPrinter()

    #execute modules
    seeds = seeding_module.execute(reference_index, query)
    print("seed set (q,r,l):")
    for seed in seeds.extract_seeds(reference_index, 100, 0, True):
        print(seed.start, seed.start_ref, seed.size)
    socs = soc_module.execute(seeds, query, reference_pack, reference_index)
    harm_seeds = harm_module.execute(socs, query)
    alignments = nw_module.execute(harm_seeds, query, reference_pack)

    for index, alignment in enumerate(alignments):
        print("Alignment ", index, ":")
        printer_module.execute(alignment, query, reference_pack)
    
    print("[Test Successful]")
