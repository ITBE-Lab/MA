from .aligner import *
from .alignmentPrinter import *
from random import choice

def test_aligner():
    parameter_manager = ParameterSetManager()
    parameter_manager.set_selected("PacBio")

    # create a reference string
    reference = ""
    for _ in range(1000000):
        reference += choice(['C', 'T', 'G'])

    # create pack and fmd index from that string
    reference_pack = Pack()
    reference_pack.append("sequence name", "sequence description", NucSeq(reference))
    reference_index = FMIndex(reference_pack)

    query = NucSeq(reference[100:125] + "A" + reference[126:150] + "A" + reference[151:175] + "A" + reference[176:200])

    # initialize modules
    seeding_module = BinarySeeding(parameter_manager)
    soc_module = StripOfConsideration(parameter_manager)
    harm_module = Harmonization(parameter_manager)
    nw_module = NeedlemanWunsch(parameter_manager)
    printer_module = AlignmentPrinter(parameter_manager)

    #execute modules
    seeds = seeding_module.execute(reference_index, query)
    print("seed set (q,r,l):")
    for seed in seeds.extract_seeds(reference_index, 100, 0, True):
        print(seed.start, seed.start_ref, seed.size)
    socs = soc_module.execute(seeds, query, reference_pack, reference_index)
    harm_seeds = harm_module.execute(socs, query)
    for i, seeds in enumerate(harm_seeds):
        print("harmonization seed set ", i, ": (q,r,l):")
        for seed in seeds:
            print(seed.start, seed.start_ref, seed.size)
    alignments = nw_module.execute(harm_seeds, query, reference_pack)
    for index, alignment in enumerate(alignments):
        print("Alignment ", index, ":")
        printer_module.execute(alignment, query, reference_pack)
    
    print("[Test Successful]")
