
from LAuS import *

try:
    # load the fm_index for the human genome
    print "loading fm_index for the human genome..."
    fm_index = FM_index()
    fm_index.load("../BioSolution/assemblies/fm-indices/GCA_000001405.22") 
    print "done"

    # load the reversed fm_index for the human genome
    print "loading fm_index for the reversed human genome..."
    rev_fm_index = FM_index()
    rev_fm_index.load("../BioSolution/Application/rev_GCA_000001405.22")
    print "done"


    # load the nucleotide sequence pack for the human genome
    print "loading the human genome pack..."
    refSeq = BWAPack()
    refSeq.load("../BioSolution/Application/rev_GCA_000001405.22")
    print "done"

    #create simulated querry
    print "creating querry sequence..."
    querrySeq = NucSeq()
    querrySeq.append("cgtaactatagaatgatagaatgatatagaatgacgtaactatagaatgatagaatgataatgatagaatgatagaatgacgtaactatagaatgatagaatgattgatagaatgactatagaatgatagaatgatagaataatgatagatga") 
    print "done"

    #create a container for all the required input
    print "setting up input vector..."
    step1 = ContainerVector()

    step1.append(fm_index)
    step1.append(rev_fm_index)
    step1.append(querrySeq)
    step1.append(refSeq)
    print "done"

    print "running the segmentation step..."
    seg = Segmentation(True, True, 10, 100000)
    seg.bSkipLongBWTIntervals = False

    step2 = seg.execute(step1)
    print "done"

    iterator = step2.begin()
    while iterator.exits():
        start = iterator.get().start()
        end = iterator.get().end()
        sequence = ""
        for i in range(start, end):
            sequence = sequence + querrySeq.at(i)
        print "segment: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence
        iterator.next()

    #print "printing segment List:"
    #print step2.toString()
    #print "something in the print function is segfaulting... TODO: fixme"


    print "searching anchor matches..."
    anchors = NlongestIntervalsAsAnchors(2)
    step3 = anchors.execute(step2)
    print "done"

    
    iterator = step3.begin()
    while iterator.exits():
        start = iterator.get().start()
        end = iterator.get().end()
        sequence = ""
        for i in range(start, end):
            sequence = sequence + querrySeq.at(i)
        print "anchor: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence
        iterator.next()

    

    print "test successful"
except Exception as ex:
    print "An Exception occured: "
    print ex
