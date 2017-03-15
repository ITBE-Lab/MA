
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
    vector1 = ContainerVector()

    vector1.append(fm_index)
    vector1.append(rev_fm_index)
    vector1.append(querrySeq)
    vector1.append(refSeq)
    print "done"

    print "running the segmentation step..."
    seg = Segmentation(True, True, 10, 100000)
    seg.bSkipLongBWTIntervals = False

    segments = seg.execute(vector1)
    print "done"

    iterator = segments.begin()
    while iterator.exits():
        start = iterator.get().start()
        end = iterator.get().end()
        sequence = ""
        for i in range(start, end):
            sequence = sequence + querrySeq.at(i)
        print "segment: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence
        iterator.next()

    #print "printing segment List:"
    #print segments.toString()
    #print "something in the print function is segfaulting... TODO: fixme"


    print "searching anchor matches..."
    anc = NlongestIntervalsAsAnchors(2)
    anchors = anc.execute(segments)
    print "done"

    
    iterator = anchors.begin()
    while iterator.exits():
        start = iterator.get().start()
        end = iterator.get().end()
        sequence = ""
        for i in range(start, end):
            sequence = sequence + querrySeq.at(i)
        print "anchor: (" + str(start) + "," + str(end) + ";" + str(end - start) + ") := " + sequence
        iterator.next()

    
    print "setting up input vector..."
    vector2 = ContainerVector()
    vector2.append(segments)
    vector2.append(anchors)
    vector2.append(querrySeq)
    vector2.append(refSeq)
    vector2.append(fm_index)
    vector2.append(rev_fm_index)
    print "done"

    bucketing = Bucketing()

    print "collecting strips of consideration..."
    stripsOfConsideration = bucketing.execute(vector2)
    print "done"

    for i in range(0,stripsOfConsideration.size()):
        strip = stripsOfConsideration.at(i)
        print "strip of consideration (" + str(i) + "): " + str(strip.getScore())

    print "test successful"
except Exception as ex:
    print "An Exception occured: "
    print ex
