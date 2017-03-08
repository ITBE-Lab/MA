
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
    querrySeq.append("cgtaactatagaatgatagaatgatagaatgatagaatgatagaatga") 
    print "done"

    #create a container for all the required input
    print "setting up input vector..."
    input1 = ContainerVector()

    input1.append(fm_index)
    input1.append(rev_fm_index)
    input1.append(querrySeq)
    input1.append(refSeq)
    print "done"

    print "running the segmentation step..."
    seg = Segmentation(True, True, 3, 10000)
    seg.bSkipLongBWTIntervals = False

    output1 = seg.execute(input1)
    print "done"

    iterator = output1.begin()
    while iterator.exits():
        start = iterator.get().start()
        end = iterator.get().end()
        sequence = ""
        for i in range(start, end):
            sequence = sequence + querrySeq.at(i)
        print "segment: (" + str(start) + "," + str(end) + ") = " + sequence
        iterator.next()

    print "printing segment List:"
    print output1.toString()

    print output1.getTypeInfo()

    print "test successful"
except Exception as ex:
    print "An Exception occured: "
    print ex
