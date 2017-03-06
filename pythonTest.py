
from LAuS import *

# load the fm_index for the human genome
fm_index = FM_index()
fm_index.load("../BioSolution/assemblies/fm-indices/GCA_000001405.22") 

# load the reversed fm_index for the human genome
rev_fm_index = FM_index()
rev_fm_index.load("../BioSolution/Application/rev_GCA_000001405.22")


# load the nucleotide sequence pack for the human genome
refSeq = BWAPack()
refSeq.load("../BioSolution/Application/rev_GCA_000001405.22")

#create simulated querry
querrySeq = NucSeq()
querrySeq.append("cgtaactatagaatga") 

#create a container for all the required input
input = ContainerVector()


a = Aligner()

print "test successful"
