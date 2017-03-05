
from LAuS import *

a = Aligner()

a.addModule(Printer())

a.step()


fm_index = FM_index()
fm_index.load("../BioSolution/assemblies/fm-indices/GCA_000001405.22")

for i in range(100, 200):
    print fm_index.at(i) 

