
from LAuS import *

a = Aligner()

a.test()

p = Printer()

a.addModule(p)

a.step()

p.getInputType()

print "generating FM_index object"
fm_index = FM_index()
try:
    fm_index.load("blub")
except:
    print "error - as expected"


print "loading FM_index"
if FM_index.exists("../BioSolution/assemblies/fm-indices/GCA_000001405.22"):
    print "pack availiable"
    fm_index.load("../BioSolution/assemblies/fm-indices/GCA_000001405.22")

    #for i in range(100, 200):
    #    print fm_index.at(i) 
else:
    print "pack not availiable"

print "test successful"
