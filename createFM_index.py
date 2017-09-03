from LAuS import *

seq = NucSeq("cgtaactatagaatgatagaatgatatagaatgacgtaactatagaatga")
fm_index = FM_index(seq).save("test")

print("done")