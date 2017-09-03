from LAuS import *

seq = NucSeq("cgtaactatagaatgatagaatgatatagaatgacgtaactatagaatga")
fm_index = FM_index(seq).store("test")

print("done")