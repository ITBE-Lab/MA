from LAuS import *
import random

seq = ""

for _ in range(1000):
    char = random.randint(1,4)
    if char == 1:
        seq += "a"
    elif char == 2:
        seq += "c"
    elif char == 3:
        seq += "t"
    else:
        seq += "g"


ref = NucSeq(seq)
fm_index = FM_index(ref).store("test")

print("done")