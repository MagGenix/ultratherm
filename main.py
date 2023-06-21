from nucl import nucl_acid, mutate
from Bio.Seq import Seq
from blist import blacklist
from params import design_parameters
from des import design
import time

blist = blacklist(path="blacklist.fasta")
x = design_parameters(blacklist=blist, target=52.5, offset=20, program="NUPACK",
    weights=[16, 16, 16, 16, 16, 16, 16], factor=0.5)
y = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNTAGTGATTCAGGAGGTACTA'),
    no_indel=[0]*17+[1]*17, no_mod=[0]*17+[1]*17, score_region=[0]*17+[1]*17, design_parameters=x, is_rna=False)
#print(y.score)

print(time.time())
for i in range(0, 127):
    z = mutate(y, x)
    #print(str(z.score) + "\t" + str(z.sequence))
    #print(z.score)
    #del z

print(time.time())
#Takes around 7 seconds for 128 mutants
#design(design_parameters=x, )