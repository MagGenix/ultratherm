from nucl import nucl_acid, mutate
from Bio.Seq import Seq
from blist import blacklist
from des import design_parameters

blist = blacklist(path="blacklist.fasta")
x = design_parameters(blacklist=blist, target=52.5, offset=20, program="NUPACK",
    weights=[16, 16, 16, 16, 16, 16, 16], pool_size=32)
y = nucl_acid(sequence=Seq('TAGTAACTTTTGAATAGTGATTCAGGAGGTACTA'),
    no_indel=[0]*17+[1]*17, no_mod=[0]*17+[1]*17, score_region=[0]*17+[1]*17, design_parameters=x, is_rna=False)
#print(y.score)

for i in range(0, 4):
    z = mutate(nucl_acid(sequence=Seq('TAGTAACTTTTGAATAGTGATTCAGGAGGTACTA'),
    no_indel=[0]*17+[1]*17, no_mod=[0]*17+[1]*17, score_region=[0]*17+[1]*17, design_parameters=x, is_rna=False), x)
    print(z.sequence)
    print(z.score)
    #del z