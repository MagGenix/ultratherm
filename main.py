from nucl import nucl_acid
from Bio.Seq import Seq
from blist import blacklist
from des import design_parameters

blist = blacklist(path="blacklist.fasta")
x = design_parameters(blacklist=blist, target=52.5, offset=20, program="NUPACK", weights=[16, 16, 16, 16, 16, 16, 16, 16])
y = nucl_acid(sequence=Seq('TAGTAACTTTTGAATAGTGATTCAGGAGGTACTA'), no_indel=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], no_mod=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], score_region=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0], design_parameters=x)
print(y.score)