
from nucl import nucl_acid, scoring_parameters
from Bio.Seq import Seq
from blist import blacklist

blist = blacklist(path="blacklist.fasta")
x = scoring_parameters(blacklist=blist, target=50, offset=10, program="NUPACK")
y = nucl_acid(sequence=Seq('ATTATATAA'), no_indel=[0,0,0,0,0,0,0,1,1], no_mod=[0,0,0,0,0,0,0,1,1], scoring_parameters=x)
print(y.score)