from nucl import nucl_acid, nucl_set
from Bio.Seq import Seq
from blist import blacklist
from params import design_parameters
from des import design

blist = blacklist(path="blacklist.fasta")
x = design_parameters(blacklist=blist, target=52.5, offset=20, program="NUPACK",
    weights=[16, 16, 16, 16, 16, 16, 16], factor=0.5)
y = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNUAAGGAGGNNNNNNAUGNNN'),
    no_indel=[0]*14+[1]*20, no_mod=[0]*14+[1]*8+[0]*6+[1]*3+[0]*3, score_region=[0]*14+[1]*8+[0]*6+[1]*3+[0]*3, design_parameters=x, is_rna=True)
#TODO: raise an exception if sequence contains RNA bases if not supposed to, vice versa

pool = nucl_set(nucls = [])
for i in range(0, 32):
    pool.append(nucl_acid(sequence=Seq('NNNNNNNNNNNNNNUAAGGAGGNNNNNNAUGNNN'),
        no_indel=[0]*14+[1]*20,
        no_mod=[0]*14+[1]*8+[0]*6+[1]*3+[0]*3,
        score_region=[0]*14+[1]*8+[0]*6+[1]*3+[0]*3,
        design_parameters=x, is_rna=True))

pool.save('test.fasta')

design(design_parameters=x, max_reps=8, current_rep=0, pool=pool, prev_min=4, iter_count=0)