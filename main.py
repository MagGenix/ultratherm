from nucl import nucl_acid, nucl_set, blacklist
from Bio.Seq import Seq
from params import design_parameters    
from des import design
import time

#Configure design parameters
blist = blacklist(path="blacklist.fasta")
des_params = design_parameters(blacklist=blist, target_temp=60, temp_offset=5, program="NUPACK",
    weights=[8, 8, 8, 8, 8, 8, 16], weight_factor=1, num_mutants=8, target_energy=-5.0,
    free_energy_max_score=0.25 , nucl_max_score=1.0, max_dimer_monomer_factor=1.0)
des_params.save('PARAMS_' + time.asctime() + '.yml')
#TODO: raise an exception if sequence contains RNA bases if not supposed to, vice versa

#Create nucleotide set
pool = nucl_set(nucls = [])
for i in range(0, 4):
    pool.append(nucl_acid(sequence=Seq('NNNNNNNNNNNNNNUAAGGAGGNNNNNNAUGNNNNNNNNNNNNNN'),
        no_indel =      [0]*14+[1]*17+[0]*14,
        no_mod =        [0]*14+[1]*8+[0]*6+[1]*3+[0]*14,
        score_region =  [0]*14+[1]*8+[0]*6+[1]*3+[0]*14,
        design_parameters=des_params, is_rna=True))

pool.save("START_" + time.asctime() + '.fasta')

#Start design loop
design(design_parameters=des_params, max_reps=16, current_rep=0, pool=pool, prev_min=4, iter_count=0)