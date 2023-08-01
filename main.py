from nucl import nucl_acid, nucl_set, blacklist
from Bio.Seq import Seq
from params import design_parameters    
from des import design
import time

#Configure design parameters
blist = blacklist(path="blacklist.fasta")
des_params = design_parameters(blacklist=blist, target_temp=70, temp_offset=5, program="NUPACK",
    weights=[8, 8, 8, 8, 8, 8, 16], weight_factor=1, num_mutants=8, target_energy=-12.0, # based on FourU Hairpin 2
    free_energy_max_score=1.0 , nucl_max_score=1.0, max_dimer_monomer_factor=1.0)
des_params.save('PARAMS_' + time.asctime() + '.yml')
#TODO: raise an exception if sequence contains RNA bases if not supposed to, vice versa

#Create nucleotide set
pool = nucl_set(nucls = [])
for i in range(0, 32):
    pool.append(nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNNUAAGGAGGNNNNNNAUG'),
        no_indel =      [0]*20+[1]*17,
        no_mod =        [0]*20+[1]*8+[0]*6+[1]*3,
        score_region =  [0]*20+[1]*8+[0]*6+[0]*3,
        design_parameters=des_params, is_rna=True))

pool.save("START_" + time.asctime() + '.fasta')

#Start design loop
design(design_parameters=des_params, max_reps=16, current_rep=0, pool=pool, prev_min=4, iter_count=0)