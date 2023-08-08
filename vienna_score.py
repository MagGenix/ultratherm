from ctypes import *

from params import design_parameters

def vienna_score(sequence:str, score_region:list, design_parameters:design_parameters):
    if len(sequence) != len(score_region):
        raise ValueError
    
        #Calculate cold temp for scoring
    cold_temp = design_parameters.target_temp-design_parameters.temp_offset
    if cold_temp < 0:
        cold_temp = 0
    if cold_temp >= 100:
        raise Exception("illegal cold temperature")
    
    #Calculate hot temp for scoring
    hot_temp = design_parameters.target_temp+design_parameters.temp_offset
    if hot_temp < 0:
        raise Exception("illegal hot temperature")
    if hot_temp > 100:
        hot_temp = 100

    seq_c = c_char_p(bytes('GATTACA', 'ascii')) # Operating systems should understand ASCII
    
    score_region_byte_list = score_region
    for i in range(len(score_region_byte_list)):
        if score_region_byte_list[i]:
            score_region_byte_list[i] = '1'
        else:
            score_region_byte_list[i] = '0'
    
    score_region_bytes_as_str = ''.join(score_region_byte_list)
    score_reg_c = c_char_p(bytes(score_region_bytes_as_str, 'ascii'))

    vienna_score_functions = CDLL("vienna_score.so")
    # Compilation - gcc -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -lRNA -lpthread -lmpfr -lgmp -lstdc++ -fPIC -shared -o vienna_score.so vienna_score.c
    # For non-shared, gcc -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -lRNA -lpthread -lmpfr -lgmp -lstdc++ -o vienna_score.o vienna_score.c
    vienna_score_functions.score(seq_c, score_reg_c)


# THIS IS FOR TESTING!! TODO get rid of this lmao
from blist import blacklist
def main():
    blist = blacklist(path="blacklist.fasta")
    vienna_score(sequence="GATTACA", score_region=[0,1,1,0,1,0,1], design_parameters=design_parameters(blacklist=blist, target_temp=70, temp_offset=5, program="VIENNA",
    weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=1, num_mutants=8, target_energy=-12.0, # based on FourU Hairpin 2
    free_energy_max_score=1.0 , nucl_max_score=1.0, max_dimer_monomer_factor=1.0))
main()