#from ctypes import *
import RNA

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

    #seq_c = c_char_p(bytes('GATTACA', 'ascii')) # Operating systems should understand ASCII
    
    score_region_byte_list = score_region
    for i in range(len(score_region_byte_list)):
        if score_region_byte_list[i]:
            score_region_byte_list[i] = '1'
        else:
            score_region_byte_list[i] = '0'
    
    #score_region_bytes_as_str = ''.join(score_region_byte_list)
    #score_reg_c = c_char_p(bytes(score_region_bytes_as_str, 'ascii'))

    vienna_score_temp(seq=sequence, score_region=score_region, temp=hot_temp)
    vienna_score_temp(seq=sequence, score_region=score_region, temp=cold_temp)

    vienna_score_dimer(seq=sequence, temp=hot_temp)
    vienna_score_dimer(seq=sequence, temp=cold_temp)

    vienna_score_energy(seq=sequence, temp=design_parameters.thermo_score_temp)
    #vienna_score_functions = CDLL("vienna_score.so")
    # Compilation - gcc -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -lRNA -lpthread -lmpfr -lgmp -lstdc++ -fPIC -shared -o vienna_score.so vienna_score.c
    # For non-shared, gcc -pthread -I/usr/local/include -I/usr/local/include/ViennaRNA -L/usr/local/lib -lRNA -lpthread -lmpfr -lgmp -lstdc++ -o vienna_score.o vienna_score.c
    #vienna_score_functions.score(seq_c, score_reg_c)

def vienna_score_temp(seq:str, score_region:list, temp: float):
    model = RNA.md()
    model.temperature = temp
    model.compute_bpp = 1
    fc = RNA.fold_compound(seq, model)
    (ss, ensemble_energy) = fc.pf()
    
    basepair_probs = list()
    for i in range(1, len(seq) + 1):
        basepair_probs.append(0.0)
        for j in range(1, len(seq) + 1):
            basepair_probs[i - 1] += fc.bpp()[i][j] + fc.bpp()[j][i]

    test = fc.bpp()


    if len(basepair_probs) != len(score_region):
        raise ValueError
    
    score_nucl = 0
    count_scored_nuc = 0
    for i, x in enumerate(score_region):
        if x == '1':
            score_nucl += basepair_probs[i]
            count_scored_nuc+=1
    score_nucl = score_nucl / count_scored_nuc

    print(score_nucl)
    print(ensemble_energy)
    
def vienna_score_dimer(seq:str, temp:float):
    model = RNA.md()
    model.temperature = temp
    fc = RNA.fold_compound(seq + "&" + seq, model)
    (ss, complex_energy) = fc.pf()
    print("COMPLEX:\t" + str(complex_energy))

def vienna_score_energy(seq:str, temp:float):
    model = RNA.md()
    model.temperature = temp
    fc = RNA.fold_compound(seq, model)
    (ss, ensemble_energy) = fc.pf()
    print(ensemble_energy)

# THIS IS FOR TESTING!! TODO get rid of this lmao
from blist import blacklist
def main():
    blist = blacklist(path="blacklist.fasta")
    vienna_score(sequence="CCCCCAAAAAGGGGG",
        score_region=[1]*5+[0]*5+[1]*5, design_parameters=design_parameters(blacklist=blist, target_temp=60,
        temp_offset=5, program="VIENNA", weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=1, num_mutants=8,
        target_energy=-12.0, free_energy_max_score=1.0 , nucl_max_score=1.0, max_dimer_monomer_factor=1.0))
    #vienna_score(sequence="GAATTC",
    #   score_region=[1]*6, design_parameters=design_parameters(blacklist=blist, target_temp=60,
    #   temp_offset=5, program="VIENNA", weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=1, num_mutants=8,
    #   target_energy=-12.0, free_energy_max_score=1.0 , nucl_max_score=1.0, max_dimer_monomer_factor=1.0))
main()