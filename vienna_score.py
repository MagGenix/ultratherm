import RNA
from math import exp, log10
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

    score_region_byte_list = score_region
    for i in range(len(score_region_byte_list)):
        if score_region_byte_list[i]:
            score_region_byte_list[i] = '1'
        else:
            score_region_byte_list[i] = '0'
    
    scores_hot = vienna_score_temp(seq=sequence, score_region=score_region, temp=hot_temp,
        max_dimer_monomer_factor=design_parameters.max_dimer_monomer_factor,
        nucl_max_score=design_parameters.nucl_max_score)
    scores_cold = vienna_score_temp(seq=sequence, score_region=score_region, temp=cold_temp,
        max_dimer_monomer_factor=design_parameters.max_dimer_monomer_factor,
        nucl_max_score=design_parameters.nucl_max_score)

    score_energy = vienna_score_energy(seq=sequence, temp=design_parameters.thermo_score_temp,
        target_energy=design_parameters.target_energy,
        free_energy_max_score=design_parameters.free_energy_max_score)

    return score_energy + sum(scores_hot) + sum(scores_cold)

# Returns (float: score_nucl, float: ensemble_energy)
def vienna_score_temp(seq:str, score_region:list, temp: float, max_dimer_monomer_factor: float, nucl_max_score: float) -> tuple[float, float]:
    model = RNA.md()
    model.temperature = temp
    model.compute_bpp = 1
    fc = RNA.fold_compound(seq, model)
    (ss, monomer_energy) = fc.pf()
    
    #NOTE! Vienna bpp array reports pairs from and onto in different positions of array.
    #TODO: verify that the NUPACK array doesn't differentiate this, or scoring may be off!
    basepair_probs = list()
    for i in range(1, len(seq) + 1):
        basepair_probs.append(0.0)
        for j in range(1, len(seq) + 1):
            #Perhaps this can be sped up by copying fc.bpp() to another structure?
            #May be unnecessary. TODO consider changing that
            basepair_probs[i - 1] += fc.bpp()[i][j] + fc.bpp()[j][i]

    if len(basepair_probs) != len(score_region):
        raise ValueError
    
    score_nucl = 0
    count_scored_nuc = 0
    for i, x in enumerate(score_region):
        if x == '1':
            score_nucl += basepair_probs[i]
            count_scored_nuc+=1
    score_nucl = score_nucl / count_scored_nuc

    if score_nucl > nucl_max_score:
        score_nucl = nucl_max_score

    dimer_energy = vienna_dimer_energy(seq=seq, temp=temp)

    # Reaction: 2 Monomer <--> Dimer
    # Energies are in kcal / mol
    delta_g = dimer_energy - (2 * monomer_energy)

    # 1 kcal = 4184 J
    # R: 8.31446261815324 J / (mol * K)
    # Temp: in C, to K is + 273.15
    # delG = RT lnQ --> Q = e^(delG / RT)
    # log10(Q) + 2 to produce dimer_monomer_factor
    dimer_monomer_factor = log10(exp((delta_g * 4184) / (8.31446261815324 * (273.15 + temp)) )) + 2
    if dimer_monomer_factor < 0:
        dimer_monomer_factor = 0 #0 is the best possible factor, indicates limited dimer formation
    if dimer_monomer_factor > max_dimer_monomer_factor:
        dimer_monomer_factor = max_dimer_monomer_factor #cap cost of having a poor monomer formation

    return (dimer_monomer_factor, score_nucl)
    
def vienna_dimer_energy(seq:str, temp:float) -> float:
    model = RNA.md()
    model.temperature = temp
    fc = RNA.fold_compound(seq + "&" + seq, model)
    (ss, dimer_energy) = fc.pf()
    return dimer_energy

def vienna_score_energy(seq:str, temp:float, target_energy: float, free_energy_max_score: float) -> float:
    model = RNA.md()
    model.temperature = temp
    fc = RNA.fold_compound(seq, model)
    (ss, ensemble_energy) = fc.pf()

    score_free_energy = (target_energy - ensemble_energy) / target_energy
    if score_free_energy < 0:
        score_free_energy = 0
    if score_free_energy > free_energy_max_score:
        score_free_energy = free_energy_max_score
    return score_free_energy

# THIS IS FOR TESTING!! TODO get rid of this lmao
from blist import blacklist
def main():
    blist = blacklist(path="blacklist.fasta")
    score = vienna_score(sequence="CCCCCAAAAAGGGGG",
        score_region=[1]*5+[0]*5+[1]*5, design_parameters=design_parameters(blacklist=blist, target_temp=60,
        temp_offset=5, program="VIENNA", weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=1, num_mutants=8,
        target_energy=-12.0, free_energy_max_score=1.0 , nucl_max_score=1.0, max_dimer_monomer_factor=1.0))
    print(score)
main()