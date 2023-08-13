from nupack import Tube, Complex, complex_analysis, complex_concentrations, Strand, SetSpec, Model
from math import log10

from params import design_parameters

#TODO consider making function more accessible to scoring user provided seq's by breaking out args?
def nupack_score(sequence:str, score_region:list, design_parameters:design_parameters):
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

    # Does NUPACK tolerate float temperatures? TODO remove if not the case
    cold_temp = int(cold_temp)
    hot_temp = int(hot_temp)

    strand_nucl = Strand(name='A', string=sequence)

    #Create NUPACK complexes for monomer and homodimer
    complex_nucl_single = Complex(strands=[strand_nucl], name='A')
    complex_nucl_double = Complex(strands=[strand_nucl, strand_nucl], name='AA')

    #Create NUPACK Tube and track both the monomer and homodimer complexes
    tube_nucl = Tube(strands={strand_nucl:1e-6}, complexes=SetSpec(max_size=2,
        include=(complex_nucl_single, complex_nucl_double)), name='tube_nucl')
    
    scores_cold = nupack_score_temp(score_region, temp=cold_temp, tube_nucl=tube_nucl,
        complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot=False,
        max_dimer_monomer_factor=design_parameters.max_dimer_monomer_factor,
        nucl_max_score=design_parameters.nucl_max_score)
    
    scores_hot = nupack_score_temp(score_region, temp=hot_temp, tube_nucl=tube_nucl,
        complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot = True,
        max_dimer_monomer_factor=design_parameters.max_dimer_monomer_factor,
        nucl_max_score=design_parameters.nucl_max_score)

    score_energy = nupack_score_energy(temp=design_parameters.thermo_score_temp, energy=design_parameters.target_energy, tube_nucl=tube_nucl, complex_nucl_single=complex_nucl_single, free_energy_max_score=design_parameters.free_energy_max_score)

    return score_energy + sum(scores_cold) + sum(scores_hot)

def nupack_score_energy(temp: int, energy: float, tube_nucl: Tube, complex_nucl_single: Complex, free_energy_max_score:float):
    model_nucl=Model(celsius=temp)
    results_nucl = complex_analysis(complexes = tube_nucl, model=model_nucl, compute=['pairs'])
    #concentrations_nucl = complex_concentrations(tube=tube_nucl, data = results_nucl)

    #Changing max value for free energy score since it can make the resulting RNATs terrible at RBS occlusion
    score_free_energy = (energy - results_nucl.complexes[complex_nucl_single].free_energy) / energy
    if score_free_energy < 0:
        score_free_energy = 0
    if score_free_energy > free_energy_max_score:
        score_free_energy = free_energy_max_score

    return score_free_energy

def nupack_score_temp(score_region: list, temp: int, tube_nucl: Tube, complex_nucl_single: Complex, complex_nucl_double: Complex, hot:bool, max_dimer_monomer_factor:float, nucl_max_score:float):
    #Make NUPACK model
    model_nucl=Model(celsius=temp)

    #Run NUPACK model for temp (Complex analysis)
    results_nucl = complex_analysis(complexes = tube_nucl, model=model_nucl, compute=['pfunc', 'pairs'])
    concentrations_nucl = complex_concentrations(tube=tube_nucl, data = results_nucl)

    #Calculate ratio of AA to A and take log10. lower is better
    if concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_double] == 0:
        dimer_monomer_factor=0
    elif concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_single] == 0:
        dimer_monomer_factor=1
    else:
        dimer_monomer_factor = log10(concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_double] /
            concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_single]) + 2 #+2 means dimer must be 2 factors of 10 less abundant to avoid score penalty
        if dimer_monomer_factor < 0:
            dimer_monomer_factor = 0 #0 is the best possible factor, indicates limited dimer formation
        if dimer_monomer_factor > max_dimer_monomer_factor:
            dimer_monomer_factor = max_dimer_monomer_factor #cap cost of having a poor monomer formation
    
    if len(results_nucl.complexes[complex_nucl_single].pairs.diagonal) != len(score_region):
        raise ValueError
    
    score_nucl = 0
    count_scored_nuc = 0
    for i, x in enumerate(score_region):
        if x:
            score_nucl += results_nucl.complexes[complex_nucl_single].pairs.diagonal[i]
            count_scored_nuc+=1
    score_nucl = score_nucl / count_scored_nuc

    if score_nucl > nucl_max_score:
        score_nucl = nucl_max_score

    if hot:
        score_nucl = nucl_max_score - score_nucl

    return (dimer_monomer_factor, score_nucl)