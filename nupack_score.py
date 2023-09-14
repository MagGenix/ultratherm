from nupack import Tube, Complex, complex_analysis, complex_concentrations, Strand, SetSpec, Model
from math import log10

from params import design_parameters

def nupack_score(sequence:str, score_region:list, is_rna: bool, design_parameters:design_parameters) -> float:
    """Scores a sequence (using NUPACK) on score region accessibility, free energy, and dimer formation.

    Args:
        sequence (str): the nucleic acid sequence.
        score_region (list): a list of 0 or 1 (int) indicating which region to be assessed for accessibility.
        is_rna (bool): whether the sequence is RNA or DNA.
        design_parameters (design_parameters): The design parameters.

    Raises:
        ValueError: _description_
        Exception: _description_
        Exception: _description_

    Returns:
        float: A total score, lower is better.
    """
    if len(sequence) != len(score_region):
        raise ValueError
    
    if is_rna:
        material = 'rna'
    else:
        material = 'dna'

        #Calculate cold temp for scoring
    cold_temp = design_parameters.target_temp-design_parameters.temp_offset
    if cold_temp < 0:
        cold_temp = 0
    elif cold_temp >= 100:
        raise Exception("illegal cold temperature")
    
    #Calculate hot temp for scoring
    hot_temp = design_parameters.target_temp+design_parameters.temp_offset
    if hot_temp < 0:
        raise Exception("illegal hot temperature")
    elif hot_temp > 100:
        hot_temp = 100

    strand_nucl = Strand(name='A', string=sequence)

    #Create NUPACK complexes for monomer and homodimer
    complex_nucl_single = Complex(strands=[strand_nucl], name='A')
    complex_nucl_double = Complex(strands=[strand_nucl, strand_nucl], name='AA')

    #Create NUPACK Tube and track both the monomer and homodimer complexes
    tube_nucl = Tube(strands={strand_nucl:design_parameters.nucl_concentration}, complexes=SetSpec(max_size=2,
        include=(complex_nucl_single, complex_nucl_double)), name='tube_nucl')
    
    scores_cold = nupack_score_temp(score_region, temp=cold_temp, tube_nucl=tube_nucl,
        complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot=False,
        parasitic_complex_max_score=design_parameters.parasitic_complex_max_score,
        dimer_max_order_magnitude=design_parameters.dimer_max_order_magnitude,
        accessibility_max_score=design_parameters.accessibility_max_score, material=material)
    
    scores_hot = nupack_score_temp(score_region, temp=hot_temp, tube_nucl=tube_nucl,
        complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot = True,
        parasitic_complex_max_score=design_parameters.parasitic_complex_max_score,
        dimer_max_order_magnitude=design_parameters.dimer_max_order_magnitude,
        accessibility_max_score=design_parameters.accessibility_max_score, material=material)

    score_energy = nupack_score_energy(temp=design_parameters.thermo_score_temp, energy=design_parameters.target_energy,
        tube_nucl=tube_nucl, complex_nucl_single=complex_nucl_single,
        free_energy_max_score=design_parameters.free_energy_max_score, material=material)

    return score_energy + sum(scores_cold) + sum(scores_hot)

def nupack_score_energy(
        temp: float, energy: float, tube_nucl: Tube, complex_nucl_single: Complex, free_energy_max_score:float, material: str
    ) -> float:
    """Generate a free energy score for a sequence at a given temperature.
    Generally reserved for usage by nupack_score().
    
    Args:
        temp (float): the temperature to score at.
        energy (float): the target free energy in kcal/mol.
        tube_nucl (Tube): NUPACK Tube containing the strand to be scored.
        complex_nucl_single (Complex): A defined complex of the monomeric strand to be assessed for concentration.
        free_energy_max_score (float): the maximum score penalty for having a free energy greater than target.
        material (str): 'dna' or 'rna'.

    Returns:
        float: score_free_energy
    """
    model_nucl=Model(kelvin=temp + 273.15, material=material)
    results_nucl = complex_analysis(complexes = tube_nucl, model=model_nucl, compute=['pairs'])
    #concentrations_nucl = complex_concentrations(tube=tube_nucl, data = results_nucl)

    #Changing max value for free energy score since it can make the resulting RNATs terrible at RBS occlusion
    score_free_energy = (energy - results_nucl.complexes[complex_nucl_single].free_energy) / energy
    if score_free_energy < 0:
        score_free_energy = 0
    elif score_free_energy > free_energy_max_score:
        score_free_energy = free_energy_max_score

    return score_free_energy

def nupack_score_temp(
        score_region: list, temp: float, dimer_max_order_magnitude:float,
        tube_nucl: Tube, complex_nucl_single: Complex, complex_nucl_double: Complex,
        hot:bool, parasitic_complex_max_score:float, accessibility_max_score:float, material: str
    ) -> tuple[float, float]:
    """Generate a tuple containing the dimerization score and the score region accessibility score respectively.
    Generally reserved for usage by nupack_score().

    Args:
        score_region (list): a list of 0 or 1 (int) indicating which region to be assessed for accessibility.
        temp (float): the temperature to score at.
        dimer_max_order_magnitude (float): The threshold at which to penalize dimer formation, as -log10([DIMER] / [MONOMER]).
        tube_nucl (Tube): NUPACK Tube containing the strand to be scored.
        complex_nucl_single (Complex): A defined complex of the monomeric strand to be assessed for concentration.
        complex_nucl_double (Complex): A defined complex of the dimerized strand to be assessed for concentration.
        hot (bool): whether to invert the score region accessibility score.
        parasitic_complex_max_score (float): The maximum score penalty for dimer formation.
        accessibility_max_score (float): The maximum score penalty for score region accessibility.
        material (str): 'dna' or 'rna'.

    Raises:
        ValueError: _description_

    Returns:
        tuple[float, float]: (dimer_monomer_factor, accessibility_score)
    """
    #Make NUPACK model
    model_nucl=Model(kelvin=temp + 273.15, material=material)

    #Run NUPACK model for temp (Complex analysis)
    results_nucl = complex_analysis(complexes = tube_nucl, model=model_nucl, compute=['pfunc', 'pairs'])
    concentrations_nucl = complex_concentrations(tube=tube_nucl, data = results_nucl)

    nucl_dimer_conc = concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_double]
    nucl_monomer_conc = concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_single]

    #Calculate ratio of AA to A and take log10. lower is better
    if nucl_dimer_conc == 0:
        dimer_monomer_factor=0
    elif nucl_monomer_conc == 0:
        dimer_monomer_factor=1
    else:
        dimer_monomer_factor = log10(nucl_dimer_conc / nucl_monomer_conc) + dimer_max_order_magnitude # +2 means dimer must be 2 factors of 10 less abundant to avoid score penalty
        if dimer_monomer_factor < 0:
            dimer_monomer_factor = 0 #0 is the best possible factor, indicates limited dimer formation
        elif dimer_monomer_factor > parasitic_complex_max_score:
            dimer_monomer_factor = parasitic_complex_max_score #cap cost of having a poor monomer formation
    
    if len(results_nucl.complexes[complex_nucl_single].pairs.diagonal) != len(score_region):
        raise ValueError
    
    accessibility_score = 0
    count_scored_nuc = 0
    diagonal = results_nucl.complexes[complex_nucl_single].pairs.diagonal
    for i, x in enumerate(score_region):
        if x:
            accessibility_score += diagonal[i]
            count_scored_nuc+=1
    accessibility_score = accessibility_score / count_scored_nuc

    if accessibility_score > accessibility_max_score:
        accessibility_score = accessibility_max_score

    if hot:
        accessibility_score = accessibility_max_score - accessibility_score

    return (dimer_monomer_factor, accessibility_score)
