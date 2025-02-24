import RNA
from math import exp, log10
from Ultratherm.params import design_parameters

def vienna_score(sequence:str, score_region:list, is_rna:bool, concentration: float, design_parameters:design_parameters) -> float:
    """Scores a sequence (using ViennaRNA) on score region accessibility, free energy, hairpin count (using Boltzmann ensemble centroid structure) and dimer formation.

    Args:
        sequence (str): the nucleic acid sequence.
        score_region (list): a list of 0 or 1 (int) indicating which region to be assessed for accessibility.
        is_rna (bool): whether the sequence is RNA or DNA.
        concentration (float): The concentration of the strand to be analyzed.
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

    scores_hot = vienna_score_temp(seq=sequence, score_region=score_region, temp=hot_temp,
        nucl_concentration=concentration,
        parasitic_max_order_magnitude=design_parameters.parasitic_max_order_magnitude,
        parasitic_complex_max_score=design_parameters.parasitic_complex_max_score,
        accessibility_max_score=design_parameters.accessibility_max_score, hot=True, is_rna=is_rna)
    scores_cold = vienna_score_temp(seq=sequence, score_region=score_region, temp=cold_temp,
        nucl_concentration=concentration,
        parasitic_max_order_magnitude=design_parameters.parasitic_max_order_magnitude,
        parasitic_complex_max_score=design_parameters.parasitic_complex_max_score,
        accessibility_max_score=design_parameters.accessibility_max_score, hot=False, is_rna=is_rna)

    scores_energy = vienna_score_energy(seq=sequence, temp=design_parameters.thermo_score_temp,
        target_energy=design_parameters.target_energy, max_hairpins=design_parameters.max_hairpins,
        free_energy_max_score=design_parameters.free_energy_max_score, num_hairpins_max_score=design_parameters.num_hairpins_max_score, is_rna=is_rna)

    return sum(scores_energy) + sum(scores_hot) + sum(scores_cold)

# Returns (float: accessibility_score, float: ensemble_energy)
def vienna_score_temp(seq:str, score_region:list, temp: float, nucl_concentration:float, parasitic_max_order_magnitude:float, parasitic_complex_max_score: float, accessibility_max_score: float, hot: bool, is_rna: bool) -> tuple[float, float]:
    """Generate a tuple containing the dimerization score and the score region accessibility score respectively.
    Generally reserved for usage by vienna_score().

    Args:
        seq (str): the nucleic acid sequence.
        score_region (list): a list of 0 or 1 (int) indicating which region to be assessed for accessibility.
        temp (float): the temperature to score at.
        nucl_concentration (float): the concentration of the monomeric nucleic acid (provide the initial concentration).
        parasitic_max_order_magnitude (float): The threshold at which to penalize dimer formation, as -log10([DIMER] / [MONOMER]).
        parasitic_complex_max_score (float): The maximum score penalty for dimer formation.
        accessibility_max_score (float): The maximum score penalty for score region accessibility.
        hot (bool): whether to invert the score region accessibility score.
        is_rna (bool): whether the nucleic acid is RNA or DNA.

    Raises:
        ValueError: _description_

    Returns:
        tuple[float, float]: (parasitic_score, accessibility_score)
    """
    # If DNA needed, needs to be selected before RNA.md() called!
    # NOTE - threadsafety WARNING!
    # inclusion of both RNA/DNA in the same set is not threadsafe as parameters are stored in GLOBALS.
    if not is_rna:
        RNA.params_load_DNA_Mathews1999()
    else:
        RNA.params_load_RNA_Turner2004() # Default parameters for Vienna
    model = RNA.md()
    model.temperature = temp
    model.compute_bpp = 1
    #model.gquad = 1 # ViennaRNA [Bug]: pf_dimer bpp matrix contains values >> 1 #209
    fc = RNA.fold_compound(seq, model)
    monomer_energy = fc.pf()[1]
    
    # NOTE! Vienna bpp array reports pairs from and onto in different positions of array.
    # This is unlike the NUPACK pairs array which is symmetric.
    basepair_probs_diagonal = list()
    bpp = fc.bpp()
    for i in range(1, len(seq) + 1):
        basepair_probs_diagonal.append(0.0)
        for j in range(1, len(seq) + 1):
            basepair_probs_diagonal[i - 1] += bpp[i][j] + bpp[j][i]
    # TODO delete unnecessary check
    if len(basepair_probs_diagonal) != len(score_region):
        raise ValueError
    
    # Array was generated by summing probabilities of each nucl to be paired to another.
    # Needs to reflect the probability it is unpaired (i.e. smaller number results from more pairing)
    for i in range(len(basepair_probs_diagonal)):
        basepair_probs_diagonal[i] = 1 - basepair_probs_diagonal[i]
    
    accessibility_score = 0
    count_scored_nuc = 0
    for i, x in enumerate(score_region):
        if x:
            accessibility_score += basepair_probs_diagonal[i]
            count_scored_nuc+=1
    accessibility_score = accessibility_score / count_scored_nuc

    if accessibility_score > 1.0:
        accessibility_score = 1.0

    dimer_energy = vienna_dimer_energy(seq=seq, temp=temp, is_rna=is_rna)

    # Reaction: 2 Monomer <--> Dimer
    # Energies are in kcal / mol
    delta_g = dimer_energy - (2 * monomer_energy)

    # 1 kcal = 4184 J
    # R: 8.31446261815324 J / (mol * K)
    # Temp: in C, to K is + 273.15
    # delG = -RT lnQ --> Q = e^(-delG / RT)
    # Q = [DIMER] / [MONOMER]^2
    # [DIMER] / [MONOMER] = Q * [MONOMER]
    # log10(Q*[MONOMER]) + parasitic_max_order_magnitude to produce parasitic_score
    # NOTE! This assumes that the FINAL concentration of the [MONOMER] is as provided, not the INITIAL.
    # This means it will not correct [MONOMER] for [DIMER], so it is an approximation for small [DIMER].
    # If this is only used for scoring this is appropriate,
    # as for [DIMER] > [MONOMER] * 10^(1-parasitic_max_order_magnitude) the score is the max score.
    if delta_g >= 0:
        parasitic_score = 1.0
    else:
        parasitic_score = log10(exp((delta_g * -4184) / (8.31446261815324 * (273.15 + temp)) ) * nucl_concentration) + parasitic_max_order_magnitude
        if parasitic_score < 0:
            parasitic_score = 0 #0 is the best possible factor, indicates limited dimer formation
        elif parasitic_score > 1.0:
            parasitic_score = 1.0 #cap cost of having a poor monomer formation

    if hot:
        accessibility_score = 1.0 - accessibility_score
    
    parasitic_score = parasitic_complex_max_score * parasitic_score
    accessibility_score = accessibility_max_score * accessibility_score

    return (parasitic_score, accessibility_score)
    
def vienna_dimer_energy(seq:str, temp:float, is_rna: bool) -> float: # TODO collapse this function into vienna_score_temp - we can have more than one fold_compond!
    """Determines the ensemble free energy of a dimerized sequence.
    Generally reserved for usage by vienna_score().

    Args:
        seq (str): the nucleic acid sequence.
        temp (float): the temperature to score at.
        is_rna (bool): whether the nucleic acid is RNA or DNA.

    Returns:
        float: dimer_energy
    """
    if not is_rna:
        RNA.params_load_DNA_Mathews1999()
    else:
        RNA.params_load_RNA_Turner2004() # Default parameters for Vienna
    model = RNA.md()
    model.temperature = temp
    model.compute_bpp = 0
    #model.gquad = 1 # ViennaRNA [Bug]: pf_dimer bpp matrix contains values >> 1 #209
    fc = RNA.fold_compound(seq + "&" + seq, model)
    dimer_energy = fc.pf()[1]
    return dimer_energy

def vienna_score_energy(seq:str, temp:float, target_energy: float, max_hairpins:int, free_energy_max_score: float, num_hairpins_max_score:float, is_rna: bool) -> tuple[float, float]:
    """Generate a free energy score for a sequence at a given temperature.
    Generally reserved for usage by vienna_score().

    Args:
        seq (str): the nucleic acid sequence.
        temp (float): the temperature to score at.
        target_energy (float): the target free energy in kcal/mol.
        max_hairpins (int): the maximum number of acceptable hairpins in the ensemble centroid structure.
        free_energy_max_score (float): the maximum score penalty for having a free energy greater than target.
        num_hairpins_max_score (float): the maximum score penalty for having more hairpins in the ensemble centroid than the desired max.
        is_rna (bool): whether the nucleic acid is RNA or DNA.

    Returns:
        tuple[float, float]: (score_free_energy, num_hairpins_score)
    """
    if not is_rna:
        RNA.params_load_DNA_Mathews1999()
    else:
        RNA.params_load_RNA_Turner2004() # Default parameters for Vienna
    model = RNA.md()
    model.temperature = temp
    #model.compute_bpp = 0 # comment out to enable centroid calculation
    #model.gquad = 1 # ViennaRNA [Bug]: pf_dimer bpp matrix contains values >> 1 #209
    fc = RNA.fold_compound(seq, model)
    ensemble_energy = fc.pf()[1]
    #(structure, mfe) = fc.mfe()
    # Using centroid, not mfe! centroid better resembles the overall structure of the ensemble
    structure = fc.centroid()[0]

    # Calculate number of hairpins in the dot-bracket structure
    pair_sum = 0
    num_hairpins = 0

    for char in structure:
        if char == ".":
            continue

        if char == "(":
            if pair_sum == 0:
                num_hairpins+=1
            pair_sum += 1
        elif char == ")":
            pair_sum -= 1

    if num_hairpins > max_hairpins:
        num_hairpins_score = num_hairpins_max_score
    else:
        num_hairpins_score = 0

    score_free_energy = (target_energy - ensemble_energy) / target_energy
    if score_free_energy < 0:
        score_free_energy = 0
    elif score_free_energy > 1.0:
        score_free_energy = 1.0
    
    score_free_energy = free_energy_max_score * score_free_energy
    return (score_free_energy, num_hairpins_score)
