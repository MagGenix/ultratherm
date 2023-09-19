import RNA
from math import log10
from params import design_parameters
import numpy

def vienna_score_hybrid(sequence_1:str, score_region_1:list, is_rna_1:bool, score_strand_1:bool, concentration_1:float,
                        sequence_2:str, score_region_2:list, is_rna_2:bool, score_strand_2:bool, concentration_2:float,
                        design_parameters:design_parameters) -> float:
    if len(sequence_1) != len(score_region_1) or len(sequence_2) != len(score_region_2):
        raise ValueError
    if is_rna_1 ^ is_rna_2: # Can only fold hybrids with the same strand types!
        raise ValueError
    if is_rna_1 and is_rna_2:
        is_rna = True
    else:
        is_rna = False

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
    
    scores_cold = vienna_score_temp(
        seq1=sequence_1, seq2=sequence_2, is_rna=is_rna, temp=hot_temp,
        parasitic_complex_max_score=design_parameters.parasitic_complex_max_score,
        accessibility_max_score=design_parameters.accessibility_max_score,
        hybrid_max_score=design_parameters.hybrid_max_score,
        parasitic_max_order_magnitude=design_parameters.parasitic_max_order_magnitude,
        score_region_1=score_region_1, score_region_2=score_region_2,
        score_strand_1=score_strand_1, score_strand_2=score_strand_2, hot=False,
        concentration_1=concentration_1, concentration_2=concentration_2)

    scores_hot = vienna_score_temp(
        seq1=sequence_1, seq2=sequence_2, is_rna=is_rna, temp=hot_temp,
        parasitic_complex_max_score=design_parameters.parasitic_complex_max_score,
        accessibility_max_score=design_parameters.accessibility_max_score,
        hybrid_max_score=design_parameters.hybrid_max_score,
        parasitic_max_order_magnitude=design_parameters.parasitic_max_order_magnitude,
        score_region_1=score_region_1, score_region_2=score_region_2,
        score_strand_1=score_strand_1, score_strand_2=score_strand_2, hot=True,
        concentration_1=concentration_1, concentration_2=concentration_2)

    score_free_energy = vienna_score_energy(
        seq1=sequence_1, seq2=sequence_2,
        temp=design_parameters.thermo_score_temp,
        target_energy=design_parameters.target_energy,
        free_energy_max_score=design_parameters.free_energy_max_score,
        is_rna=is_rna)
    
    return sum(scores_cold) + sum(scores_hot) + score_free_energy

def vienna_score_temp(seq1: str, seq2: str,
                      is_rna:bool, temp: float,
                      parasitic_complex_max_score: float, accessibility_max_score: float, hybrid_max_score: float,
                      parasitic_max_order_magnitude: float,
                      score_region_1:list, score_region_2:list,
                      score_strand_1: bool, score_strand_2:bool, hot:bool,
                      concentration_1:float, concentration_2:float) -> tuple[float, float, float]:
    if not is_rna:
        RNA.params_load_DNA_Mathews1999()
    else:
        RNA.params_load_RNA_Turner2004()
    model = RNA.md()
    model.temperature = temp
    model.compute_bpp = 1
    model.gquad = 1

    ab = seq1 + "&" + seq2
    fc_ab = RNA.fold_compound(ab, model)
    (energy_a, energy_b, energy_ab) = fc_ab.pf_dimer()[1:4]
    bpp_tuple = fc_ab.bpp()
    bpp = numpy.array(bpp_tuple)[1:,1:]

    aa = seq1 + "&" + seq1
    fc_aa = RNA.fold_compound(aa, model)
    energy_aa = fc_aa.pf()[1]

    bb = seq2 + "&" + seq2
    fc_bb = RNA.fold_compound(bb, model)
    energy_bb = fc_bb.pf()[1]

    RNA.co_pf_fold("C") # Bug requires pf calculation
    (hybrid_concentration, aa_final, bb_final, a_final, b_final) = RNA.get_concentrations(energy_ab, energy_aa, energy_bb, energy_a, energy_b, concentration_1, concentration_2)

    total_unbound_concentration = a_final + b_final
    total_parasitic_concentration = aa_final + bb_final

    parasitic_score = 1.0
    hybrid_score = 1.0
    accessibility_score = 1.0

    if total_unbound_concentration == 0 and hybrid_concentration == 0: # Worst case - all parasitic, no unbound and no hybrid
        parasitic_score = 1.0
    elif total_parasitic_concentration == 0:
        parasitic_score = 0.0
    else:
        parasitic_score = log10(
            total_parasitic_concentration / (total_unbound_concentration + hybrid_concentration)
            ) + parasitic_max_order_magnitude
        if parasitic_score < 0:
            parasitic_score = 0 #0 is the best possible factor, indicates limited dimer formation
        elif parasitic_score > 1.0:
            parasitic_score = 1.0 #cap cost of having a poor monomer formation
    
    if hybrid_concentration == 0:
        hybrid_score = 1.0
    elif total_unbound_concentration == 0:
        hybrid_score = 0.0
    else:
        # NOTE - in this implementation, hybrid_score and accessibility_score are given the same caps.
        # This is because they vary together.
        hybrid_score = (log10(total_unbound_concentration / hybrid_concentration) + 0.5) * 31.62278 # For min score, needs a change of 10^2
        if hybrid_score < 0:
            hybrid_score = 0
        elif hybrid_score > 1.0:
            hybrid_score = 1.0 #cap cost

    sub_pairs_arr_1 = bpp[0:len(score_region_1), len(score_region_1):]
    sub_pairs_arr_2 = bpp[len(score_region_1):, 0:len(score_region_1)]

    sub_pairs_arr = sub_pairs_arr_1
    
    for i, row in enumerate(sub_pairs_arr_2):
        for j, elem in enumerate(row):
            sub_pairs_arr[j,i] += elem

    paired_strand_1 = numpy.sum(sub_pairs_arr, axis=0)
    paired_strand_2 = numpy.sum(sub_pairs_arr, axis=1)
    
    total_bound_1 = 0.0
    count_scored_nuc_1 = 0

    total_bound_2 = 0.0
    count_scored_nuc_2 = 0

    if score_strand_1:
        for i, x in enumerate(score_region_1): # NOTE! Higher pair probs mean MORE pairing, not less as in the diagonal!
            if x:
                total_bound_1 += paired_strand_1[i]
                count_scored_nuc_1+=1
    if score_strand_2:
        for i, x in enumerate(score_region_2):
            if x:
                total_bound_2 += paired_strand_2[i]
                count_scored_nuc_2+=1
    
    accessibility_score = (total_bound_1 + total_bound_2) / float(count_scored_nuc_1 + count_scored_nuc_2)
    if accessibility_score > 1.0:
        accessibility_score = 1.0

    if hot:
        hybrid_score = 1.0 - hybrid_score
        monomeric_accessibility_score_1 = 0.0
        monomeric_accessibility_score_2 = 0.0

        if score_strand_1:
            monomeric_accessibility_score_1 = vienna_score_monomeric_accessibility(seq=seq1, model=model, score_region=score_region_1)
        
        if score_strand_2:
            monomeric_accessibility_score_2 = vienna_score_monomeric_accessibility(seq=seq2, model=model, score_region=score_region_2)
        
        if monomeric_accessibility_score_1 > 1.0:
            monomeric_accessibility_score_1 = 1.0
        if monomeric_accessibility_score_2 > 1.0:
            monomeric_accessibility_score_2 = 1.0
        
        monomeric_accessibility_score_1 = 1.0 - monomeric_accessibility_score_1
        monomeric_accessibility_score_2 = 1.0 - monomeric_accessibility_score_2

        accessibility_score += monomeric_accessibility_score_1 + monomeric_accessibility_score_2

    else:
        accessibility_score = 1.0 - accessibility_score

    parasitic_score = parasitic_complex_max_score * parasitic_score
    hybrid_score = hybrid_max_score * hybrid_score
    accessibility_score = accessibility_max_score * accessibility_score
    return (parasitic_score, hybrid_score, accessibility_score)

def vienna_score_energy(seq1:str, seq2:str, temp:float, target_energy: float, free_energy_max_score: float, is_rna: bool) -> float:
    if not is_rna:
        RNA.params_load_DNA_Mathews1999()
    else:
        RNA.params_load_RNA_Turner2004() # Default parameters for Vienna
    model = RNA.md()
    model.temperature = temp
    model.compute_bpp = 0
    model.gquad = 1
    seq = seq1 + "&" + seq2
    fc = RNA.fold_compound(seq, model)
    ensemble_energy = fc.pf()[1]

    score_free_energy = (target_energy - ensemble_energy) / target_energy
    if score_free_energy < 0:
        score_free_energy = 0
    elif score_free_energy > 1.0:
        score_free_energy = 1.0
    
    score_free_energy = free_energy_max_score * score_free_energy
    return score_free_energy

def vienna_score_monomeric_accessibility(seq: str, model: RNA.md, score_region:list):
    fc = RNA.fold_compound(seq, model)
    basepair_probs_diagonal = list()
    fc.pf()
    bpp = fc.bpp()
    for i in range(1, len(seq) + 1):
        basepair_probs_diagonal.append(0.0)
        for j in range(1, len(seq) + 1):
            basepair_probs_diagonal[i - 1] += bpp[i][j] + bpp[j][i]
    for i in range(len(basepair_probs_diagonal)):
        basepair_probs_diagonal[i] = 1 - basepair_probs_diagonal[i]
    monomeric_accessibility_score = 0
    count_scored_nuc = 0
    for i, x in enumerate(score_region):
        if x:
            monomeric_accessibility_score += basepair_probs_diagonal[i]
            count_scored_nuc+=1
    
    monomeric_accessibility_score = monomeric_accessibility_score / count_scored_nuc

    return monomeric_accessibility_score
