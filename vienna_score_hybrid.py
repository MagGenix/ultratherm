import RNA
from math import exp, log10
from params import design_parameters
import ctypes

def vienna_score_hybrid(sequence_1:str, score_region_1:list, is_rna_1:bool, score_strand_1:bool,
                        sequence_2:str, score_region_2:list, is_rna_2:bool, score_strand_2:bool,
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
    
    vienna_score_temp(seq1=sequence_1, seq2=sequence_2, is_rna=is_rna, temp=hot_temp)

    return 6.0

def vienna_score_temp(seq1: str, seq2: str, is_rna:bool, temp: float):
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
    bpp = fc_ab.bpp()

    aa = seq1 + "&" + seq1
    fc_aa = RNA.fold_compound(aa, model)
    energy_aa = fc_aa.pf()[1]

    bb = seq2 + "&" + seq2
    fc_bb = RNA.fold_compound(bb, model)
    energy_bb = fc_bb.pf()[1]

    RNA.co_pf_fold("C") # Bug requires pf calculation
    (ab_final, aa_final, bb_final, a_final, b_final) = RNA.get_concentrations(energy_ab, energy_aa, energy_bb, energy_a, energy_b, 1e-6, 1e-6) # TODO put concentrations here
    # TODO make strand concentration a parameter of nucl_acid

    return 6.0

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
    elif score_free_energy > free_energy_max_score:
        score_free_energy = free_energy_max_score
    return score_free_energy
