import RNA
from math import exp, log10
from params import design_parameters

def vienna_score_hybrid(sequence_1:str, score_region_1:list, is_rna_1:bool, score_strand_1:bool,
                        sequence_2:str, score_region_2:list, is_rna_2:bool, score_strand_2:bool,
                        design_parameters:design_parameters) -> float:
    return 6.0
