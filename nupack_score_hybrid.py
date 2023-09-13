from nupack import Tube, Complex, complex_analysis, complex_concentrations, Strand, SetSpec, Model
from math import log10
from params import design_parameters

def nupack_score_hybrid(sequence_1:str, score_region_1:list, is_rna_1:bool, sequence_2:str, score_region_2:list, is_rna_2:bool, design_parameters:design_parameters) -> float:
    return 6.0
