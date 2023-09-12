import RNA
from math import exp, log10
from params import design_parameters

def vienna_score(sequence:str, score_region:list, is_rna:bool, design_parameters:design_parameters) -> float:
    return 6.0