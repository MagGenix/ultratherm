from nupack import Tube, Complex, complex_analysis, complex_concentrations, Strand, SetSpec, Model
from math import log10

from params import design_parameters

def nupack_score(sequence:str, score_region:list, is_rna: bool, design_parameters:design_parameters) -> float:
    return 6.0