from nupack import Tube, Complex, complex_analysis, complex_concentrations, Strand, SetSpec, Model
from math import log10
from params import design_parameters

def nupack_score_hybrid(sequence_1:str, score_region_1:list, is_rna_1:bool, sequence_2:str, score_region_2:list, is_rna_2:bool, design_parameters:design_parameters) -> float:
    if len(sequence_1) != len(score_region_1) or len(sequence_2) != len(score_region_2):
        raise ValueError
    if not(is_rna_1 and is_rna_2): # Can only fold hybrids with the same strand types!
        raise ValueError
    if is_rna_1 and is_rna_2:
        material = 'rna'
    else:
        material = 'dna'

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


    strand_A = Strand(name='A', string=sequence_1)
    strand_B = Strand(name='B', string=sequence_2)

    complex_A = Complex(strands=[strand_A], name='A')
    complex_B = Complex(strands=[strand_B], name='B')
    complex_AA = Complex(strands=[strand_A, strand_A], name='AA')
    complex_AB = Complex(strands=[strand_A, strand_B], name='AB')
    complex_BB = Complex(strands=[strand_B, strand_B], name='BB')

    unbound_complexes = list()
    unbound_complexes.append(complex_A)
    unbound_complexes.append(complex_B)

    parasitic_complexes = list()
    parasitic_complexes.append(complex_AA)
    parasitic_complexes.append(complex_BB)

    # TODO make strand concentration a parameter of each nucl_acid and have it passed in as an arg to this fxn
    tube_nucl = Tube(strands={strand_A:design_parameters.nucl_concentration, strand_B:design_parameters.nucl_concentration},
                     complexes=SetSpec(max_size=2, include=(complex_A, complex_B, complex_AA, complex_AB, complex_BB)), name='tube_nucl')
    
    scores_cold = nupack_score_temp(
        material=material, temp=cold_temp, tube_nucl=tube_nucl,
        unbound_complexes=unbound_complexes, parasitic_complexes=parasitic_complexes,
        hybrid_complex=complex_AB,
        hot=False, parasitic_max_order_magnitude=design_parameters.dimer_max_order_magnitude)

    return sum(scores_cold)
def nupack_score_temp(
        material: str, temp:float, tube_nucl: Tube,
        unbound_complexes: list[Complex], parasitic_complexes: list[Complex], hybrid_complex: Complex,
        hot: bool, parasitic_max_order_magnitude:float
) -> tuple[float, float, float]:
    
    model_nucl = Model(kelvin=temp + 273.15, material=material)
    results_nucl = complex_analysis(complexes=tube_nucl, model=model_nucl, compute=['pairs'])
    concentrations_nucl = complex_concentrations(tube=tube_nucl, data = results_nucl)

    total_unbound_concentration = 0.0
    total_parasitic_concentration = 0.0
    
    for complex in unbound_complexes:
        total_unbound_concentration+=concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex]
    for complex in parasitic_complexes:
        total_parasitic_concentration+=concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex]
    hybrid_concentration=concentrations_nucl.tubes[tube_nucl].complex_concentrations[hybrid_complex]

    parasitic_score = 1.0
    hybrid_score = 1.0
    accessibility_score = 1.0

    if total_unbound_concentration == 0 and hybrid_concentration == 0: # Worst case - all parasitic, no unbound and no hybrid
        return(1.0, 1.0, 1.0)
    
    if total_parasitic_concentration == 0:
        parasitic_score = 0.0
    else:
        parasitic_score = log10(
            total_parasitic_concentration / (total_unbound_concentration + hybrid_concentration)
            ) + parasitic_max_order_magnitude
    
    if hybrid_concentration == 0:
        hybrid_score = 1.0
    elif total_unbound_concentration == 0:
        hybrid_score = 0.0
    else:
        hybrid_score = log10(
            total_unbound_concentration / hybrid_concentration
        ) + parasitic_max_order_magnitude # TODO Need sep. variable!

    total_bound = 0.0
    count_scored_nuc = 0
    diagonal = results_nucl.complexes[hybrid_complex].pairs.diagonal

    # TODO finish implementation
    # TODO add caps for score penalties


    return(parasitic_score, hybrid_score, accessibility_score)
    

    