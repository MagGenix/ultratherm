from nucl import nucl_set, mutate, nucl_acid, nucl_hybrid
from params import design_parameters
import time
import os
import multiprocessing
from typing import Union
import random
import math

def design(design_parameters:design_parameters, nucl_pool:nucl_set) -> None:
    """A design loop that retains the best nucl_acid's and decrements mutation weights as it runs
    Repetitions are counted in current_rep.
    If max_reps is hit and the weights can be decremented, the weights will be decremented and current_rep reset.
    If max_reps is hit and the weights cannot be decremented, the loop ends.

    Args:
        design_parameters (design_parameters): the design parameters.
        nucl_pool (nucl_set): a nucl_set that will be modified by the loop.
    """
    #Note - now that the temp offset decrementing code is gone, the nucl_acid's are NEVER rescored (the temp never changes). This is less expensive.
    #If for some reason that becomes necessary in the future (I doubt it), it will have to be added back.
    
    ABSOLUTE_MAX_REPS = (design_parameters.max_reps + 1) * (min(design_parameters.weights[0:6])) # type: ignore

    current_rep = 0
    prev_min = 0.0
    iter_count = 0

    parallel_pool = None

    if design_parameters.parallel and design_parameters.program == 'VIENNA':
        parallel_pool = multiprocessing.Pool()

    while (iter_count <= ABSOLUTE_MAX_REPS):
        print("iter_count:\t" + str(iter_count))
        print("current_rep:\t" + str(current_rep))
        print("min weight:\t" + str(min(design_parameters.weights[0:6])))
        print('')

        if iter_count == 0:
            if not os.path.exists("RESULTS"):
                os.makedirs('RESULTS')
            design_parameters.save("RESULTS/" + 'PARAMS_' + time.asctime().replace(":", "_") + '.yml')
            nucl_pool.save("RESULTS/" + "START_" + time.asctime().replace(":", "_") + '.fastq')
        if current_rep == design_parameters.max_reps:
            if design_parameters.can_decrement_weights():
                design_parameters.decrement_weights()
                current_rep = 0
            else:
                break
        
        nucl_pool_size = len(nucl_pool)
        list_new_nucl = list()
        if parallel_pool != None:
            for nucl in nucl_pool.nucls:
                list_nucl = [(nucl, design_parameters)] * design_parameters.num_mutants
                parallel_results = parallel_pool.starmap(func=mutate_and_score, iterable=list_nucl)
                for result in parallel_results:
                    list_new_nucl.append(result)

        else:
            for nucl in nucl_pool.nucls:
                for i in range(0, design_parameters.num_mutants):
                    list_new_nucl.append(mutate_and_score(nucl=nucl, design_parameters=design_parameters))
        
        for mutant in list_new_nucl:
            nucl_pool.append(mutant)
        
        best_nucls = list()
        num_best_nucls = int(nucl_pool_size * 0.15) # Top ~15%
        if num_best_nucls == 0:
            num_best_nucls = 1
        for i in range(0, num_best_nucls): # Only necessarily retain top ~15% of scores (elitist GA)
            best_nucl = nucl_pool.pop(nucl_pool.scores.index(min(nucl_pool.scores)))
            best_nucls.append(best_nucl)

        num_nucls_to_delete = len(nucl_pool) - nucl_pool_size + len(best_nucls)

        weights_list = [math.e**(design_parameters.optimization_rate * x) for x in nucl_pool.scores] # Concave up function and derivative always > 0, higher scores always more likely to be deleted

        for i in range(0, num_nucls_to_delete):
            elem_to_delete = random.choices(nucl_pool.nucls, weights=weights_list)[0]
            index = nucl_pool.nucls.index(elem_to_delete)
            nucl_pool.remove(index)
            weights_list.pop(index)

        for nucl in best_nucls:
            nucl_pool.append(nucl)

        current_min = min(nucl_pool.scores)
        
        if current_min >= prev_min: # type: ignore
            current_rep+=1
        else:
            current_rep = 0
            if design_parameters.can_decrement_weights():
                design_parameters.decrement_weights()
        iter_count+=1
        # TODO consider making the number of rounds to save an intermediate file a parameter
        if iter_count % 50 == 0:
            nucl_pool.save("RESULTS/" + "MID_" +time.asctime().replace(":", "_") + "_w" + str(min(design_parameters.weights[0:6])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fastq')

        prev_min = current_min # type: ignore
    

    nucl_pool.save("RESULTS/" + "END_" +time.asctime().replace(":", "_") + "_w" + str(min(design_parameters.weights[0:7])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fastq')

    if parallel_pool != None:
        parallel_pool.close()
        parallel_pool.terminate()

def mutate_and_score(nucl: Union[nucl_acid, nucl_hybrid], design_parameters: design_parameters) -> Union[nucl_acid, nucl_hybrid]:
    """Useful for parallelization in design(). Accepts a nucl_acid or nucl_hybrid and returns a mutated, scored nucl_acid or nucl_hybrid respectively.

    Args:
        nucl (Union[nucl_acid, nucl_hybrid]): nucl_acid or nucl_hybrid to mutate.
        design_parameters (design_parameters): the design parameters.

    Returns:
        Union[nucl_acid, nucl_hybrid]: new_nucl
    """
    new_nucl = mutate(nucl = nucl, design_parameters = design_parameters)
    new_nucl.fitness_score(design_parameters = design_parameters)
    return new_nucl
