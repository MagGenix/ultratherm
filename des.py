from nucl import nucl_set, mutate
from params import design_parameters
import time
import os

def design(design_parameters:design_parameters, pool:nucl_set, current_rep:int = 0, prev_min:float = 0.0, iter_count:int = 0) -> None:
    """A recursive design loop that retains the best nucl_acid's and decrements mutation weights as it runs
    Repetitions are counted in current_rep.
    If max_reps is hit and the weights can be decremented, the weights will be decremented and current_rep reset.
    If max_reps is hit and the weights cannot be decremented, the loop ends.

    Args:
        design_parameters (design_parameters): the design parameters.
        current_rep (int): the repetition the preceding loop (or main) started on.
        pool (nucl_set): a nucl_set that will be modified by the loop.
        prev_min (float): the minimum pool score from the preceding loop.
        iter_count (int): the total number of repetitions.
    """
    #Note - now that the temp offset decrementing code is gone, the nucl_acid's are NEVER rescored (the temp never changes). This is less expensive.
    #If for some reason that becomes necessary in the future (I doubt it), it will have to be added back.
    
    if iter_count == 0:
        if not os.path.exists("RESULTS"):
            os.makedirs('RESULTS')
        design_parameters.save("RESULTS/" + 'PARAMS_' + time.asctime() + '.yml')
        pool.save("RESULTS/" + "START_" + time.asctime() + '.fastq')
    if current_rep == design_parameters.max_reps:
        if design_parameters.can_decrement_weights():
            design_parameters.decrement_weights()
            current_rep = 0
        else:
            pool.save("RESULTS/" + "END_" +time.asctime() + "_w" + str(min(design_parameters.weights[0:7])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fastq')
            return
    for nucl in pool.nucls:
        for i in range(0, design_parameters.num_mutants):
            pool.append(mutate(nucl=nucl, design_parameters=design_parameters))
        
        for i in range(0, design_parameters.num_mutants):
            pool.remove(pool.scores.index(max(pool.scores)))

    current_min = min(pool.scores)
    
    if current_min >= prev_min: # type: ignore
        current_rep+=1
    else:
        current_rep = 0
        if design_parameters.can_decrement_weights():
            design_parameters.decrement_weights()
    iter_count+=1
    # TODO consider making the number of rounds to save an intermediate file a parameter
    if iter_count % 50 == 0:
        pool.save("RESULTS/" + "MID_" +time.asctime() + "_w" + str(min(design_parameters.weights[0:6])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fastq')
    
    print("iter_count:\t" + str(iter_count))
    print("current_rep:\t" + str(current_rep))
    print("min weight:\t" + str(min(design_parameters.weights[0:6])))
    print('')
    design(design_parameters=design_parameters, current_rep=current_rep, pool=pool, prev_min=current_min, iter_count=iter_count) # type: ignore
