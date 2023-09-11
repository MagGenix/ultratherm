from nucl import nucl_set, mutate, score
from params import design_parameters
import time
import os
import multiprocessing

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
            design_parameters.save("RESULTS/" + 'PARAMS_' + time.asctime() + '.yml')
            nucl_pool.save("RESULTS/" + "START_" + time.asctime() + '.fastq')
        if current_rep == design_parameters.max_reps:
            if design_parameters.can_decrement_weights():
                design_parameters.decrement_weights()
                current_rep = 0
            else:
                break
        
        if parallel_pool != None:
            for nucl in nucl_pool.nucls:
                list_mut_nucl = list()
                for i in range(design_parameters.num_mutants):
                    list_mut_nucl.append((mutate(nucl=nucl, design_parameters=design_parameters), design_parameters))

                list_new_nucl = parallel_pool.starmap(func=score, iterable=list_mut_nucl)
                
                for mutant in list_new_nucl:
                    nucl_pool.append(mutant)
                    
                for i in range(0, design_parameters.num_mutants):
                    nucl_pool.remove(nucl_pool.scores.index(max(nucl_pool.scores)))

        else:
            for nucl in nucl_pool.nucls:
                for i in range(0, design_parameters.num_mutants):
                    nucl_pool.append(mutate(nucl=nucl, design_parameters=design_parameters))
                
                for i in range(0, design_parameters.num_mutants):
                    nucl_pool.remove(nucl_pool.scores.index(max(nucl_pool.scores)))

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
            nucl_pool.save("RESULTS/" + "MID_" +time.asctime() + "_w" + str(min(design_parameters.weights[0:6])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fastq')

        prev_min = current_min # type: ignore
    

    nucl_pool.save("RESULTS/" + "END_" +time.asctime() + "_w" + str(min(design_parameters.weights[0:7])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fastq')

    if parallel_pool != None:
        parallel_pool.close()
        parallel_pool.terminate()
