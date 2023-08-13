from nucl import nucl_set, mutate
from params import design_parameters
import time

def design(design_parameters:design_parameters, max_reps:int, current_rep:int, pool:nucl_set, prev_min:float, iter_count:int):
    #Note - now that the temp offset decrementing code is gone, the nucl_acid's are NEVER rescored (the temp never changes). This is less expensive.
    #If for some reason that becomes necessary in the future (I doubt it), it will have to be added back.
    if current_rep == max_reps:
        if design_parameters.can_decrement_weights():
            design_parameters.decrement_weights()
            current_rep = 0
        else:
            pool.save("END_" +time.asctime() + "_w" + str(min(design_parameters.weights[0:7])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fasta')
            return
    for nucl in pool.nucls:
        for i in range(0, design_parameters.num_mutants):
            #TODO consider changing this to bracket-indexed replacement
            #Probably numpy array
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
    if iter_count % 50 == 0:
        pool.save("MID_" +time.asctime() + "_w" + str(min(design_parameters.weights[0:7])) + "_o" + str(design_parameters.temp_offset) + "_i" + str(iter_count) + '.fasta')
    
    print("iter_count:\t" + str(iter_count))
    print("current_rep:\t" + str(current_rep))
    print("min weight:\t" + str(min(design_parameters.weights[0:6])))
    print('')
    design(design_parameters=design_parameters, max_reps=max_reps, current_rep=current_rep, pool=pool, prev_min=current_min, iter_count=iter_count) # type: ignore