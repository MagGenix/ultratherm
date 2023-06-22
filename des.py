from nucl import nucl_set, mutate
from params import design_parameters
import time

def design(design_parameters:design_parameters, max_reps:int, current_rep:int, pool:nucl_set, prev_min:float, iter_count:int):
    if current_rep == max_reps:
        if design_parameters.can_decrement_weights() or design_parameters.can_decrement_offset():
            if design_parameters.can_decrement_weights():
                design_parameters.decrement_weights()
            if design_parameters.can_decrement_offset():
                design_parameters.decrement_offset()
            current_rep = 0
        else:
            pool.save("END_" +time.asctime() + "_" + str(iter_count) + '.fasta')
            return
    for nucl in pool.nucls:
        for i in range(0, design_parameters.num_mutants):
            pool.append(mutate(nucl=nucl, design_parameters=design_parameters))

        if current_rep == 0:
            nucl.fitness_score(design_parameters=design_parameters)
        
        for i in range(0, design_parameters.num_mutants):
            pool.remove(pool.scores.index(max(pool.scores)))

    current_min = min(pool.scores)
    
    #pool.save(time.asctime() + "_" + str(iter_count) + '.fasta')
    if current_min >= prev_min:
        current_rep+=1
    else:
        current_rep = 0
        if design_parameters.can_decrement_weights():
            design_parameters.decrement_weights()
        if design_parameters.can_decrement_offset():
            design_parameters.decrement_offset()

    iter_count+=1
    print("iter_count:\t" + str(iter_count))
    print("current_rep:\t" + str(current_rep))
    print("offset:\t" + str(design_parameters.offset))
    print("min weight:\t" + str(min(design_parameters.weights[0:7])))
    print('')
    design(design_parameters=design_parameters, max_reps=max_reps, current_rep=current_rep, pool=pool, prev_min=current_min, iter_count=iter_count)