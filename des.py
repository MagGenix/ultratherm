from nucl import nucl_set, mutate
from params import design_parameters

def design(design_parameters:design_parameters, max_reps:int, current_rep:int, pool:nucl_set, prev_max:float, iter_count:int):
    if current_rep == max_reps:
        return
    if design_parameters.can_decrement():
        design_parameters.decrement_offset()
    else:
        return
    for nucl in pool.nucls:
        pool.append(mutate(nucl=nucl, design_parameters=design_parameters))
        pool.remove(pool.scores.index(min(pool.scores)))
    current_max = max(pool.scores)
    
    pool.save(str(iter_count) + '.fasta')
    if current_max >= prev_max:
        current_rep+=1
    else:
        current_rep = 0

    iter_count+=1
    design(design_parameters=design_parameters, max_reps=max_reps, current_rep=current_rep, pool=pool, prev_max=current_max)