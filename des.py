from nucl import nucl_set, mutate
from params import design_parameters

def design(design_parameters:design_parameters, max_reps:int, current_rep:int, pool:nucl_set, prev_min:float, iter_count:int):
    if current_rep == max_reps:
        return
    if not design_parameters.can_decrement():
        return
    for nucl in pool.nucls:
        pool.append(mutate(nucl=nucl, design_parameters=design_parameters))
        pool.remove(pool.scores.index(max(pool.scores)))
    current_min = min(pool.scores)
    
    pool.save(str(iter_count) + '.fasta')
    if current_min >= prev_min:
        current_rep+=1
    else:
        current_rep = 0
        design_parameters.decrement_offset()

    iter_count+=1
    print(iter_count)
    print(current_rep)
    print(design_parameters.offset)
    design(design_parameters=design_parameters, max_reps=max_reps, current_rep=current_rep, pool=pool, prev_min=current_min, iter_count=iter_count)