from blist import blacklist
from nucl import nucl_set

class design_parameters():
    def __init__(self, pool_size:int, blacklist: blacklist, target: float, offset: float, program: str, weights:list):
        self.pool_size = int
        self.blacklist = blacklist
        self.target = target
        self.offset = offset
        self.program = program
        if len(weights) != 7:
            raise TypeError
        for weight in weights:
            if type(weights) != int:
                weight = int(weight)
            if weight > 16 or weight < 0:
                raise TypeError
        self.weights = weights
    
    def decrement_offset(self, factor:float):
        if self.offset > factor:
            self.offset -= factor
        else:
            raise ValueError
        
def design(design_parameters: design_parameters, max_cycles: int, current_cycle:int, factor:float, old_nucl_set: nucl_set):
    pass
