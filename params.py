from blist import blacklist

class design_parameters():
    def __init__(self, blacklist: blacklist, target: float, offset: float, factor: float, program: str, weights:list):
        self.blacklist = blacklist
        self.target = target
        self.offset = offset
        self.factor = factor
        self.program = program
        if len(weights) != 7:
            raise TypeError
        for weight in weights:
            if type(weights) != int:
                weight = int(weight)
            if weight > 16 or weight < 0:
                raise TypeError
        self.weights = weights
    
    def can_decrement(self):
        return self.offset > self.factor

    def decrement_offset(self):
        if self.offset > self.factor:
            self.offset -= self.factor
        else:
            raise ValueError
        
#def design(design_parameters: design_parameters, max_cycles: int, current_cycle:int, factor:float, old_nucl_set: nucl_set):
#    pass
    #Only increase current_cycle if no score improvement happened since the last iteration!
    #Stop when a certain number of cycles with no improvement has been reached, or if the offset cannot be decremented anymore