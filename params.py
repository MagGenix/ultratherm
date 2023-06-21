from blist import blacklist

class design_parameters():
    def __init__(self, blacklist: blacklist, target: int, offset: int, temp_factor: int, weight_factor: int, num_mutants: int, program: str, weights:list):
        self.blacklist = blacklist
        self.target = target
        self.offset = offset
        self.temp_factor = temp_factor
        self.weight_factor = weight_factor
        self.program = program
        self.num_mutants = num_mutants

        if len(weights) != 7:
            raise TypeError
        for weight in weights:
            if type(weights) != int:
                weight = int(weight)
            if weight > 16 or weight < 0:
                raise TypeError
        self.weights = weights
    
    def can_decrement_offset(self):
        return self.offset > self.temp_factor

    def decrement_offset(self):
        if self.offset > self.temp_factor:
            self.offset -= self.temp_factor
        else:
            raise ValueError
    
    def can_decrement_weights(self):
        return min(self.weights[0:7]) > self.weight_factor

    def decrement_weights(self):
        if min(self.weights[0:7]) > self.weight_factor:
            for i in range(0, 7):
                self.weights[i] -= self.weight_factor
        else:
            raise ValueError