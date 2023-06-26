from blist import blacklist

class design_parameters():
    def __init__(self, blacklist: blacklist, target_temp: int, temp_offset: int, weight_factor: int, num_mutants: int, program: str, weights:list, target_energy:float=-10.0):
        self.blacklist = blacklist
        self.target_temp = target_temp
        self.temp_offset = temp_offset
        self.weight_factor = weight_factor
        self.program = program
        self.num_mutants = num_mutants
        self.target_energy = target_energy

        if target_energy >= 0:
            raise ValueError

        if len(weights) != 7:
            raise TypeError
        for weight in weights:
            if type(weights) != int:
                weight = int(weight)
            if weight > 16 or weight < 0:
                raise TypeError
        self.weights = weights
    
    def can_decrement_weights(self):
        return min(self.weights[0:7]) > self.weight_factor and self.weight_factor != 0

    def decrement_weights(self):
        if min(self.weights[0:7]) > self.weight_factor:
            for i in range(0, 7):
                self.weights[i] -= self.weight_factor
        else:
            raise ValueError
        
    def save(self, path:str):
        try:
            with open(path, 'w') as handle:
                handle.write('blacklist: ' + self.blacklist.blacklist_path + '\n')
                handle.write('target_temp: ' + str(self.target_temp) + '\n')
                handle.write('temp_offset: ' + str(self.temp_offset) + '\n')
                handle.write('weight_factor: ' + str(self.weight_factor) + '\n')
                handle.write('program: ' + self.program + '\n')
                handle.write('num_mutants: ' + str(self.num_mutants) + '\n')
                handle.write('weights:\n')
                for weight in self.weights:
                    handle.write(' - ' + str(weight)+"\n")
        except IOError:
            print("Warning: could not save parameters")
