from blist import blacklist

class design_parameters():
    def __init__(self, blacklist: blacklist, target: float, offset: float, program: str, weights:list):
        self.blacklist = blacklist
        self.target = target
        self.offset = offset
        self.program = program
        if len(weights) != 8:
            raise TypeError
        for weight in weights:
            if type(weights) != int:
                weight = int(weight)
            if weight > 16 or weight < 0:
                raise TypeError