from blist import blacklist
import yaml

#NOTE design_parameters() is the ONLY function that should have default values!!
#Adding default values to other functions may conflict with the ones specified here.
#Any and all arguments to other functions should be required, and if unspecified, that function should error.
class design_parameters():
    """_summary_
    """
    def __init__(
            self, target_energy:float, target_temp: float, temp_offset: float = 5.0, thermo_score_temp:int=37,
            nucl_concentration:float = 1e-6, dimer_max_order_magnitude:float = 2.0,
            blacklist: blacklist = blacklist(''), num_mutants: int = 16, program: str = 'VIENNA',
            weights:list = [16, 16, 16, 16, 16, 16, 16], weight_factor: int = 1,
            free_energy_max_score:float=1.0, nucl_max_score:float=1.0, max_dimer_monomer_factor:float=1.0,
        ):
        
        self.target_energy = target_energy
        self.target_temp = target_temp
        self.temp_offset = temp_offset
        self.thermo_score_temp = thermo_score_temp

        self.nucl_concentration = nucl_concentration
        self.dimer_max_order_magnitude = dimer_max_order_magnitude

        self.blacklist = blacklist
        self.num_mutants = num_mutants
        self.program = program
        
        self.weight_factor = weight_factor

        self.free_energy_max_score = free_energy_max_score
        self.nucl_max_score = nucl_max_score
        self.max_dimer_monomer_factor = max_dimer_monomer_factor

        if free_energy_max_score < 0 or free_energy_max_score > 1:
            raise ValueError
        if nucl_max_score < 0 or nucl_max_score > 1:
            raise ValueError
        if max_dimer_monomer_factor < 0 or max_dimer_monomer_factor > 1:
            raise ValueError

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
    
    #Only decrement the weights for mutations, not the no-mutation weight
    def can_decrement_weights(self) -> bool:
        """_summary_

        Returns:
            bool: _description_
        """
        return min(self.weights[0:6]) > self.weight_factor and self.weight_factor != 0 # type: ignore

    def decrement_weights(self) -> None:
        """_summary_

        Raises:
            ValueError: _description_
        """
        if min(self.weights[0:6]) > self.weight_factor: # type: ignore
            for i in range(0, 6):
                self.weights[i] -= self.weight_factor
        else:
            raise ValueError
        
    def save(self, path:str) -> None:
        """_summary_

        Args:
            path (str): _description_
        """
        yml_dict = {
            'target_energy':            self.target_energy,
            'target_temp':              self.target_temp,
            'temp_offset':              self.temp_offset,
            'thermo_score_temp':        self.thermo_score_temp,

            'nucl_concentration':       self.nucl_concentration,
            'dimer_max_order_magnitude':self.dimer_max_order_magnitude,

            'blacklist':                self.blacklist.blacklist_path,
            'num_mutants':              self.num_mutants,
            'program':                  self.program,

            'weights':                  self.weights,
            'weight_factor':            self.weight_factor,

            'free_energy_max_score':    self.free_energy_max_score,
            'nucl_max_score':           self.nucl_max_score,
            'max_dimer_monomer_factor': self.max_dimer_monomer_factor
        }

        stream = open(path, 'w')
        yaml.safe_dump(data=yml_dict, stream=stream, sort_keys=False)
        stream.close()

def read_parameters(path:str) -> design_parameters:
    """_summary_

    Args:
        path (str): _description_

    Returns:
        design_parameters: _description_
    """
    stream = open(path, "r")
    yml = yaml.safe_load(stream=stream)
    stream.close()

    yml['blacklist'] = blacklist(yml['blacklist'])

    return design_parameters(**yml)
