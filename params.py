from blist import blacklist
import yaml

#NOTE design_parameters() is the ONLY function that should have default values!!
#Adding default values to other functions may conflict with the ones specified here.
#Any and all arguments to other functions should be required, and if unspecified, that function should error.
class design_parameters():
    """Stores parameters used in nucleic acid design.
    """
    def __init__(
            self, target_energy:float, target_temp: float, temp_offset: float = 5.0, thermo_score_temp:int=37,
            nucl_concentration:float = 1e-6, dimer_max_order_magnitude:float = 2.0,
            blacklist: blacklist = blacklist(''), num_mutants: int = 16, program: str = 'VIENNA', parallel:bool = True,
            weights:list = [16, 16, 16, 16, 16, 16, 16], weight_factor: int = 1, max_reps:int = 16,
            free_energy_max_score:float=1.0, nucl_max_score:float=1.0, max_dimer_monomer_factor:float=1.0,
        ):
        """Create a new design_parameters object. THIS IS THE ONLY FUNCTION WITH DEFAULTS!

        Args:
            target_energy (float): the target ensemble free energy (kcal/mol). If higher than target, penalties are incurred (higher score). If lower, no penalty.
            target_temp (float): Target temp for 50% pairing of the score_region. At target_temp - temp_offset, pairing should be 0% and 100% at target_temp + temp_offset.
            temp_offset (float, optional): The offset temperature for score_region pair probability assessment. See target_temp. Defaults to 5.0.
            thermo_score_temp (int, optional): The temperature at which the ensemble energy is scored. Defaults to 37.
            nucl_concentration (float, optional): The concentration of the nucleic acids. Defaults to 1e-6.
            dimer_max_order_magnitude (float, optional): The threshold at which to penalize dimer formation, as -log10([DIMER] / [MONOMER]). Defaults to 2.0.
            blacklist (blacklist, optional): A blacklist object. Defaults to blacklist('').
            num_mutants (int, optional): The number of mutants to generate per nucl_acid in the nucl_set. Defaults to 16.
            program (str, optional): 'NUPACK' or 'VIENNA'. Defaults to 'VIENNA'.
            parallel (bool, optional): Whether or not to parallelize mutation by computing mutants with different workers. Only applies to VIENNA.
            weights (list, optional): Mutation weights. All weights besides the no modification weight are decremented in the design loop. The higher the weight, the higher the probabiltiy said mutation is chosen. Each vary from 0 to 16. [A, T/U, G, C, insertion, deletion, no modification]. Defaults to [16, 16, 16, 16, 16, 16, 16].
            weight_factor (int, optional): How much to decrement the weights by. Defaults to 1.
            max_reps (int, optional): the maximum number of loops to perform without a decrease in minimum pool score.
            free_energy_max_score (float, optional): The maximum score penalty for having a free energy greater than target. Defaults to 1.0.
            nucl_max_score (float, optional): The maximum score penalty for score region accessibility. Defaults to 1.0.
            max_dimer_monomer_factor (float, optional): The maximum score penalty for dimer formation. Defaults to 1.0.

        Raises:
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            TypeError: _description_
            TypeError: _description_
        """
        
        self.target_energy = target_energy
        self.target_temp = target_temp
        self.temp_offset = temp_offset
        self.thermo_score_temp = thermo_score_temp

        self.nucl_concentration = nucl_concentration
        self.dimer_max_order_magnitude = dimer_max_order_magnitude

        self.blacklist = blacklist
        self.num_mutants = num_mutants
        self.program = program
        self.parallel = parallel
        
        self.weight_factor = weight_factor
        self.max_reps = max_reps

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

        if max_reps < 1:
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
        """Returns whether the mutation weights can be decremented by the weight_factor.

        Returns:
            bool: mutation weights can be decremented (True / False)
        """
        return min(self.weights[0:6]) > self.weight_factor and self.weight_factor != 0 # type: ignore

    def decrement_weights(self) -> None:
        """Decrements the first 6 weights (not the no-mutation weight) by the weight_factor.

        Raises:
            ValueError: Attempted to decrement weights where impossible.
        """
        if min(self.weights[0:6]) > self.weight_factor: # type: ignore
            for i in range(0, 6):
                self.weights[i] -= self.weight_factor
        else:
            raise ValueError
        
    def save(self, path:str) -> None:
        """Writes the design parameters to a .yml file.

        Args:
            path (str): Path to save the file. If no directory specified will save to source directory.
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
            'parallel':                 self.parallel,

            'weights':                  self.weights,
            'weight_factor':            self.weight_factor,
            'max_reps':                 self.max_reps,

            'free_energy_max_score':    self.free_energy_max_score,
            'nucl_max_score':           self.nucl_max_score,
            'max_dimer_monomer_factor': self.max_dimer_monomer_factor
        }

        stream = open(path, 'w')
        yaml.safe_dump(data=yml_dict, stream=stream, sort_keys=False)
        stream.close()

def read_parameters(path:str) -> design_parameters:
    """Reads design parameters from a .yml file and generates a design_parameters object.
    If optargs are omitted, defaults will be used.

    Args:
        path (str): Path to read the file. If no directory specified will assume source directory.

    Returns:
        design_parameters: the design parameters read from the file.
    """
    stream = open(path, "r")
    yml = yaml.safe_load(stream=stream)
    stream.close()

    yml['blacklist'] = blacklist(yml['blacklist'])

    return design_parameters(**yml)
