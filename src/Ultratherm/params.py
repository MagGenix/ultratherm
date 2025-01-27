import os
from blist import blacklist
import yaml

#NOTE design_parameters() is the ONLY function that should have default values!!
#Adding default values to other functions may conflict with the ones specified here.
#Any and all arguments to other functions should be required, and if unspecified, that function should error.
class design_parameters():
    """Stores parameters used in nucleic acid design.
    """
    def __init__(
            self, target_energy:float, target_temp: float, temp_offset: float = 5.0, thermo_score_temp:float=37, max_hairpins:int=1,
            parasitic_max_order_magnitude:float = 2.0,
            blacklist: blacklist = blacklist(''), num_mutants: int = 16, program: str = 'VIENNA', parallel:bool = True,
            weights:list = [16, 16, 16, 16, 16, 16, 16], weight_factor: int = 1, max_reps:int = 32,
            free_energy_max_score:float=1.0, accessibility_max_score:float=1.0, parasitic_complex_max_score:float=1.0, num_hairpins_max_score:float=1.0,
            optimization_rate:float = 5.0, result_save_path:str = 'RESULTS'
        ):
        """Create a new design_parameters object. THIS IS THE ONLY FUNCTION WITH DEFAULTS!

        Args:
            target_energy (float): the target ensemble free energy (kcal/mol). If higher than target, penalties are incurred (higher score). If lower, no penalty.
            target_temp (float): Target temp for 50% pairing of the score_region. At target_temp - temp_offset, pairing should be 0% and 100% at target_temp + temp_offset.
            temp_offset (float, optional): The offset temperature for score_region pair probability assessment. See target_temp. Defaults to 5.0.
            thermo_score_temp (int, optional): The temperature at which the ensemble energy is scored. Defaults to 37.
            max_hairpins (int, optional): the maximum number of acceptable hairpins in the ensemble centroid structure.
            parasitic_max_order_magnitude (float, optional): The threshold at which to penalize dimer formation, as -log10([DIMER] / [MONOMER]). Defaults to 2.0.
            blacklist (blacklist, optional): A blacklist object. Defaults to blacklist('').
            num_mutants (int, optional): The number of mutants to generate per nucl_acid in the nucl_set. Defaults to 16.
            program (str, optional): 'NUPACK' or 'VIENNA'. Defaults to 'VIENNA'.
            parallel (bool, optional): Whether or not to parallelize mutation by computing mutants with different workers. Only applies to VIENNA.
            weights (list, optional): Mutation weights. All weights besides the no modification weight are decremented in the design loop. The higher the weight, the higher the probabiltiy said mutation is chosen. Each are above 0. [A, T/U, G, C, insertion, deletion, no modification]. Defaults to [16, 16, 16, 16, 16, 16, 16].
            weight_factor (int, optional): How much to decrement the weights by. Defaults to 1.
            max_reps (int, optional): the maximum number of loops to perform without a decrease in minimum pool score.
            free_energy_max_score (float, optional): The maximum score penalty for having a free energy greater than target. Defaults to 1.0.
            accessibility_max_score (float, optional): The maximum score penalty for score region accessibility. Defaults to 1.0.
            parasitic_complex_max_score (float, optional): The maximum score penalty for dimer formation. Defaults to 1.0.
            num_hairpins_max_score (float, optional): the maximum score penalty for having more hairpins in the ensemble centroid than the desired max.
            optimization_rate (float, optional): A representation of how quickly to remove suboptimal structures from the pool, 1 meaning slowest and >10 meaning very fast. Minimum 0.0. Maximum 14.0. Defaults to 5.0.
            results_save_path (str, optional): The path to save results. Defaults to 'RESULTS/'.

        Raises:
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            TypeError: _description_
            TypeError: _description_
        """

        try:
            result_save_path = os.path.normpath(result_save_path)
            if not os.path.exists(result_save_path):
                os.makedirs(result_save_path)
            result_save_path += os.sep
        except OSError:
            print("results_save_path not writable!")
            raise ValueError
        
        self.result_save_path = result_save_path

        self.target_energy = target_energy
        self.target_temp = target_temp
        self.temp_offset = temp_offset
        self.thermo_score_temp = thermo_score_temp
        self.max_hairpins = max_hairpins

        self.parasitic_max_order_magnitude = parasitic_max_order_magnitude

        self.blacklist = blacklist
        self.num_mutants = num_mutants
        self.program = program
        self.parallel = parallel
        
        self.weight_factor = weight_factor
        self.max_reps = max_reps

        self.free_energy_max_score          = free_energy_max_score
        self.accessibility_max_score        = accessibility_max_score
        self.parasitic_complex_max_score    = parasitic_complex_max_score
        self.num_hairpins_max_score         = num_hairpins_max_score
        
        self.optimization_rate              = optimization_rate

        if max_hairpins < 0:
            raise ValueError

        if free_energy_max_score < 0:
            raise ValueError
        if accessibility_max_score < 0:
            raise ValueError
        if parasitic_complex_max_score < 0:
            raise ValueError
        if num_hairpins_max_score < 0:
            raise ValueError
        
        if optimization_rate <= 0.0 or optimization_rate >=14: # 14.78 will hit float32 limit
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
            if weight < 0:
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
            'target_energy':                self.target_energy,
            'target_temp':                  self.target_temp,
            'temp_offset':                  self.temp_offset,
            'thermo_score_temp':            self.thermo_score_temp,
            'max_hairpins':                 self.max_hairpins,

            'parasitic_max_order_magnitude':    self.parasitic_max_order_magnitude,

            'blacklist':                    self.blacklist.blacklist_path,
            'num_mutants':                  self.num_mutants,
            'program':                      self.program,
            'parallel':                     self.parallel,

            'weights':                      self.weights,
            'weight_factor':                self.weight_factor,
            'max_reps':                     self.max_reps,

            'free_energy_max_score':        self.free_energy_max_score,
            'accessibility_max_score':      self.accessibility_max_score,
            'parasitic_complex_max_score':  self.parasitic_complex_max_score,
            'num_hairpins_max_score':       self.num_hairpins_max_score,
            
            'optimization_rate':            self.optimization_rate,
            'result_save_path':             self.result_save_path
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
