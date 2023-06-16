from nupack import *
from Bio.Seq import Seq
#from ViennaRNA import RNA

from blist import blacklist

class scoring_parameters():
    def __init__(self, blacklist: blacklist, target: float, offset: float, program: str):
        self.blacklist = blacklist
        self.target = target
        self.offset = offset
        self.program = program

def fitness_score(self, scoring_parameters: scoring_parameters):
        if scoring_parameters.blacklist.is_blacklisted(self):
            return 2
        if scoring_parameters.program == "NUPACK":
            #Create NUPACK strand using sequence of nucl
            strand_nucl = Strand(name='A', string=str(self.sequence))

            #Create NUPACK complexes for monomer and homodimer
            complex_nucl_single = Complex(strands=[strand_nucl], name='A',)
            complex_nucl_double = Complex(strands=[strand_nucl, strand_nucl], name='AA')

            #Create NUPACK Tube and track both the monomer and homodimer complexes
            tube_nucl = Tube(strands={strand_nucl:1e-6}, complexes=SetSpec(max_size=2,
                include=[complex_nucl_single, complex_nucl_double]), name='tube_nucl')
            
            #Calculate cold temp for scoring
            cold_temp = scoring_parameters.target-scoring_parameters.offset
            if cold_temp < 0:
                cold_temp = 0
            if cold_temp >= 100:
                raise Exception("illegal cold temperature")
            
            #Make NUPACK model for cold temp
            model_nucl_cold=Model(celsius=cold_temp)

            #Run NUPACK model for cold temp (Tube analysis)
            results_nucl_cold = tube_analysis(tubes=[tube_nucl], model=model_nucl_cold, compute=['pairs'])
            nucl_cold_complex_concentrations = results_nucl_cold.tubes[tube_nucl].complex_concentrations
            
            #Calculate ratio of AA to A
            dimer_monomer_ratio_cold = nucl_cold_complex_concentrations[complex_nucl_double] / nucl_cold_complex_concentrations[complex_nucl_single]

            print(dimer_monomer_ratio_cold)
            return 2
        if scoring_parameters.program == "VIENNA":
            return 2
        raise Exception("no program specified for scoring")

class nucl_acid():
    def __init__(self, sequence: Seq, no_mod: list, no_indel: list, scoring_parameters: scoring_parameters):
        if len(no_mod) != len(sequence):
            raise Exception("no_mod length is not equal to seq length")
        if len(no_indel) != len(sequence):
            raise Exception("no_indel length is not equal to seq length")
        self.sequence = sequence
        self.no_mod = no_mod
        self.no_indel = no_indel
        self.score = fitness_score(self, scoring_parameters)
        self.score_entire_seq = no_mod.count(1) == 0