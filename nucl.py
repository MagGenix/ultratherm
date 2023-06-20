from nupack import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from math import log10
#from ViennaRNA import RNA

from blist import blacklist

class scoring_parameters():
    def __init__(self, blacklist: blacklist, target: float, offset: float, program: str):
        self.blacklist = blacklist
        self.target = target
        self.offset = offset
        self.program = program

class nucl_acid():
    def __init__(self, sequence: Seq, no_mod: list, no_indel: list, score_region: list, scoring_parameters: scoring_parameters):
        if len(no_mod) != len(sequence):
            raise Exception("no_mod length is not equal to seq length")
        if len(no_indel) != len(sequence):
            raise Exception("no_indel length is not equal to seq length")
        if len(score_region) != len(sequence):
            raise Exception("score_region length is not equal to seq length")
        self.sequence = sequence
        self.no_mod = no_mod
        self.no_indel = no_indel
        self.score_region = score_region
        self.score_entire_seq = score_region.count(1) == 0
        self.score = self.fitness_score(scoring_parameters)
    
    def __len__(self):
        return len(self.sequence)

    def fitness_score(self, scoring_parameters: scoring_parameters):
        if scoring_parameters.blacklist.is_blacklisted(self):
            return 4
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
            scores_cold = self.nupack_score_temp(temp=cold_temp, tube_nucl=tube_nucl, complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot=False)

            #Calculate hot temp for scoring
            hot_temp = scoring_parameters.target+scoring_parameters.offset
            if hot_temp < 0:
                raise Exception("illegal hot temperature")
            if hot_temp > 100:
                hot_temp = 100
            scores_hot = self.nupack_score_temp(temp=hot_temp, tube_nucl=tube_nucl, complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot = True)

            return sum(scores_cold) + sum(scores_hot)
        if scoring_parameters.program == "VIENNA":
            return 4
        raise Exception("no program specified for scoring")

    def nupack_score_temp(self, temp: float, tube_nucl: Tube, complex_nucl_single: Complex, complex_nucl_double: Complex, hot:bool):
        #Make NUPACK model
        model_nucl=Model(celsius=temp)

        #Run NUPACK model for temp (Complex analysis)
        results_nucl = complex_analysis(complexes = tube_nucl, model=model_nucl, compute=['pfunc', 'pairs'])
        concentrations_nucl = complex_concentrations(tube=tube_nucl, data = results_nucl)

        #Calculate ratio of AA to A and take log10. lower is better
        if concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_double] == 0:
            dimer_monomer_factor=0
        elif concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_single] == 0:
            dimer_monomer_factor=1
        else:
            dimer_monomer_factor = log10(concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_double] /
                concentrations_nucl.tubes[tube_nucl].complex_concentrations[complex_nucl_single]) + 2 #+2 means dimer must be 2 factors of 10 less abundant to avoid score penalty
            if dimer_monomer_factor < 0:
                dimer_monomer_factor = 0 #0 is the best possible factor, indicates limited dimer formation
            if dimer_monomer_factor > 1:
                dimer_monomer_factor = 1 #cap cost of having a poor monomer formation
        
        if self.score_entire_seq:
            score_nucl = results_nucl.fraction_bases_unpaired
        else:
            score_nucl = 0
            count_scored_nuc = 0
            for i, x in enumerate(self.score_region):
                if x==1:
                    score_nucl += results_nucl.complexes[complex_nucl_single].pairs.diagonal[i]
                    count_scored_nuc+=1
            score_nucl = score_nucl / count_scored_nuc

        if hot:
            score_nucl = 1-score_nucl
        return [dimer_monomer_factor, score_nucl]

class nucl_set():
    def __init__(self, nucls: list):
        self.nucls = nucls
        self.scores = []
        for nucl in self.nucls:
            if type(nucl) != nucl_acid:
                #Raise an error if an object in the list was not a nucl_acid
                raise TypeError
            self.scores = nucl.score
    def __len__(self):
        return len(self.nucls)

    def replace(self, index:int, new_nucl_acid:nucl_acid):
        if index < 0:
            raise IndexError
        
        #nucl_set is zero-indexed!
        if index > len(self) - 1:
            raise IndexError
        self.nucls[index] = new_nucl_acid
        self.scores[index] = new_nucl_acid.score

    def append(self, new_nucl_acid:nucl_acid):
        self.nucls.append(new_nucl_acid)
        self.scores.append(new_nucl_acid.score)

    def remove(self, index:int):
        del self.nucls[index]
        del self.scores[index]

    def save(self, path:str):
        #raises IOError if path is invalid
        with open(path) as handle:
            for i in range(0, len(self)):
                SeqIO.write(SeqRecord(seq=self.nucls(i), name=str(i),description="score=" + self.nucls(i).score), handle=handle, format='fasta')