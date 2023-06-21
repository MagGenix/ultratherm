from nupack import *
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from math import log10
import random
import copy
#from ViennaRNA import RNA

from params import design_parameters

class nucl_acid():
    def __init__(self, sequence: Seq, no_mod: list, no_indel: list, score_region: list, design_parameters: design_parameters, is_rna: bool):
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
        self.is_rna = is_rna
        if self.sequence.find('N') != -1:
            self.sequence = MutableSeq(self.sequence)
            if self.is_rna:
                nucleotides = ['A', 'U', 'G', 'C']
            else:
                nucleotides = ['A', 'T', 'G', 'C']
            while self.sequence.find('N') != -1:
                self.sequence[self.sequence.find('N')] = nucleotides[random.getrandbits(2)]
            self.sequence = Seq(self.sequence)
        #Score the entire thing if everything was set to 0. Have to score something at minimum
        if score_region.count(0) == len(score_region):
            self.score_region = [1] * len(score_region)

        self.score = self.fitness_score(design_parameters)

    
    def __len__(self):
        return len(self.sequence)

    def fitness_score(self, design_parameters: design_parameters):
        if design_parameters.blacklist.is_blacklisted(self):
            return 4
        if design_parameters.program == "NUPACK":
            #TODO: specify whether sequence to be analyzed is DNA or RNA
            #Create NUPACK strand using sequence of nucl
            strand_nucl = Strand(name='A', string=str(self.sequence))

            #Create NUPACK complexes for monomer and homodimer
            complex_nucl_single = Complex(strands=[strand_nucl], name='A',)
            complex_nucl_double = Complex(strands=[strand_nucl, strand_nucl], name='AA')

            #Create NUPACK Tube and track both the monomer and homodimer complexes
            tube_nucl = Tube(strands={strand_nucl:1e-6}, complexes=SetSpec(max_size=2,
               include=[complex_nucl_single, complex_nucl_double]), name='tube_nucl')
            
            #Calculate cold temp for scoring
            cold_temp = design_parameters.target-design_parameters.offset
            if cold_temp < 0:
                cold_temp = 0
            if cold_temp >= 100:
                raise Exception("illegal cold temperature")
            scores_cold = self.nupack_score_temp(temp=cold_temp, tube_nucl=tube_nucl, complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot=False)

            #Calculate hot temp for scoring
            hot_temp = design_parameters.target+design_parameters.offset
            if hot_temp < 0:
                raise Exception("illegal hot temperature")
            if hot_temp > 100:
                hot_temp = 100
            scores_hot = self.nupack_score_temp(temp=hot_temp, tube_nucl=tube_nucl, complex_nucl_single=complex_nucl_single, complex_nucl_double=complex_nucl_double, hot = True)

            return sum(scores_cold) + sum(scores_hot)
        if design_parameters.program == "VIENNA":
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
        

        score_nucl = 0
        count_scored_nuc = 0
        for i, x in enumerate(self.score_region):
            if x:
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
        #raises IOError if path is invalid, overwrites existing files
        with open(path, 'w') as handle:
            for i in range(0, len(self)):
                SeqIO.write(SeqRecord(seq=self.nucls[i].sequence, id=str(i), description="score=" + str(self.nucls[i].score)), handle=handle, format='fasta')

def mutate(nucl:nucl_acid, design_parameters:design_parameters):
    #I'm trying to make this function as fast as possible since it will be called once per every single
    #nucleotide in a sequence to generate one variant.
    #For larger pool sizes this can really stack up.
    
    if len(nucl.no_mod) != len(nucl.sequence) or len(nucl.no_indel) != len(nucl.score_region) or len(nucl.no_mod) != len(nucl.no_indel):
        raise ValueError

    #As of writing the design parameters function forces the weights to be whole numbers between 0 and 16.
    #This could be made a touch faster explicitly using ints.
    weights = copy.copy(design_parameters.weights)
    weights_total = sum(weights)

    #Make the weights array into an array of increasing values which can be compared to a random int
    for i in range(0, len(weights)):
        weights[i] = weights[i - 1] + weights[i]

    sequence = list(copy.copy(nucl.sequence))
    no_mod = copy.copy(nucl.no_mod)
    no_indel = copy.copy(nucl.no_indel)
    score_region = copy.copy(nucl.score_region)

    #Define nucleotide set
    #Really don't want to be repeatedly checking if RNA or DNA inside loop
    if nucl.is_rna:
        nucleotides = ['A', 'U', 'G', 'C']
    else:
        nucleotides = ['A', 'T', 'G', 'C']

    #mod_range = list(range(0, len(sequence)))

    i = 0
    while i < len(sequence):
        if no_mod[i]:
            i = i+1
            continue
        choice = random.randint(0, weights_total)

        #Comparing the generated int with values in the array of increasing ints.
        for j in range(0, len(weights)):
            if choice > weights[j]:
                continue
            selection = j
            break
        
        #If no change selected, move on
        if selection == 6:
            i = i+1
            continue
        
        #If indels not allowed for this position, continue
        if (selection == 4 or selection == 5) and no_indel[i]:
            i = i+1
            continue
        
        #Index the nucleotide array instead of checking several times. Hope it's faster
        #It is fewer lines, anyway
        if selection >= 0 and selection < 4:
            sequence[i] = nucleotides[selection]
            i = i+1
            continue

        # Deletion is faster than insertion so I've moved it here for now
        if selection == 5:
            del sequence[i]

            #The other arrays need to be the same size!
            del no_mod[i]
            del no_indel[i]
            del score_region[i]
            #If this was the last nucleotide, the while loop will see that anyway
            #Repeat this nucleotide for the next mod so no incrementing i
            continue

        # If an insertion was selected
        #Unnecessary check (selection == 5) but I'll leave it here in case the order changes
        #TODO: remove unnecessary check
        if selection == 4:
            #Generate a random nucleotide from nucleotide list (getrandbits (2) will generate an int 0-4)
            new_nt = nucleotides[random.getrandbits(2)]
            

            #I'm pretty worried about the score region changing code. It could be slow and could have bug(s)
            
            #Randomly generate bit to determine whether to insert to the left or right
            #If bit was 1, insert to the left
            if random.getrandbits(1):
                sequence.insert(i, new_nt)
                no_mod.insert(i, 0)
                no_indel.insert(i, 0)


                #If at position 0, there is nothing to the left!
                if i == 0 and score_region[0]:
                    score_region.insert(i, 1)
                    i = i+2
                    continue
                #At any other position besides 0, check to the left and right
                elif score_region[i - 1] and score_region[i]:
                    score_region.insert(i, 1)
                    i = i+2
                    continue
                
                #default behavior is to insert a 0
                score_region.insert(i, 0)
                i = i+2
                continue
            #Else, insert to the right

            #Check if we are at the end - use append instead
            if i == len(sequence) - 1:
                sequence.append(new_nt)
                no_mod.append(0)
                no_indel.append(0)

                #If we are at the end and it's scored, extend the score region
                if score_region[i]:
                    score_region.append(1)
                
                #Exit the loop because that's it
                break
            
            #Otherwise, we insert
            sequence.insert(i+1, new_nt)
            no_mod.insert(i+1, 0)
            no_indel.insert(i+1, 0)
            
            if score_region[i] and score_region[i+1]:
                score_region.insert(i+1, 1)
                i = i+2
                continue
            
            #default behavior is to insert a 0
            score_region.insert(i+1, 0)
            i = i+2
            #Unnecessary but I'll leave it here in case the order changes
            #TODO: remove unnecessary continue
            continue
    
    return nucl_acid(sequence=Seq(''.join(sequence)), no_mod=no_mod, no_indel=no_indel, score_region=score_region,
                     design_parameters=design_parameters, is_rna=nucl.is_rna)