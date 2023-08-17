from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import copy

from params import design_parameters
from blist import blacklist
from nupack_score import nupack_score
from vienna_score import vienna_score

class nucl_acid():
    """_summary_
    """
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

        #Check that definded RNAs do not contain DNA bases and vice versa
        if self.is_rna and self.sequence.find('T') != -1:
            raise ValueError
        elif (not self.is_rna) and self.sequence.find('U') != -1:
            raise ValueError

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

        self.fitness_score(design_parameters)

    
    def __len__(self) -> int:
        return len(self.sequence)
    
    def __str__(self) -> str:
        as_string = ""
        as_string += str(self.score) + "\n"
        as_string += str(self.sequence) + "\n"
        as_string += str(self.no_mod) + "\n"
        as_string += str(self.no_indel) + "\n"
        as_string += str(self.score_region) + "\n"
        return as_string

    def fitness_score(self, design_parameters: design_parameters) -> None:
        """_summary_

        Args:
            design_parameters (design_parameters): _description_

        Raises:
            Exception: _description_
        """
        if self.is_blacklisted(blacklist=design_parameters.blacklist):
            self.score = 6.0
        if design_parameters.program == "NUPACK":  
            self.score = nupack_score(sequence=str(self.sequence), score_region=self.score_region, is_rna=self.is_rna, design_parameters=design_parameters)
        if design_parameters.program == "VIENNA":
            self.score = vienna_score(sequence=str(self.sequence), score_region=self.score_region, is_rna=self.is_rna, design_parameters=design_parameters)
        raise Exception("no program specified for scoring")

    #As of right now, this function is only accessed from within fitness_score, but could be useful to external processes.
    #If used in external processes it would be more efficent to store it as a boolean member variable.
    #However, since member variables are public in Python, this would make is_blacklisted modifiable and potentially inaccurate.
    #For that reason, and since it should only be performed once but is also generally important,
    #this stays a function for now.
    def is_blacklisted(self, blacklist: blacklist) -> bool:
        """_summary_

        Args:
            blacklist (blacklist): _description_

        Returns:
            bool: _description_
        """
        if blacklist.is_empty:
            return False
        if self.sequence in blacklist.blacklist_sequences:
            return True
        return False

class nucl_set():
    """_summary_
    """
    def __init__(self, nucls: list):
        self.nucls = nucls
        self.scores = []
        for i in range (0, len(self.nucls)):
            if type(self.nucls[i]) != nucl_acid:
                #Raise an error if an object in the list was not a nucl_acid
                raise TypeError
            self.scores[i] = self.nucls[i].score
    def __len__(self) -> int:
        return len(self.nucls)

    def __str__(self) -> str:
        as_string = ""
        for nucl in self.nucls:
            as_string = as_string + str(nucl) + "\n"
        return as_string

    def replace(self, index:int, new_nucl_acid:nucl_acid) -> None:
        """_summary_

        Args:
            index (int): _description_
            new_nucl_acid (nucl_acid): _description_

        Raises:
            IndexError: _description_
            IndexError: _description_
        """
        if index < 0:
            raise IndexError
        
        #nucl_set is zero-indexed!
        if index > len(self) - 1:
            raise IndexError
        self.nucls[index] = new_nucl_acid
        self.scores[index] = new_nucl_acid.score

    def append(self, new_nucl_acid:nucl_acid) -> None:
        """_summary_

        Args:
            new_nucl_acid (nucl_acid): _description_
        """
        self.nucls.append(new_nucl_acid)
        self.scores.append(new_nucl_acid.score)

    def remove(self, index:int) -> None:
        """_summary_

        Args:
            index (int): _description_
        """
        del self.nucls[index]
        del self.scores[index]

    def save(self, path:str) -> None:
        """_summary_

        Args:
            path (str): _description_
        """
        #raises IOError if path is invalid, overwrites existing files
        #This saves to a .fastq, where the quality scores are used instead for bitwise indicating nomod, noindel, scoreregion
        #This will be done with ASCII 0-7, p-w  [XXX011Y]
        #These were chosen to avoid any extraneous ASCII characters.
        #b1: no_mod
        #b2: no_indel
        #b3: score_region
        #b7: DNA / RNA
        # AS SANGER SCORES:
        # 0: 15, 7: 22
        # p: 79, w: 86
        with open(path, 'w') as handle:
            for i in range(0, len(self)):
                nucl = self.nucls[i]
                if nucl.is_rna:
                    offset = 79
                else:
                    offset = 15
                quals = list()
                for j in range (0, len(nucl)):
                    bitString = str(nucl.no_mod[j]) + str(nucl.no_indel[j]) + str(nucl.score_region[j])
                    bitsAsInt = int(bitString, 2)
                    quals.append(bitsAsInt + offset)
                record = SeqRecord(seq=self.nucls[i].sequence, id=str(i), description="score=" + str(self.nucls[i].score))
                record.letter_annotations["phred_quality"] = quals
                SeqIO.write(record, handle=handle, format='fastq')
                del quals
                del record
                del nucl

    def read(self, path:str, design_parameters:design_parameters) -> None:
        """_summary_

        Args:
            path (str): _description_
            design_parameters (design_parameters): _description_

        Raises:
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
        """
        #If any other characters are encountered, an error should be raised.
        for record in SeqIO.parse(path, "fastq"):
            quals = record.letter_annotations["phred_quality"]
            if max(quals) > 86 or min(quals) < 15:
                raise ValueError("NOT SPSS")
            elif max(quals) > 22: # Either RNA or invalid
                if max(quals) < 79:
                    raise ValueError("NOT SPSS")
                if min(quals) < 79:
                    raise ValueError("NOT SPSS")
                offset = 79
                is_rna = True
            elif min(quals) <  79: # Either DNA or invalid
                if min(quals) > 22:
                    raise ValueError("NOT SPSS")
                if max(quals) > 22:
                    raise ValueError("NOT SPSS")
                offset = 15
                is_rna = False
            else:
                raise ValueError("NOT SPSS")
            no_mod = list()
            no_indel = list()
            score_region = list()
            for qual in quals:
                #Convert qual score to int. Backfill 0s to make 3 bit int
                bits_string = str(bin(qual - offset))[2:].zfill(3)
                no_mod.append(int(bits_string[0] == '1'))
                no_indel.append(int(bits_string[1] == '1'))
                score_region.append(int(bits_string[2] == '1'))
            self.append(new_nucl_acid=nucl_acid(sequence=record.seq, no_mod=no_mod, no_indel=no_indel, score_region=score_region, is_rna=is_rna, design_parameters=design_parameters))

def mutate(nucl:nucl_acid, design_parameters:design_parameters) -> nucl_acid:
    """_summary_

    Args:
        nucl (nucl_acid): _description_
        design_parameters (design_parameters): _description_

    Raises:
        ValueError: _description_

    Returns:
        nucl_acid: _description_
    """
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

    sequence = list(copy.copy(str(nucl.sequence)))
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
        selection = 6
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
                #Othwerise, if we are at the end and it isn't scored, don't extend the score region
                else:
                    score_region.append(0)
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
