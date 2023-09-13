from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import copy
from typing import Union

from params import design_parameters
from blist import blacklist
from nupack_score import nupack_score
from vienna_score import vienna_score
from nupack_score_hybrid import nupack_score_hybrid
from vienna_score_hybrid import vienna_score_hybrid

class nucl_acid():
    """nucleic acid. Stores its sequence, score, no_mod, no_no_indel, score_region, and whether it is RNA or DNA.
    """
    def __init__(self, sequence: Seq, no_mod: list, no_indel: list, score_region: list, is_rna: bool):
        """Create a new nucl_acid.

        Args:
            sequence (Seq): the nucleic acid sequence.
            no_mod (list): A list() of 0 and 1 (ints), where 1 denotes nucleotides that CANNOT be modified (neither indels nor substitutions).
            no_indel (list): A list() of 0 and 1 (ints), where 1 denotes nucleotides that CANNOT have indel mutations adjacent to them.
            score_region (list): A list() of 0 and 1 (ints), where 1 denotes nucleotides that will be scored for high pair probability at temp - offset and low pair probability at temp + offset.
            design_parameters (design_parameters): the design parameters.
            is_rna (bool): whether the nucleic acid is DNA or RNA.

        Raises:
            Exception: _description_
            Exception: _description_
            Exception: _description_
            ValueError: _description_
            ValueError: _description_
        """
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

        self.score = None
    
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
        """Called during initialization to score the nucl_acid.

        Args:
            design_parameters (design_parameters): the design parameters.

        Raises:
            Exception: _description_
        """
        if self.is_blacklisted(blacklist=design_parameters.blacklist):
            self.score = 6.0
        elif design_parameters.program == "NUPACK":  
            self.score = nupack_score(sequence=str(self.sequence), score_region=self.score_region, is_rna=self.is_rna, design_parameters=design_parameters)
        elif design_parameters.program == "VIENNA":
            self.score = vienna_score(sequence=str(self.sequence), score_region=self.score_region, is_rna=self.is_rna, design_parameters=design_parameters)
        else:
            raise Exception("no program specified for scoring")

    #As of right now, this function is only accessed from within fitness_score, but could be useful to external processes.
    #If used in external processes it would be more efficent to store it as a boolean member variable.
    #However, since member variables are public in Python, this would make is_blacklisted modifiable and potentially inaccurate.
    #For that reason, and since it should only be performed once but is also generally important,
    #this stays a function for now.
    def is_blacklisted(self, blacklist: blacklist) -> bool:
        """Called during initialization to determine whether the nucl_acid is blacklisted.

        Args:
            blacklist (blacklist): the blacklist object.

        Returns:
            bool: whether or not the nucl_acid is blacklisted.
        """
        if blacklist.is_empty:
            return False
        if self.sequence in blacklist.blacklist_sequences:
            return True
        return False

class nucl_hybrid():
    def __init__(self, nucl_1: nucl_acid, nucl_2:nucl_acid):
        self.nucl_1 = nucl_1
        self.nucl_2 = nucl_2
        self.score = None
    def fitness_score(self, design_parameters: design_parameters):
        if self.is_blacklisted(blacklist=design_parameters.blacklist):
            self.score = 6.0
        
        elif design_parameters.program == "NUPACK":  
            self.score = nupack_score_hybrid(sequence_1=str(self.nucl_1.sequence),
                                             score_region_1=self.nucl_1.score_region,
                                             is_rna_1=self.nucl_1.is_rna,
                                             sequence_2=str(self.nucl_2.sequence),
                                             score_region_2=self.nucl_2.score_region,
                                             is_rna_2=self.nucl_2.is_rna,
                                             design_parameters=design_parameters)
        elif design_parameters.program == "VIENNA":
            self.score = vienna_score_hybrid(sequence_1=str(self.nucl_1.sequence),
                                             score_region_1=self.nucl_1.score_region,
                                             is_rna_1=self.nucl_1.is_rna,
                                             sequence_2=str(self.nucl_2.sequence),
                                             score_region_2=self.nucl_2.score_region,
                                             is_rna_2=self.nucl_2.is_rna,
                                             design_parameters=design_parameters)
        else:
            raise Exception("no program specified for scoring")

    def __str__(self) -> str:
        as_string = ""
        as_string += str(self.score) + "\n"
        as_string += str(self.nucl_1.sequence)      + "\t" + str(self.nucl_2.sequence)      + "\n"
        as_string += str(self.nucl_1.no_mod)        + "\t" + str(self.nucl_2.no_mod)        + "\n"
        as_string += str(self.nucl_1.no_indel)      + "\t" + str(self.nucl_2.no_indel)      + "\n"
        as_string += str(self.nucl_1.score_region)  + "\t" + str(self.nucl_2.score_region)  + "\n"
        return as_string

    def is_blacklisted(self, blacklist: blacklist) -> bool:
        if blacklist.is_empty:
            return False
        if self.nucl_1.is_blacklisted or self.nucl_2.is_blacklisted:
            return True
        return False
    
class nucl_set():
    """A data type to store nucl_acid's for design. Stores ordered lists of nucl_acid and their scores.
    Use .append(), .remove(), and .replace() to modify.
    Use .save() and .read() for I/O.

    DO NOT DIRECTLY MODIFY MEMBERS! They are intended for read-only. Use the above methods for modification.
    """
    def __init__(self, nucls: list):
        """Create a new nucl_set.

        Args:
            nucls (list): a list of nucl_acid and/or nucl_hybrid.

        Raises:
            TypeError: A member of the list was not a nucl_acid.
        """
        self.nucls = nucls
        self.scores = []
        for i in range (0, len(self.nucls)):
            if type(self.nucls[i]) != nucl_acid and type(self.nucls[i]) != nucl_hybrid:
                #Raise an error if an object in the list was not a nucl_acid or nucl_hybrid
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
        """Replace a nucl_acid in the nucl_set.

        Args:
            index (int): the index of the nucl_acid in nucls[].
            new_nucl_acid (nucl_acid): a new nucl_acid object to insert at the index.

        Raises:
            IndexError: _description_
            IndexError: _description_
        """
        if index < 0:
            raise IndexError
        
        #nucl_set is zero-indexed!
        if index > len(self) - 1:
            raise IndexError
        
        if new_nucl_acid.score == None:
            raise ValueError # scores array MUST support comparison!
        
        self.nucls[index] = new_nucl_acid
        self.scores[index] = new_nucl_acid.score

    def append(self, new_nucl: Union[nucl_acid, nucl_hybrid]) -> None:
        """Append a nucl_acid or nucl_hybrid to the end of the nucl_set.

        Args:
            new_nucl_acid (nucl_acid): a nucl_acid object to append.
        """
        if new_nucl.score == None:
            raise ValueError # scores array MUST support comparison!
        
        self.nucls.append(new_nucl)
        self.scores.append(new_nucl.score)

    def remove(self, index:int) -> None:
        """Remove a nucl_acid from the nucl_set.

        Args:
            index (int): the index of the nucl_acid object in nucls[] to be removed.
        """
        del self.nucls[index]
        del self.scores[index]

    def save(self, path:str) -> None:
        """Write the nucl_set to an SPSS formatted .fastq file.

        Args:
            path (str): Path to save the file. If no directory specified will save to source directory.
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

                # If a nucl_acid, follow the basic protocol to write to a FASTQ record
                if type(nucl) == nucl_acid:
                    no_mod = nucl.no_mod
                    no_indel = nucl.no_indel
                    score_region = nucl.score_region
                    
                    sequence = nucl.sequence

                    length = len(nucl)

                    if nucl.is_rna:
                        offset = [79] * length
                    else:
                        offset = [15] * length

                    score = nucl.score
                # Combine and represent the nucl_hybrid so it can be printed to one FASTQ record
                elif type(nucl) == nucl_hybrid:
                    no_mod =        nucl.nucl_1.no_mod +        [-1] + nucl.nucl_2.no_mod # TODO fix this fucking shit
                    no_indel =      nucl.nucl_1.no_indel +      [-1] + nucl.nucl_2.no_mod
                    score_region =  nucl.nucl_2.score_region +  [-1] + nucl.nucl_2.score_region
                    
                    if nucl.nucl_1.is_rna:
                        offset1 = 79
                    else:
                        offset1 = 15
                    
                    if nucl.nucl_2.is_rna:
                        offset2 = 79
                    else:
                        offset2 = 15
                    
                    sequence = nucl.nucl_1.sequence + "&" + nucl.nucl_2.sequence

                    length = len(nucl.nucl_1) + 1 + len(nucl.nucl_2)
                    offset = [offset1] * len(nucl.nucl_1) + [-1] + [offset2] * len(nucl.nucl_2)

                    score = nucl.score
                    
                # Element was not a nucl_acid or nucl_hybrid!
                else:
                    raise TypeError

                quals = list()
                for j in range (0, length):
                    if no_mod[j] == -1:
                        quals.append(5) # The '&' character, denoting a hybrid of 2 strands
                        continue
                    bitString = str(no_mod[j]) + str(no_indel[j]) + str(score_region[j])
                    bitsAsInt = int(bitString, 2)
                    quals.append(bitsAsInt + offset[j])
                record = SeqRecord(seq=sequence, id=str(i), description=str(score))
                record.letter_annotations["phred_quality"] = quals
                SeqIO.write(record, handle=handle, format='fastq')

                # TODO do these variables need to be cleared every time ? ? ?
                del quals
                del record
                del nucl

    def read(self, path:str, design_parameters:design_parameters) -> None:
        """Reads an SPSS formatted .fastq file and appends its nucl_acid's to the nucl_set.

        Args:
            path (str): Path to read the file. If no directory specified will assume source directory.
            design_parameters (design_parameters): the design parameters.

        Raises:
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
        """
        # TODO implement reading for nucl_hybrid's
        #If any other characters are encountered, an error should be raised.
        for record in SeqIO.parse(path, "fastq"):
            quals = record.letter_annotations["phred_quality"]
            sequence = record.seq
            strands = list()
            nucls = list()
            if sequence.count("&") > 0:
                if sequence.count("&") > 1:
                    raise ValueError("MORE THAN 2 STRANDS DETECTED")
                if quals.count(5) != sequence.count("&"):
                    raise ValueError("NOT SPSS")
                
                split_point = quals.index(5)

                strands.append((sequence[0:split_point],    quals[0:split_point]))
                strands.append((sequence[split_point+1:],   quals[split_point+1:]))

            else:
                strands.append((sequence, quals))
            
            for strand in strands:
                strand_seq =    strand[0]
                strand_quals =  strand[1]

                if max(strand_quals) > 86 or min(strand_quals) < 15:
                    raise ValueError("NOT SPSS")
                elif max(strand_quals) > 22: # Either RNA or invalid
                    if max(strand_quals) < 79:
                        raise ValueError("NOT SPSS")
                    if min(strand_quals) < 79:
                        raise ValueError("NOT SPSS")
                    offset = 79
                    is_rna = True
                elif min(strand_quals) <  79: # Either DNA or invalid
                    if min(strand_quals) > 22:
                        raise ValueError("NOT SPSS")
                    if max(strand_quals) > 22:
                        raise ValueError("NOT SPSS")
                    offset = 15
                    is_rna = False
                else:
                    raise ValueError("NOT SPSS")
                no_mod = list()
                no_indel = list()
                score_region = list()
                for qual in strand_quals:
                    #Convert qual score to int. Backfill 0s to make 3 bit int
                    bits_string = str(bin(qual - offset))[2:].zfill(3)
                    no_mod.append(int(bits_string[0] == '1'))
                    no_indel.append(int(bits_string[1] == '1'))
                    score_region.append(int(bits_string[2] == '1'))
                
                nucls.append(nucl_acid(sequence=strand_seq, no_mod=no_mod, no_indel=no_indel, score_region=score_region, is_rna=is_rna))
            if len(strands) == 1:
                nucls[0].fitness_score(design_parameters=design_parameters)
                self.append(new_nucl=nucls[0])
            elif len(strands) == 2:
                new_hybrid = nucl_hybrid(nucls[0], nucls[1])
                new_hybrid.fitness_score(design_parameters=design_parameters)
                self.append(new_hybrid)

def mutate(nucl: Union[nucl_acid, nucl_hybrid], design_parameters:design_parameters) -> Union[nucl_acid, nucl_hybrid]:
    """Mutates a nucl_acid given design parameters and returns a mutated, scored nucl_acid. Does not modify the original.

    Args:
        nucl (nucl_acid): the nucl_acid to make a mutant of.
        design_parameters (design_parameters): the design parameters.

    Raises:
        ValueError: _description_

    Returns:
        nucl_acid: a new nucl_acid.
    """
    #I'm trying to make this function as fast as possible since it will be called once per every single
    #nucleotide in a sequence to generate one variant.
    #For larger pool sizes this can really stack up.
    
    if type(nucl) == nucl_hybrid:
        nucl1 = mutate(nucl.nucl_1, design_parameters=design_parameters)
        nucl2 = mutate(nucl.nucl_2, design_parameters=design_parameters)
        return nucl_hybrid(nucl_1=nucl1, nucl_2=nucl2) # type: ignore
    
    elif type(nucl) == nucl_acid:
        if len(nucl.no_mod) != len(nucl.sequence) or len(nucl.no_indel) != len(nucl.score_region) or len(nucl.no_mod) != len(nucl.no_indel):
            raise ValueError

        #As of writing the design parameters function forces the weights to be ints above 0.
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
        
        new_nucl_acid = nucl_acid(sequence=Seq(''.join(sequence)), no_mod=no_mod, no_indel=no_indel, score_region=score_region,
                                is_rna=nucl.is_rna)
        return new_nucl_acid
    
    else:
        raise TypeError
