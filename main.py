from nucl import nucl_acid, nucl_set, nucl_hybrid
from blist import blacklist
from Bio.Seq import Seq
from params import design_parameters, read_parameters
from des import design

from signal import signal, SIGPIPE, SIG_IGN

# NOTE Customize these!
def rna_thermometer_prok():
    signal(SIGPIPE, SIG_IGN) # Ignore broken pipe (usually ssh) and continue program
    #Configure design parameters
    blist = blacklist(path="blacklist.fasta")
    des_params = design_parameters(blacklist=blist, target_temp=52, program='VIENNA',
        num_mutants=8, target_energy=-9.75, # based on FourU Hairpin 2
        weights=[32, 32, 32, 32, 32, 32, 16]
        )

    #Create nucleotide set
    nucl_pool = nucl_set(nucls = [])
    for i in range(0, 8):
        new_nucl = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNNUAAGGAGGNNNNNNAUG'),
            no_indel =      [0]*20+[1]*17,
            no_mod =        [0]*20+[1]*8+[0]*6+[1]*3,
            score_region =  [0]*20+[1]*8+[0]*6+[0]*3,
            is_rna=True)
        new_nucl.fitness_score(design_parameters=des_params)
        nucl_pool.append(new_nucl=new_nucl)

    #Start design loop
    design(design_parameters=des_params, nucl_pool=nucl_pool)

def rna_thermometer_euk_fiveprime():
    signal(SIGPIPE, SIG_IGN) # Ignore broken pipe (usually ssh) and continue program
    #Configure design parameters
    blist = blacklist(path="blacklist.fasta")
    des_params = design_parameters(blacklist=blist, target_temp=62, program='VIENNA',
        num_mutants=8, target_energy=-13.15,
        weights=[32, 32, 32, 32, 32, 32, 16]
        )

    #Create nucleotide set
    nucl_pool = nucl_set(nucls = [])
    for i in range(0, 8):
        new_nucl = nucl_acid(sequence=Seq('GAANNNNNNNNNNNNNNNNNNNN'), # Optimal for non-5'capped mRNAs (i.e. T7)
            no_indel =      [1]*3+[0]*20,
            no_mod =        [1]*3+[0]*20,
            score_region =  [1]*3+[0]*20,
            is_rna=True)
        new_nucl.fitness_score(design_parameters=des_params)
        nucl_pool.append(new_nucl=new_nucl)

    #Start design loop
    design(design_parameters=des_params, nucl_pool=nucl_pool)


def heteroduplex():
    signal(SIGPIPE, SIG_IGN) # Ignore broken pipe (usually ssh) and continue program
    #Configure design parameters
    blist = blacklist(path="blacklist.fasta")
    des_params = design_parameters(blacklist=blist, target_temp=62, program='VIENNA',
        num_mutants=8, target_energy=-13.15,
        weights=[32, 32, 32, 32, 32, 32, 16]
        )

    #Create nucleotide set
    nucl_pool = nucl_set(nucls = [])
    for i in range(0, 8):
        new_nucl_1 = nucl_acid(sequence=Seq('GAANNNNNNNNNNNNNNNNNNNN'), # Optimal for non-5'capped mRNAs (i.e. T7)
            no_indel =      [1]*3+[0]*20,
            no_mod =        [1]*3+[0]*20,
            score_region =  [1]*3+[0]*20,
            is_rna=True) # Limitation - heteroduplices not supported. Standard is to model as RNA
        new_nucl_2 = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNN'),
            no_indel =      [0]*20,
            no_mod =        [0]*20,
            score_region =  [0]*20,
            is_rna=True) # Limitation - heteroduplices not supported. Standard is to model as RNA
        new_nucl = nucl_hybrid(new_nucl_1, new_nucl_2, True, True)
        new_nucl.fitness_score(design_parameters=des_params)
        nucl_pool.append(new_nucl=new_nucl)

    #Start design loop
    design(design_parameters=des_params, nucl_pool=nucl_pool)

#####

###   #   #  #  #
#  #  #   #  ## #
##    #   #  # ##
#  #   ###   #  #

#####
if __name__ == '__main__':
    #rna_thermometer_prok()
    #rna_thermometer_euk_fiveprime()
    heteroduplex()
    pass
