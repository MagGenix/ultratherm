from nucl import nucl_acid, nucl_set, nucl_hybrid
from blist import blacklist
from Bio.Seq import Seq
from params import design_parameters, read_parameters
from des import design

# NOTE Customize these!
def rna_thermometer_prok():
    #Configure design parameters
    blist = blacklist(path="blacklist.fasta")
    des_params = design_parameters(blacklist=blist, target_temp=42, program='NUPACK',
        num_mutants=8, target_energy=-6, # based on http://dx.doi.org/10.1101/017269
        weights=[32, 32, 32, 32, 32, 32, 16],
        thermo_score_temp=30,
        max_reps=16,
        max_hairpins=1,
        accessibility_max_score=3,
        parasitic_complex_max_score=0,
        optimization_rate=3,
        temp_offset=10
        )

    #Create nucleotide set
    nucl_pool = nucl_set(nucls = [])
    for i in range(0, 64):
        new_nucl = nucl_acid(sequence=Seq('GNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAAGGAGGNNNNNNAUGAGUAAAGGC'),
            no_indel =      [0]+[0]*40+[1]*26,
            no_mod =        [1]+[0]*40+[1]*8+[0]*6+[1]*12,
            score_region =  [0]+[0]*40+[1]*8+[0]*6+[0]*12,
            is_rna=True)
        new_nucl.fitness_score(design_parameters=des_params)
        nucl_pool.append(new_nucl=new_nucl)

    #Start design loop
    design(design_parameters=des_params, nucl_pool=nucl_pool)

def rna_thermometer_euk_fiveprime():
    #Configure design parameters
    blist = blacklist(path="blacklist.fasta")
    des_params = design_parameters(blacklist=blist, target_temp=62, program='VIENNA',
        num_mutants=8, target_energy=-13.15,
        weights=[32, 32, 32, 32, 32, 32, 16], optimization_rate=10
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
    #Configure design parameters
    blist = blacklist(path="blacklist.fasta")
    des_params = design_parameters(blacklist=blist, target_temp=62, program='VIENNA',
        num_mutants=8, target_energy=-13.15,
        weights=[32, 32, 32, 32, 32, 32, 16], optimization_rate=10
        )

    #Create nucleotide set
    nucl_pool = nucl_set(nucls = [])
    for i in range(0, 8):
        new_nucl_1 = nucl_acid(sequence=Seq('GAANNNNNNNNNNNNNNNNN'), # Optimal for non-5'capped mRNAs (i.e. T7)
            no_indel =      [1]*3+[0]*17,
            no_mod =        [1]*3+[0]*17,
            score_region =  [1]*3+[0]*17,
            is_rna=True) # Limitation - heteroduplices not supported. Standard is to model as RNA
        new_nucl_2 = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNN'),
            no_indel =      [0]*20,
            no_mod =        [0]*20,
            score_region =  [0]*20,
            is_rna=True) # Limitation - heteroduplices not supported. Standard is to model as RNA
        new_nucl = nucl_hybrid(new_nucl_1, new_nucl_2, True, False)
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
    rna_thermometer_prok()
    #rna_thermometer_euk_fiveprime()
    #heteroduplex()
    pass
