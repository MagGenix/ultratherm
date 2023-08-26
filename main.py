from nucl import nucl_acid, nucl_set
from blist import blacklist
from Bio.Seq import Seq
from params import design_parameters, read_parameters
from des import design
from vienna_score import vienna_score

# NOTE Customize this!
def main():
    #Configure design parameters
    blist = blacklist(path="blacklist.fasta")
    des_params = design_parameters(blacklist=blist, target_temp=55, program='NUPACK',
        num_mutants=2, target_energy=-8.0, # based on FourU Hairpin 2
        )

    #Create nucleotide set
    pool = nucl_set(nucls = [])
    for i in range(0, 4):
        pool.append(nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNNUAAGGAGGNNNNNNAUG'),
            no_indel =      [0]*20+[1]*17,
            no_mod =        [0]*20+[1]*8+[0]*6+[1]*3,
            score_region =  [0]*20+[1]*8+[0]*6+[0]*3,
            design_parameters=des_params, is_rna=True))

    #Start design loop
    design(design_parameters=des_params, pool=pool)

def test():
    blist = blacklist(path="blacklist.fasta")
    # score = vienna_score(sequence="CGAAAUCCCAACAGUGAAAACUUCCUCCAUGUUACAUAAUAGUAAGGAGGAAACAAAUG",
    #     score_region=[0]*42+[1]*8+[0]*9, design_parameters=design_parameters(blacklist=blist, target_temp=70,
    #     temp_offset=5, program="VIENNA", weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=1, num_mutants=8,
    #     target_energy=-12.0, free_energy_max_score=1.0 , nucl_max_score=1.0, max_dimer_monomer_factor=1.0), is_rna=True)
    # print(score)

    des_params = design_parameters(blacklist=blist, target_temp=55, temp_offset=4, program="NUPACK",
        weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=2, num_mutants=8, target_energy=-8.0, # based on FourU Hairpin 2
        free_energy_max_score=0.9 , nucl_max_score=0.9, max_dimer_monomer_factor=0.9, thermo_score_temp=36)
    
    des_params.save('PARAMS.yml')
    test_parameters = read_parameters(path='PARAMS.yml')
    test_parameters.save('PARAMS2.yml')

    pool = nucl_set(nucls = [])
    for i in range(0, 16):
        pool.append(nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNNTAAGGAGGNNNNNNATG'),
            no_indel =      [0]*20+[1]*17,
            no_mod =        [0]*20+[1]*8+[0]*6+[1]*3,
            score_region =  [0]*20+[1]*8+[0]*6+[0]*3,
            design_parameters=des_params, is_rna=False))

    pool.save("TEST" + '.fastq')
    del pool
    
    new_pool = nucl_set(nucls = [])
    new_pool.read("TEST.fastq", design_parameters=des_params)
    print(new_pool)

#####

###   #   #  #  #
#  #  #   #  ## #
##    #   #  # ##
#  #   ###   #  #

#####

main()
#test()
