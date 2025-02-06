from Ultratherm.blist import blacklist
from Ultratherm.params import design_parameters, read_parameters
from Ultratherm.nucl import nucl_set, nucl_acid, nucl_hybrid
from Bio.Seq import Seq

import hashlib

def file_hash(filename: str) -> str:
    CHUNK_SIZE = 65536
    hash = hashlib.sha256()
    with open(filename, 'rb') as file:
        while next_bytes := file.read(CHUNK_SIZE):
            hash.update(next_bytes)
    return hash.hexdigest()

def test():
    blist = blacklist(path="blacklist.fasta")
    # score = vienna_score(sequence="CGAAAUCCCAACAGUGAAAACUUCCUCCAUGUUACAUAAUAGUAAGGAGGAAACAAAUG",
    #     score_region=[0]*42+[1]*8+[0]*9, design_parameters=design_parameters(blacklist=blist, target_temp=70,
    #     temp_offset=5, program="VIENNA", weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=1, num_mutants=8,
    #     target_energy=-12.0, free_energy_max_score=1.0 , accessibility_max_score=1.0, parasitic_complex_max_score=1.0), is_rna=True)
    # print(score)

    des_params = design_parameters(blacklist=blist, target_temp=55, temp_offset=4, program="VIENNA",
        weights=[8, 8, 8, 8, 10, 10, 16], weight_factor=2, num_mutants=8, target_energy=-8.0, # based on FourU Hairpin 2
        free_energy_max_score=0.9 , accessibility_max_score=0.9, parasitic_complex_max_score=0.9, thermo_score_temp=36)
    
    des_params.save('PARAMS.yml')
    test_parameters = read_parameters(path='PARAMS.yml')
    test_parameters.save('PARAMS2.yml')

    nucl_pool = nucl_set(nucls = [])
    for i in range(0, 16):
        new_nucl = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNNTAAGGAGGNNNNNNATG'),
            no_indel =      [0]*20+[1]*17,
            no_mod =        [0]*20+[1]*8+[0]*6+[1]*3,
            score_region =  [0]*20+[1]*8+[0]*6+[0]*3,
            is_rna=False)
        new_nucl.fitness_score(design_parameters=des_params)
        nucl_pool.append(new_nucl=new_nucl)
    
    nucl_pool.save("TEST" + '.fastq')
    del nucl_pool
    
    new_pool = nucl_set(nucls = [])
    new_pool.read("TEST.fastq", design_parameters=des_params)
    print(new_pool)

    new_nucl1 = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNNUAAGGAGGNNNNNNAUG'),
            no_indel =      [0]*20+[1]*17,
            no_mod =        [0]*20+[1]*8+[0]*6+[1]*3,
            score_region =  [0]*20+[1]*8+[0]*6+[0]*3,
            is_rna=True)
    new_nucl2 = nucl_acid(sequence=Seq('NNNNNNNNNNNNNNNNNNNNUAAGGAGGNNNNNN'),
            no_indel =      [0]*20+[1]*14,
            no_mod =        [0]*20+[1]*8+[0]*6,
            score_region =  [0]*20+[1]*8+[0]*6,
            is_rna=True)
    new_nucl_dimer = nucl_hybrid(new_nucl1, new_nucl2, True, False)
    new_nucl_dimer.fitness_score(design_parameters=des_params)
    new_pool.append(new_nucl_dimer)
    
    print(new_pool)
    print(new_nucl_dimer)

    new_pool.save("TEST2.fastq")
    newer_pool = nucl_set(nucls=[])
    newer_pool.read("TEST2.fastq", design_parameters=des_params)
    print(newer_pool)
    newer_pool.save("TEST3.fastq")

    assert file_hash("PARAMS.yml") == file_hash("PARAMS2.yml")
    assert file_hash("TEST2.fastq") == file_hash("TEST3.fastq")
    
test()
