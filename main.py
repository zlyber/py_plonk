import time
from load import read_pk_data,read_pp_data,read_cs_data
from composer import StandardComposer
import gen_proof
from transcript import transcript

if __name__ == "__main__":

    pp_file = "params.txt"
    pk_file = "pk.txt"
    cs_file = "cs.txt"

    start_time = time.time()
    pp = read_pp_data(pp_file)
    pk = read_pk_data(pk_file)
    csdata = read_cs_data(cs_file)
    end_time = time.time()
    load_time = end_time - start_time
    print(f"load time: {load_time} s")
    cs=StandardComposer(n=csdata["n"],q_m=csdata["q_m"],q_l=csdata["q_l"],q_r=csdata["q_r"],
                        q_o=csdata["q_o"],q_4=csdata["q_4"],q_c=csdata["q_c"],q_hl=csdata["q_hl"],
                        q_hr=csdata["q_hr"],q_h4=csdata["q_h4"],q_arith=csdata["q_arith"],
                        q_range=csdata["q_range"],q_logic=csdata["q_logic"],
                        q_fixed_group_add=csdata["q_fixed"],public_inputs=csdata["public_inputs"],
                        q_variable_group_add=csdata["q_variable"],
                        q_lookup=csdata["q_lookup"],intended_pi_pos=csdata["intended_pi_pos"],
                        w_l=csdata["w_l"],w_r=csdata["w_r"],w_o=csdata["w_o"],w_4=csdata["w_4"],
                        lookup_table=csdata["lookup_table"],zero_var=csdata["zero_var"])

    
    
    transcript_init = b"Merkle tree"
    preprocessed_transcript = transcript.Transcript.new(transcript_init)
    start_time = time.time()
    pi = gen_proof.gen_proof(pp,pk,cs,preprocessed_transcript)
    end_time = time.time()
    print("Generate proof successfully\n")
    execution_time = end_time - start_time
    print(f"execution time: {execution_time} s")
