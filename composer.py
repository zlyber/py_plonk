from dataclasses import dataclass
from typing import List, Dict

def next_power_of_2(x):
    return 1<<(x-1).bit_length()

@dataclass
class StandardComposer:
    # Number of arithmetic gates in the circuit
    n: int

    # Selector vectors
    q_m: List
    q_l: List
    q_r: List
    q_o: List
    q_4: List
    q_c: List

    # New selectors for poseidon hashes
    q_hl: List
    q_hr: List
    q_h4: List

    q_arith: List
    q_range: List
    q_logic: List
    q_fixed_group_add: List
    q_variable_group_add: List
    q_lookup: List

    # Sparse representation of Public Inputs
    intended_pi_pos: List[int]
    public_inputs: List
    # Witness vectors
    w_l: List
    w_r: List
    w_o: List
    w_4: List

    # Public lookup table
    lookup_table: List

    # A zero Variable
    zero_var: int


    def total_size(self):
        return max(self.n,len(self.lookup_table))
    
    def circuit_bound(self):
        return next_power_of_2(self.total_size())
