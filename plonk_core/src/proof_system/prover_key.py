from dataclasses import dataclass
from field import field
from typing import List, Tuple
from plonk_core.src.proof_system.widget.arithmetic import Arith
from plonk_core.src.proof_system.widget.lookup import Lookup
from plonk_core.src.proof_system.permutation import Permutation
@dataclass 
class Prover_Key:
    arithmetic: Arith

    range_selector: Tuple[List[field],List[field]]

    logic_selector: Tuple[List[field],List[field]]

    lookup: Lookup

    fixed_group_add_selector: Tuple[List[field],List[field]]

    variable_group_add_selector: Tuple[List[field],List[field]]

    permutation: Permutation

    v_h_coset_8n: List[field]