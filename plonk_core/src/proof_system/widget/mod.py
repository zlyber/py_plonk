from bls12_381 import fr
from dataclasses import dataclass

@dataclass
class WitnessValues:
    a_val: fr.Fr  # Left Value
    b_val: fr.Fr  # Right Value
    c_val: fr.Fr  # Output Value
    d_val: fr.Fr  # Fourth Value


def delta(f:fr.Fr):
    one = f.one()
    two = fr.Fr.from_repr(2)
    three = fr.Fr.from_repr(3)
    f_1 = f.sub(one)
    f_2 = f.sub(two)
    f_3 = f.sub(three)
    mid1 = f_1.mul(f_2)
    mid2 = mid1.mul(f_3)
    res = f.mul(mid2)
    return res