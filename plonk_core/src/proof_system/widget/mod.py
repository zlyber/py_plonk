from field import field
from dataclasses import dataclass

@dataclass
class WitnessValues:
    a_val: field  # Left Value
    b_val: field  # Right Value
    c_val: field  # Output Value
    d_val: field  # Fourth Value


def delta(f:field):
    one = f.one()
    two = field.from_repr(2,f.params)
    three = field.from_repr(3,f.params)
    f_1 = f.sub(one)
    f_2 = f.sub(two)
    f_3 = f.sub(three)
    mid1 = f_1.mul(f_2)
    mid2 = mid1.mul(f_3)
    res = f.mul(mid2)
    return res