from dataclasses import dataclass
from typing import List
from field import field
from transcript import flags
@dataclass
class G2Coordinate:
    c0: any
    c1: any

@dataclass
class AffinePointG1:
    x: field
    y: field

    @classmethod
    def new(cls,x,y):
        return cls(x,y)

    #Returns the point at infinity, which always has Z = 0.
    @classmethod
    def zero(cls,params):
        x=field.zero(params)
        y=field.zero(params)
        return cls(x,y)
    
    def is_zero(self):
        return self.x.value == 0 and self.y.value == 0
    
    def serialize(self,writer):
        if self.is_zero():
            flag = flags.SWFlags.infinity()
            zero = field.zero(self.x.params)
            writer = zero.serialize_with_flags(writer,flag)
            return writer
        else:
            neg_y = self.y.neg()
            flag = flags.SWFlags.from_y_sign(self.y.value > neg_y.value)
            writer = self.x.serialize_with_flags(writer, flag)
            return writer
@dataclass
class AffinePointG2:
    x: G2Coordinate
    y: G2Coordinate

@dataclass
class UniversalParams:
    powers_of_g: List[AffinePointG1]
    powers_of_gamma_g: List[AffinePointG1]
    h: any
    beta_h: any


@dataclass
class OpenProof:
    # This is a commitment to the witness polynomial; see [KZG10] for more details.
    w: AffinePointG1
    # This is the evaluation of the random polynomial at the point for which
    # the evaluation proof was produced.
    random_v: any





