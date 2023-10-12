from field import field
from plonk_core.src.proof_system.mod import CustomEvaluations
from plonk_core.src.proof_system.widget.mod import WitnessValues,delta
from arithmetic import poly_mul_const

class RangeValues:
    def __init__(self, d_next_val:field):
        self.d_next_val = d_next_val
        
    @staticmethod
    def from_evaluations(custom_vals:CustomEvaluations):
        d_next_val = custom_vals.get("d_next_eval")
        return RangeValues(d_next_val)

class RangeGate:

    @staticmethod
    def constraints(separation_challenge:field, wit_vals:WitnessValues, custom_vals:RangeValues):
        four = field.from_repr(4,separation_challenge.params)
        kappa = separation_challenge.square()
        kappa_sq = kappa.square()
        kappa_cu = kappa_sq.mul(kappa)

        b_1_1 = four.mul(wit_vals.d_val)
        f_b1 = wit_vals.c_val.sub(b_1_1)
        b_1 = delta(f_b1)

        b_2_1 = four.mul(wit_vals.c_val)
        b_2_2 = wit_vals.b_val.sub(b_2_1)
        f_b2 = delta(b_2_2)
        b_2 = f_b2.mul(kappa)

        b_3_1 = four.mul(wit_vals.b_val)
        b_3_2 = wit_vals.a_val.sub(b_3_1)
        f_b3 = delta(b_3_2)
        b_3 = f_b3.mul(kappa_sq)

        b_4_1 = four.mul(wit_vals.a_val)
        b_4_2 = custom_vals.d_next_val.sub(b_4_1)
        f_b4 = delta(b_4_2)
        b_4 = f_b4.mul(kappa_cu)

        mid1 = b_1.add(b_2)
        mid2 = mid1.add(b_3)
        mid3 = mid2.add(b_4)
        res = mid3.mul(separation_challenge)

        return res
    
    @staticmethod
    def quotient_term(selector, separation_challenge, wit_vals, custom_vals):
        temp = RangeGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = selector.mul(temp)
        return res
    
    @staticmethod
    def linearisation_term(selector_poly, separation_challenge, wit_vals, custom_vals):
        temp = RangeGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = poly_mul_const(selector_poly,temp)
        return res