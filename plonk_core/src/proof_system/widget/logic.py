from dataclasses import dataclass
from bls12_381 import fr
from plonk_core.src.proof_system.mod import CustomEvaluations
from plonk_core.src.proof_system.widget.mod import WitnessValues,delta
from arithmetic import poly_mul_const
@dataclass
class LogicValues:
    # Left wire value in the next position
    a_next_val: fr.Fr
    # Right wire value in the next position
    b_next_val: fr.Fr
    # Fourth wire value in the next position
    d_next_val: fr.Fr
    # Constant selector value
    q_c_val: fr.Fr

    @staticmethod
    def from_evaluations(custom_evals:CustomEvaluations):
        a_next_val = custom_evals.get("a_next_eval")
        b_next_val = custom_evals.get("b_next_eval")
        d_next_val = custom_evals.get("d_next_eval")
        q_c_val = custom_evals.get("q_c_eval")
        return LogicValues(a_next_val,b_next_val,d_next_val,q_c_val)
    
class LogicGate:
    @staticmethod
    def constraints(separation_challenge:fr.Fr, wit_vals:WitnessValues, custom_vals:LogicValues):
        four = fr.Fr.from_repr(4)
        kappa = separation_challenge.square()
        kappa_sq = kappa.square()
        kappa_cu = kappa_sq.mul(kappa)
        kappa_qu = kappa_cu.mul(kappa)

        a_1 = four.mul(wit_vals.a_val)
        a = custom_vals.a_next_val.sub(a_1)
        c_0 = delta(a)

        b_1 = four.mul(wit_vals.b_val)
        b = custom_vals.b_next_val.sub(b_1)
        c_1 = delta(b)

        d_1 = four.mul(wit_vals.d_val)
        d = custom_vals.d_next_val.sub(d_1)
        c_2 = delta(d)

        w = wit_vals.c_val
        w_1 = a.mul(b)
        w_2 = w.sub(w_1)
        c_3 = w_2.mul(kappa_cu)

        c_4_1 = delta_xor_and(a,b,w,d,custom_vals.q_c_val)
        c_4 = c_4_1.mul(kappa_qu)

        mid1 = c_0.add(c_1)
        mid2 = mid1.add(c_2)
        mid3 = mid2.add(c_3)
        mid4 = mid3.add(c_4)
        res = mid4.mul(separation_challenge)
        return res
    
    @staticmethod
    def quotient_term(selector: fr.Fr, separation_challenge: fr.Fr, 
                      wit_vals: WitnessValues, custom_vals:LogicValues):
        temp = LogicGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = selector.mul(temp)
        return res
    
    @staticmethod
    def linearisation_term(selector_poly, separation_challenge, wit_vals, custom_vals):
        temp = LogicGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = poly_mul_const(selector_poly,temp)
        return res

# The identity we want to check is `q_logic * A = 0` where:
# A = B + E
# B = q_c * [9c - 3(a+b)]
# E = 3(a+b+c) - 2F
# F = w[w(4w - 18(a+b) + 81) + 18(a^2 + b^2) - 81(a+b) + 83]
def delta_xor_and(a: fr.Fr, b: fr.Fr, w: fr.Fr, c: fr.Fr, q_c: fr.Fr):
    nine = fr.Fr.from_repr(9)
    two = fr.Fr.from_repr(2)
    three = fr.Fr.from_repr(3)
    four = fr.Fr.from_repr(4)
    eighteen = fr.Fr.from_repr(18)
    eighty_one = fr.Fr.from_repr(81)
    eighty_three = fr.Fr.from_repr(83)

    f_1_1 = four.mul(w)
    f_1_2_1 = a.add(b)
    f_1_2 = eighteen.mul(f_1_2_1)
    f_1 = f_1_1.sub(f_1_2)
    f_1 = f_1.add(eighty_one)
    f_1 = f_1.mul(w)

    f_2_1_1 = a.square()
    f_2_1_2 = b.square()
    f_2_1 = f_2_1_1.add(f_2_1_2)
    f_2 = eighteen.mul(f_2_1)

    f_3_1 = a.add(b)
    f_3 = eighty_one.mul(f_3_1)

    f = f_1.add(f_2)
    f = f.sub(f_3)
    f = f.add(eighty_three)
    f = w.mul(f)

    e_1_1 = f_3_1.add(c)
    e_1 = three.mul(e_1_1)
    e_2 = two.mul(f)
    e = e_1.sub(e_2)

    b_1_1=nine.mul(c)
    b_1_2 = three.mul(f_3_1)
    b_1 = b_1_1.sub(b_1_2)
    b = q_c.mul(b_1)

    res = b.add(e)
    return res


