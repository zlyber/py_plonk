from dataclasses import dataclass
from bls12_381 import fr
from typing import List, Tuple
from plonk_core.src.proof_system.widget.mod import WitnessValues
from plonk_core.src.constraint_system.hash import SBOX_ALPHA
from arithmetic import poly_mul_const,poly_add_poly
@dataclass
class Arith:
    q_m: Tuple[List,List]
    q_l: Tuple[List,List]
    q_r: Tuple[List,List]
    q_o: Tuple[List,List]
    q_4: Tuple[List,List]
    q_hl: Tuple[List,List]
    q_hr: Tuple[List,List]
    q_h4: Tuple[List,List]
    q_c: Tuple[List,List]
    q_arith: Tuple[List,List]

    # Computes the arithmetic gate contribution to the quotient polynomial at
    # the element of the domain at the given `index`.
    def compute_quotient_i(self, index: int, wit_vals: WitnessValues):

        mult = wit_vals.a_val.mul(wit_vals.b_val)
        mult = mult.mul(self.q_m[1][index])
        left = wit_vals.a_val.mul(self.q_l[1][index])
        right = wit_vals.b_val.mul(self.q_r[1][index])
        out = wit_vals.c_val.mul(self.q_o[1][index])
        fourth = wit_vals.d_val.mul(self.q_4[1][index])
        a_high = wit_vals.a_val.pow(SBOX_ALPHA)
        b_high = wit_vals.b_val.pow(SBOX_ALPHA)
        f_high = wit_vals.d_val.pow(SBOX_ALPHA)
        a_high = a_high.mul(self.q_hl[1][index])
        b_high = b_high.mul(self.q_hr[1][index])
        f_high = f_high.mul(self.q_h4[1][index])
        
        mid1 = mult.add(left)
        mid2 = mid1.add(right)
        mid3 = mid2.add(out)
        mid4 = mid3.add(fourth)
        mid5 = mid4.add(a_high)
        mid6 = mid5.add(b_high)
        mid7 = mid6.add(f_high)
        mid8 = mid7.add(self.q_c[1][index])

        arith_val = mid8.mul(self.q_arith[1][index])
        return arith_val
    # Computes the arithmetic gate contribution to the linearisation
    # polynomial at the given evaluation points.
    def compute_linearisation(
        self, 
        a_eval: fr.Fr,
        b_eval: fr.Fr, 
        c_eval: fr.Fr, 
        d_eval: fr.Fr, 
        q_arith_eval: fr.Fr):
        mid1_1 = a_eval.mul(b_eval)
        mid1 = poly_mul_const(self.q_m[0] ,mid1_1)
        mid2 = poly_mul_const(self.q_l[0] ,a_eval)
        mid3 = poly_mul_const(self.q_r[0] ,b_eval)
        mid4 = poly_mul_const(self.q_o[0] ,c_eval)
        mid5 = poly_mul_const(self.q_4[0] ,d_eval)
        mid6_1 = a_eval.pow(SBOX_ALPHA)
        mid6 = poly_mul_const(self.q_hl[0] ,mid6_1)
        mid7_1 = b_eval.pow(SBOX_ALPHA)
        mid7 = poly_mul_const(self.q_hr[0] ,mid7_1)
        mid8_1 = d_eval.pow(SBOX_ALPHA)
        mid8 = poly_mul_const(self.q_h4[0] ,mid8_1)

        add1 = poly_add_poly(mid1, mid2)
        add2 = poly_add_poly(add1, mid3)
        add3 = poly_add_poly(add2, mid4)
        add4 = poly_add_poly(add3, mid5)
        add5 = poly_add_poly(add4, mid6)
        add6 = poly_add_poly(add5, mid7)
        add7 = poly_add_poly(add6, mid8)
        add8 = poly_add_poly(add7, self.q_c[0])

        result = poly_mul_const(add8, q_arith_eval)
        return result
