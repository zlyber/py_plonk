from dataclasses import dataclass
from bls12_381 import fr
from plonk_core.src.proof_system.mod import CustomEvaluations
from plonk_core.src.proof_system.widget.mod import WitnessValues
from bls12_381.edwards import EdwardsParameters as P
from arithmetic import poly_mul_const
@dataclass
class CAValues:
    # Left wire value in the next position
    a_next_val: fr.Fr
    # Right wire value in the next position
    b_next_val: fr.Fr
    # Fourth wire value in the next position
    d_next_val: fr.Fr

    @staticmethod
    def from_evaluations(custom_evals:CustomEvaluations):
        a_next_val = custom_evals.get("a_next_eval")
        b_next_val = custom_evals.get("b_next_eval")
        d_next_val = custom_evals.get("d_next_eval")

        return CAValues(a_next_val,b_next_val,d_next_val)
    
class CAGate:
    @staticmethod
    def constraints(separation_challenge: fr.Fr, wit_vals: WitnessValues, custom_vals: CAValues):
        x_1 = wit_vals.a_val
        x_3 = custom_vals.a_next_val
        y_1 = wit_vals.b_val
        y_3 = custom_vals.b_next_val
        x_2 = wit_vals.c_val
        y_2 = wit_vals.d_val
        x1_y2 = custom_vals.d_next_val

        kappa = separation_challenge.square()

        # Check that `x1 * y2` is correct
        x1y2 = x_1.mul(y_2)
        xy_consistency = x1y2.sub(x1_y2)

        y1_x2 = y_1.mul(x_2)
        y1_y2 = y_1.mul(y_2)
        x1_x2 = x_1.mul(x_2)

        # Check that `x_3` is correct
        x3_lhs = x1_y2.add(y1_x2)
        x_3_D = x_3.mul(P.COEFF_D)
        x_3_D_x1_y2 = x_3_D.mul(x1_y2)
        x_3_D_x1_y2_y1_x2 = x_3_D_x1_y2.mul(y1_x2)
        x3_rhs = x_3.add(x_3_D_x1_y2_y1_x2)
        x3_l_sub_r = x3_lhs.sub(x3_rhs)
        x3_consistency = x3_l_sub_r.mul(kappa)

        # Check that `y_3` is correct
        x1_x2_A = P.COEFF_A.mul(x1_x2)
        y3_lhs = y1_y2.sub(x1_x2_A)
        y_3_D = y_3.mul(P.COEFF_D)
        y_3_D_x1_y2 = y_3_D.mul(x1_y2)
        y_3_D_x1_y2_y1_x2 = y_3_D_x1_y2.mul(y1_x2)
        y3_rhs = y_3.sub(y_3_D_x1_y2_y1_x2)
        y3_l_sub_r = y3_lhs.sub(y3_rhs)
        kappa2 = kappa.square()
        y3_consistency = y3_l_sub_r.mul(kappa2)

        mid1 = xy_consistency.add(x3_consistency)
        mid2 = mid1.add(y3_consistency)
        result = mid2.mul(separation_challenge)
        return result


    @staticmethod
    def quotient_term(selector: fr.Fr, separation_challenge: fr.Fr, 
                      wit_vals: WitnessValues, custom_vals:CAValues):
        temp = CAGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = selector.mul(temp)
        return res
    
    @staticmethod
    def linearisation_term(selector_poly, separation_challenge, wit_vals, custom_vals):
        temp = CAGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = poly_mul_const(selector_poly,temp)
        return res
