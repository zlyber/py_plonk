from dataclasses import dataclass
from bls12_381 import fr
from plonk_core.src.proof_system.mod import CustomEvaluations
from plonk_core.src.proof_system.widget.mod import WitnessValues
from bls12_381.edwards import EdwardsParameters as P
from arithmetic import poly_mul_const
@dataclass
class FBSMValues:
    # Left wire value in the next position
    a_next_val: fr.Fr
    # Right wire value in the next position
    b_next_val: fr.Fr
    # Fourth wire value in the next position
    d_next_val: fr.Fr
    # Left selector value
    q_l_val: fr.Fr
    # Right selector value
    q_r_val: fr.Fr
    # Constant selector value
    q_c_val: fr.Fr

    @staticmethod
    def from_evaluations(custom_evals:CustomEvaluations):
        a_next_val = custom_evals.get("a_next_eval")
        b_next_val = custom_evals.get("b_next_eval")
        d_next_val = custom_evals.get("d_next_eval")
        q_l_val = custom_evals.get("q_l_eval")
        q_r_val = custom_evals.get("q_r_eval")
        q_c_val = custom_evals.get("q_c_eval")

        return FBSMValues(a_next_val,b_next_val,
            d_next_val,q_l_val,q_r_val,q_c_val)
    
class FBSMGate:
    @staticmethod
    def constraints(separation_challenge:fr.Fr, wit_vals:WitnessValues, custom_vals:FBSMValues):
        kappa = separation_challenge.square()
        kappa_sq = kappa.square()
        kappa_cu = kappa_sq.mul(kappa)

        x_beta_eval = custom_vals.q_l_val
        y_beta_eval = custom_vals.q_r_val

        acc_x = wit_vals.a_val
        acc_x_next = custom_vals.a_next_val
        acc_y = wit_vals.b_val
        acc_y_next = custom_vals.b_next_val

        xy_alpha = wit_vals.c_val

        accumulated_bit = wit_vals.d_val
        accumulated_bit_next = custom_vals.d_next_val
        bit = extract_bit(accumulated_bit, accumulated_bit_next)

        # Check bit consistency
        bit_consistency = check_bit_consistency(bit)

        one = y_beta_eval.one()
        y_beta_sub_one = y_beta_eval.sub(one)
        bit2 = bit.square()
        y_alpha_1 = bit2.mul(y_beta_sub_one)
        y_alpha = y_alpha_1.add(one)
        x_alpha = x_beta_eval.mul(bit)

        # xy_alpha consistency check
        bit_times_q_c_val = bit.mul(custom_vals.q_c_val)
        xy_consistency = bit_times_q_c_val.sub(xy_alpha)
        xy_consistency = xy_consistency.mul(kappa)

        # x accumulator consistency check
        x_3 = acc_x_next
        x_3_times_xy_alpha = x_3.mul(xy_alpha)
        x_3_times_xy_alpha_times_acc_x = x_3_times_xy_alpha.mul(acc_x)
        x_3_times_xy_alpha_times_acc_x_times_acc_y = x_3_times_xy_alpha_times_acc_x.mul(acc_y)
        x_3_times_xy_alpha_times_acc_x_times_acc_y_times_coeff_d = x_3_times_xy_alpha_times_acc_x_times_acc_y.mul(P.COEFF_D)
        lhs_x = x_3.add(x_3_times_xy_alpha_times_acc_x_times_acc_y_times_coeff_d)
        x_3_times_acc_y = x_alpha.mul(acc_y)
        y_alpha_times_acc_x = y_alpha.mul(acc_x)
        rhs_x = x_3_times_acc_y.add(y_alpha_times_acc_x) 
        x_acc_consistency = lhs_x.sub(rhs_x)
        x_acc_consistency = x_acc_consistency.mul(kappa_sq)

        # y accumulator consistency check
        y_3 = acc_y_next
        y_3_times_xy_alpha = y_3.mul(xy_alpha)
        y_3_times_xy_alpha_times_acc_x = y_3_times_xy_alpha.mul(acc_x)
        y_3_times_xy_alpha_times_acc_x_times_acc_y = y_3_times_xy_alpha_times_acc_x.mul(acc_y)
        y_3_times_xy_alpha_times_acc_x_times_acc_y_times_coeff_d = y_3_times_xy_alpha_times_acc_x_times_acc_y.mul(P.COEFF_D)
        lhs_y = y_3.sub(y_3_times_xy_alpha_times_acc_x_times_acc_y_times_coeff_d)
        y_alpha_times_acc_y = y_alpha.mul(acc_y)
        coeff_A_times_x_alpha = P.COEFF_A.mul(x_alpha)
        coeff_A_times_x_alpha_times_acc_x = coeff_A_times_x_alpha.mul(acc_x)
        rhs_y = y_alpha_times_acc_y.sub(coeff_A_times_x_alpha_times_acc_x)
        y_acc_consistency = lhs_y.sub(rhs_y)
        y_acc_consistency = y_acc_consistency.mul(kappa_cu)

        mid1 = bit_consistency.add(x_acc_consistency)
        mid2 = mid1.add(y_acc_consistency)
        checks = mid2.add(xy_consistency)
        res = checks.mul(separation_challenge)
        return res

    @staticmethod
    def quotient_term(selector: fr.Fr, separation_challenge: fr.Fr, 
                      wit_vals: WitnessValues, custom_vals:FBSMValues):
        temp = FBSMGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = selector.mul(temp)
        return res
    
    @staticmethod
    def linearisation_term(selector_poly, separation_challenge, wit_vals, custom_vals):
        temp = FBSMGate.constraints(separation_challenge, wit_vals, custom_vals)
        res = poly_mul_const(selector_poly,temp)
        return res
# Extracts the bit value from the accumulated bit.
def extract_bit(curr_acc: fr.Fr, next_acc: fr.Fr):
    mid1 = next_acc.sub(curr_acc)
    res = mid1.sub(curr_acc)
    return res


# Ensures that the bit is either `+1`, `-1`, or `0`.
def check_bit_consistency(bit: fr.Fr):
    one = bit.one()
    mid1 = bit.sub(one)
    mid2 = bit.add(one)
    res = bit.mul(mid1)
    res = res.mul(mid2)
    return res
