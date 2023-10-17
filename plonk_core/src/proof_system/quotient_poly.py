from domain import Radix2EvaluationDomain
from field import field
from arithmetic import NTT,INTT,coset_NTT,coset_INTT,from_coeff_vec
from plonk_core.src.proof_system.widget.mod import WitnessValues
from plonk_core.src.proof_system.widget.range import RangeGate,RangeValues
from plonk_core.src.proof_system.widget.logic import LogicGate,LogicValues
from plonk_core.src.proof_system.widget.fixed_base_scalar_mul import FBSMGate,FBSMValues
from plonk_core.src.proof_system.widget.curve_addition import CAGate,CAValues
from plonk_core.src.proof_system.mod import CustomEvaluations

# Computes the first lagrange polynomial with the given `scale` over `domain`.
def compute_first_lagrange_poly_scaled(domain: Radix2EvaluationDomain,scale:field):
    x_evals = [field.zero(scale.params) for _ in range(domain.size)]
    x_evals[0] = scale
    x_coeffs = INTT(domain,x_evals)
    result_poly = from_coeff_vec(x_coeffs)
    return result_poly

def compute_gate_constraint_satisfiability(domain, 
    range_challenge, logic_challenge, fixed_base_challenge,
    var_base_challenge, prover_key, wl_eval_8n, wr_eval_8n, 
    wo_eval_8n, w4_eval_8n, pi_poly):

    #get Fr
    params = range_challenge.params
    domain_8n = Radix2EvaluationDomain.new(8 * domain.size,params)

    pi_eval_8n = coset_NTT(pi_poly,domain_8n,params)

    gate_contributions = []

    for i in range(domain_8n.size):
        wit_vals = WitnessValues(
            a_val=wl_eval_8n[i],
            b_val=wr_eval_8n[i],
            c_val=wo_eval_8n[i],
            d_val=w4_eval_8n[i]
        )

        custom_vals = CustomEvaluations(
            vals=[
                ("a_next_eval", wl_eval_8n[i + 8]),
                ("b_next_eval", wr_eval_8n[i + 8]),
                ("d_next_eval", w4_eval_8n[i + 8]),
                ("q_l_eval", prover_key.arithmetic.q_l[1][i]),
                ("q_r_eval", prover_key.arithmetic.q_r[1][i]),
                ("q_c_eval", prover_key.arithmetic.q_c[1][i]),
                # Possibly unnecessary but included nonetheless...
                ("q_hl_eval", prover_key.arithmetic.q_hl[1][i]),
                ("q_hr_eval", prover_key.arithmetic.q_hr[1][i]),
                ("q_h4_eval", prover_key.arithmetic.q_hr[1][i])
            ]
        )

        arithmetic = prover_key.arithmetic.compute_quotient_i(i, wit_vals)
        range_term = RangeGate.quotient_term(
            prover_key.range_selector[1][i],
            range_challenge,
            wit_vals,
            custom_vals = RangeValues.from_evaluations(custom_vals)
        )
        logic_term = LogicGate.quotient_term(
            prover_key.logic_selector[1][i],
            logic_challenge,
            wit_vals,
            LogicValues.from_evaluations(custom_vals)
        )
        fixed_base_scalar_mul_term = FBSMGate.quotient_term(
            prover_key.fixed_group_add_selector[1][i],
            fixed_base_challenge,
            wit_vals,
            FBSMValues.from_evaluations(custom_vals)
        )
        curve_addition_term = CAGate.quotient_term(
            prover_key.variable_group_add_selector[1][i],
            var_base_challenge,
            wit_vals,
            CAValues.from_evaluations(custom_vals)
        )

        mid1 = arithmetic.add(pi_eval_8n[i])
        mid2 = mid1.add(range_term)
        mid3 = mid2.add(logic_term)
        mid4 = mid3.add(fixed_base_scalar_mul_term)
        gate_i = mid4.add(curve_addition_term)
        gate_contributions.append(gate_i)

    return gate_contributions

def compute_permutation_checks(
    domain:Radix2EvaluationDomain,
    prover_key,
    wl_eval_8n: list[field], wr_eval_8n: list[field],
    wo_eval_8n: list[field], w4_eval_8n: list[field],
    z_eval_8n: list[field], alpha: field, beta: field, gamma: field):

    #get Fr
    params = alpha.params
    #get NTT domain
    domain_8n:Radix2EvaluationDomain = Radix2EvaluationDomain.new(8 * domain.size,params)

    # Calculate l1_poly_alpha and l1_alpha_sq_evals
    alpha2 = alpha.square()
    l1_poly_alpha = compute_first_lagrange_poly_scaled(domain, alpha2)
    l1_alpha_sq_evals = coset_NTT(l1_poly_alpha, domain_8n,params)

    # Initialize result list
    result = []

    # Calculate permutation contribution for each index
    for i in range(domain_8n.size):
        quotient_i = prover_key.permutation.compute_quotient_i(
            i,
            wl_eval_8n[i],
            wr_eval_8n[i],
            wo_eval_8n[i],
            w4_eval_8n[i],
            z_eval_8n[i],
            z_eval_8n[i + 8],
            alpha,
            l1_alpha_sq_evals[i],
            beta,
            gamma
        )
        result.append(quotient_i)

    return result

def compute(domain: Radix2EvaluationDomain, 
            prover_key, 
            z_poly, z2_poly, 
            w_l_poly, w_r_poly, w_o_poly, w_4_poly, 
            public_inputs_poly, 
            f_poly, table_poly, h1_poly, h2_poly, 
            alpha:field, beta, gamma, delta, epsilon, zeta, 
            range_challenge, logic_challenge, 
            fixed_base_challenge, var_base_challenge, 
            lookup_challenge):
    
    #get Fr
    params = alpha.params
    #get NTT domain
    domain_8n = Radix2EvaluationDomain.new(8 * domain.size,params)
    
    l1_poly = compute_first_lagrange_poly_scaled(domain, alpha.one())
    l1_eval_8n = coset_NTT(l1_poly,domain_8n,params)

    z_eval_8n = coset_NTT(z_poly,domain_8n,params)
    z_eval_8n += z_eval_8n[:8]

    wl_eval_8n = coset_NTT(w_l_poly,domain_8n,params)
    wl_eval_8n += wl_eval_8n[:8]

    wr_eval_8n = coset_NTT(w_r_poly,domain_8n,params)
    wr_eval_8n += wr_eval_8n[:8]

    wo_eval_8n = coset_NTT(w_o_poly,domain_8n,params)

    w4_eval_8n = coset_NTT(w_4_poly,domain_8n,params)
    w4_eval_8n += w4_eval_8n[:8]

    z2_eval_8n = coset_NTT(z2_poly,domain_8n,params)
    z2_eval_8n +=z2_eval_8n[:8]

    f_eval_8n = coset_NTT(f_poly,domain_8n,params)

    table_eval_8n = coset_NTT(table_poly,domain_8n,params)
    table_eval_8n += table_eval_8n[:8]

    h1_eval_8n = coset_NTT(h1_poly,domain_8n,params)
    h1_eval_8n += h1_eval_8n[:8]

    h2_eval_8n = coset_NTT(h2_poly,domain_8n,params)

    gate_constraints = compute_gate_constraint_satisfiability(
        domain,
        range_challenge,logic_challenge,
        fixed_base_challenge,var_base_challenge,
        prover_key,
        wl_eval_8n,wr_eval_8n,wo_eval_8n,w4_eval_8n,
        public_inputs_poly,
    )

    permutation = compute_permutation_checks(
        domain,
        prover_key,
        wl_eval_8n,wr_eval_8n,wo_eval_8n,w4_eval_8n,z_eval_8n,
        alpha,beta,gamma,
    )

    lookup = prover_key.lookup.compute_lookup_quotient_term(
        domain,
        wl_eval_8n,
        wr_eval_8n,
        wo_eval_8n,
        w4_eval_8n,
        f_eval_8n,
        table_eval_8n,
        h1_eval_8n,
        h2_eval_8n,
        z2_eval_8n,
        l1_eval_8n,
        delta,
        epsilon,
        zeta,
        lookup_challenge,
    )
    quotient = []
    for i in range(domain_8n.size):
        numerator = gate_constraints[i].add(permutation[i])
        numerator = numerator.add(lookup[i])
        denominator = field.inverse(prover_key.v_h_coset_8n[i],params)
        res = numerator.mul(denominator)
        quotient.append(res)

    quotient_poly = coset_INTT(quotient,domain_8n)
    hx = from_coeff_vec(quotient_poly)

    return hx