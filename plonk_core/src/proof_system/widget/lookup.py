from field import field
from domain import Radix2EvaluationDomain
from dataclasses import dataclass
from typing import List, Tuple
from plonk_core.lookup.multiset import MultiSet
from plonk_core.src.utils import lc
from arithmetic import poly_add_poly,poly_mul_const
@dataclass
class Lookup:
    # Lookup selector
    q_lookup: Tuple[List[field],List[field]]
    # Column 1 of lookup table
    table_1: List[field]
    # Column 2 of lookup table
    table_2: List[field]
    # Column 3 of lookup table
    table_3: List[field]
    # Column 4 of lookup table
    table_4: List[field]

    # Compute lookup portion of quotient polynomial
    def compute_lookup_quotient_term(self,
        domain: Radix2EvaluationDomain,
        wl_eval_8n: List[field],
        wr_eval_8n: List[field],
        wo_eval_8n: List[field],
        w4_eval_8n: List[field],
        f_eval_8n: List[field],
        table_eval_8n: List[field],
        h1_eval_8n: List[field],
        h2_eval_8n: List[field],
        z2_eval_8n: List[field],
        l1_eval_8n: List[field],
        delta: field,
        epsilon: field,
        zeta: field,
        lookup_sep:field):

        params = zeta.params
        domain_8n:Radix2EvaluationDomain = Radix2EvaluationDomain.new(8 * domain.size,params)

        # Initialize result list
        result = []

        # Calculate lookup quotient term for each index
        for i in range(domain_8n.size):
            quotient_i = self.compute_quotient_i(
                i,
                wl_eval_8n[i],
                wr_eval_8n[i],
                wo_eval_8n[i],
                w4_eval_8n[i],
                f_eval_8n[i],
                table_eval_8n[i],
                table_eval_8n[i + 8],
                h1_eval_8n[i],
                h1_eval_8n[i + 8],
                h2_eval_8n[i],
                z2_eval_8n[i],
                z2_eval_8n[i + 8],
                l1_eval_8n[i],
                delta,
                epsilon,
                zeta,
                lookup_sep
            )
            result.append(quotient_i)

        return result
    
    def compute_quotient_i(
        self,
        index: field,
        w_l_i: field,
        w_r_i: field,
        w_o_i: field,
        w_4_i: field,
        f_i: field,
        table_i: field,
        table_i_next: field,
        h1_i: field,
        h1_i_next: field,
        h2_i: field,
        z2_i: field,
        z2_i_next: field,
        l1_i: field,
        delta: field,
        epsilon: field,
        zeta: field,
        lookup_sep: field
    ):
        # q_lookup(X) * (a(X) + zeta * b(X) + (zeta^2 * c(X)) + (zeta^3 * d(X)
        # - f(X))) * α_1
        one = delta.one()
        lookup_sep_sq = lookup_sep.square()
        lookup_sep_cu = lookup_sep_sq.mul(lookup_sep)
        one_plus_delta = delta.add(one)
        epsilon_one_plus_delta = epsilon.mul(one_plus_delta)

        q_lookup_i = self.q_lookup[1][index]
        compressed_tuple = lc([w_l_i, w_r_i, w_o_i, w_4_i], zeta)
        mid1 = compressed_tuple.sub(f_i)
        mid2 = q_lookup_i.mul(mid1)
        a = mid2.mul(lookup_sep)

        # z2(X) * (1+δ) * (ε+f(X)) * (ε*(1+δ) + t(X) + δt(Xω)) * lookup_sep^2
        b_0 = epsilon.add(f_i)
        b_1_1 = epsilon_one_plus_delta.add(table_i)
        b_1_2 = delta.mul(table_i_next)
        b_1 = b_1_1.add(b_1_2)
        mid1 = z2_i.mul(one_plus_delta)
        mid2 = mid1.mul(b_0)
        mid3 = mid2.mul(b_1)
        b = mid3.mul(lookup_sep_sq)

        # − z2(Xω) * (ε*(1+δ) + h1(X) + δ*h2(X)) * (ε*(1+δ) + h2(X) + δ*h1(Xω))
        # * lookup_sep^2
        c_0_1 = epsilon_one_plus_delta.add(h1_i)
        c_0_2 = delta.mul(h2_i)
        c_0 =c_0_1.add(c_0_2)
        c_1_1 = epsilon_one_plus_delta.add(h2_i)
        c_1_2 = delta.mul(h1_i_next)
        c_1 = c_1_1.add(c_1_2)
        neg_z2_next = z2_i_next.neg()
        mid1 = neg_z2_next.mul(c_0)
        mid2 = mid1.mul(c_1)
        c = mid2.mul(lookup_sep_cu)

        d_1 = z2_i.sub(one)
        d_2 = l1_i.mul(lookup_sep_cu)
        d = d_1.mul(d_2)

        mid1 = a.add(b)
        mid2 = mid1.add(c)
        res = mid2.add(d)
        return res
    
    def compute_linearisation(
        self,
        l1_eval: field,
        a_eval: field,
        b_eval: field,
        c_eval: field,
        d_eval: field,
        f_eval: field,
        table_eval: field,
        table_next_eval: field,
        h1_next_eval: field,
        h2_eval: field,
        z2_next_eval: field,
        delta: field,
        epsilon: field,
        zeta: field,
        z2_poly: list[field],
        h1_poly: list[field],
        lookup_sep: field
    ):
        one = delta.one()
        lookup_sep_sq = lookup_sep.mul(lookup_sep)
        lookup_sep_cu = lookup_sep.mul(lookup_sep_sq)
        one_plus_delta = delta.add(one)
        epsilon_one_plus_delta = epsilon.mul(one_plus_delta)

        compressed_tuple = lc([a_eval, b_eval, c_eval, d_eval], zeta)
        compressed_tuple_sub_f_eval = compressed_tuple.sub(f_eval)
        const1 = compressed_tuple_sub_f_eval.mul(lookup_sep)
        a = poly_mul_const(self.q_lookup[0],const1)

        # z2(X) * (1 + δ) * (ε + f_bar) * (ε(1+δ) + t_bar + δ*tω_bar) *
        # lookup_sep^2
        b_0 = epsilon.add(f_eval)
        epsilon_one_plus_delta_plus_tabel_eval = epsilon_one_plus_delta.add(table_eval)
        delta_times_table_next_eval = delta.mul(table_next_eval)
        b_1 = epsilon_one_plus_delta_plus_tabel_eval.add(delta_times_table_next_eval)
        b_2 = l1_eval.mul(lookup_sep_cu)
        one_plus_delta_b_0 = one_plus_delta.mul(b_0)
        one_plus_delta_b_0_b_1 = one_plus_delta_b_0.mul(b_1)
        one_plus_delta_b_0_b_1_lookup = one_plus_delta_b_0_b_1.mul(lookup_sep_sq)
        const2 = one_plus_delta_b_0_b_1_lookup.add(b_2)
        b = poly_mul_const(z2_poly,const2)

        # h1(X) * (−z2ω_bar) * (ε(1+δ) + h2_bar  + δh1ω_bar) * lookup_sep^2
        neg_z2_next_eval = z2_next_eval.neg()
        c_0 = neg_z2_next_eval.mul(lookup_sep_sq)
        epsilon_one_plus_delta_h2_eval = epsilon_one_plus_delta.add(h2_eval)
        delta_h1_next_eval =  delta.add(h1_next_eval)
        c_1 = epsilon_one_plus_delta_h2_eval.add(delta_h1_next_eval)
        c0_c1 = c_0.mul(c_1)
        c = poly_mul_const(h1_poly, c0_c1)

        ab = poly_add_poly(a, b)
        abc = poly_add_poly(ab, c)
        return abc
