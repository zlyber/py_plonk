from dataclasses import dataclass
from bls12_381 import fr
from typing import List, Tuple
from domain import Radix2EvaluationDomain
from arithmetic import poly_mul_const,poly_add_poly
from plonk_core.src.permutation.constants import K1,K2,K3
@dataclass
class Permutation:
    # Left Permutation
    left_sigma: Tuple[List[fr.Fr],List[fr.Fr]]

    # Right Permutation
    right_sigma: Tuple[List[fr.Fr],List[fr.Fr]]

    # Output Permutation
    out_sigma: Tuple[List[fr.Fr],List[fr.Fr]]

    # Fourth Permutation
    fourth_sigma: Tuple[List[fr.Fr],List[fr.Fr]]

    # Linear Evaluations
    linear_evaluations: List[fr.Fr]

    def compute_quotient_i(self, index,
        w_l_i: fr.Fr, w_r_i: fr.Fr, w_o_i: fr.Fr, w_4_i: fr.Fr,
        z_i: fr.Fr, z_i_next: fr.Fr,
        alpha: fr.Fr, l1_alpha_sq: fr.Fr,
        beta: fr.Fr, gamma: fr.Fr):

        a = self.compute_quotient_identity_range_check_i(
            index, w_l_i, w_r_i, w_o_i, w_4_i, z_i, alpha, beta, gamma,
        )
        b = self.compute_quotient_copy_range_check_i(
            index, w_l_i, w_r_i, w_o_i, w_4_i, z_i_next, alpha, beta, gamma,
        )
        c = self.compute_quotient_term_check_one_i(z_i, l1_alpha_sq)

        res = a.add(b)
        res = res.add(c)
        return res
    
    # Computes the following:
    # (a(x) + beta * X + gamma) (b(X) + beta * k1 * X + gamma) (c(X) + beta *
    # k2 * X + gamma)(d(X) + beta * k3 * X + gamma)z(X) * alpha
    def compute_quotient_identity_range_check_i(
        self,index,
        w_l_i: fr.Fr,w_r_i: fr.Fr,w_o_i: fr.Fr,w_4_i: fr.Fr,
        z_i: fr.Fr,alpha: fr.Fr,beta: fr.Fr,gamma: fr.Fr,):

        x = self.linear_evaluations[index]

        k1 = K1()
        k2 = K2()
        k3 = K3()
        
        mid1_1 = beta.mul(x)
        mid1_2 = w_l_i.add(mid1_1)
        mid1 = mid1_2.add(gamma)

        mid2_1_1 = beta.mul(k1)
        mid2_1 = mid2_1_1.mul(x)
        mid2_2 = w_r_i.add(mid2_1)
        mid2 = mid2_2.add(gamma)

        mid3_1_1 = beta.mul(k2)
        mid3_1 = mid3_1_1.mul(x)
        mid3_2 = w_o_i.add(mid3_1)
        mid3 = mid3_2.add(gamma)

        mid4_1_1 = beta.mul(k3)
        mid4_1 = mid4_1_1.mul(x)
        mid4_2 = w_4_i.add(mid4_1)
        mid4 = mid4_2.add(gamma)

        mid5 = mid1.mul(mid2)
        mid6 = mid5.mul(mid3)
        mid7 = mid6.mul(mid4)
        res = mid7.mul(z_i)
        res = res.mul(alpha)
        return res

    # Computes the following:
    # (a(x) + beta* Sigma1(X) + gamma) (b(X) + beta * Sigma2(X) + gamma) (c(X)
    # + beta * Sigma3(X) + gamma)(d(X) + beta * Sigma4(X) + gamma) Z(X.omega) *
    # alpha
    def compute_quotient_copy_range_check_i(
        self,index,
        w_l_i: fr.Fr,
        w_r_i: fr.Fr,
        w_o_i: fr.Fr,
        w_4_i: fr.Fr,
        z_i_next: fr.Fr,
        alpha: fr.Fr,
        beta: fr.Fr,
        gamma: fr.Fr
    ):
        left_sigma_eval = self.left_sigma[1][index]
        right_sigma_eval = self.right_sigma[1][index]
        out_sigma_eval = self.out_sigma[1][index]
        fourth_sigma_eval = self.fourth_sigma[1][index]

        mid1_1 = beta.mul(left_sigma_eval)
        mid1_2 = w_l_i.add(mid1_1)
        mid1 = mid1_2.add(gamma)

        mid2_1 = beta.mul(right_sigma_eval)
        mid2_2 = w_r_i.add(mid2_1)
        mid2 = mid2_2.add(gamma)
        
        mid3_1 = beta.mul(out_sigma_eval)
        mid3_2 = w_o_i.add(mid3_1)
        mid3 = mid3_2.add(gamma)

        mid4_1 = beta.mul(fourth_sigma_eval)
        mid4_2 = w_4_i.add(mid4_1)
        mid4 = mid4_2.add(gamma)

        mid5 = mid1.mul(mid2)
        mid5 = mid5.mul(mid3)
        mid5 = mid5.mul(mid4)
        mid5 = mid5.mul(z_i_next)
        product = mid5.mul(alpha)

        res = product.neg()
        return res

    # Computes the following:
    # L_1(X)[Z(X) - 1]
    def compute_quotient_term_check_one_i(self, z_i: fr.Fr, l1_alpha_sq: fr.Fr):
        one = z_i.one()
        z_i_sub_one = z_i.sub(one)
        res = z_i_sub_one.mul(l1_alpha_sq)
        return res
    
    # Computes the permutation term of the linearisation polynomial.
    def compute_linearisation(
        self, 
        n, 
        z_challenge: fr.Fr, 
        challengTuple: Tuple[fr.Fr,fr.Fr,fr.Fr], 
        wireTuple: Tuple[fr.Fr,fr.Fr,fr.Fr,fr.Fr], 
        sigmaTuple: Tuple[fr.Fr,fr.Fr,fr.Fr], 
        z_eval, z_poly):
        a = self.compute_lineariser_identity_range_check(
            wireTuple[0],wireTuple[1],wireTuple[2],wireTuple[3],
            z_challenge,
            challengTuple[0],challengTuple[1],challengTuple[2],
            z_poly
        )
        
        b = self.compute_lineariser_copy_range_check(
            wireTuple[0], wireTuple[1], wireTuple[2],
            z_eval,
            sigmaTuple[0],sigmaTuple[1],sigmaTuple[2],
            challengTuple[0],challengTuple[1],challengTuple[2],
            self.fourth_sigma[0]
        )

        domain = Radix2EvaluationDomain.new(n,z_challenge)
        alpha2 = challengTuple[0].square()
        c = self.compute_lineariser_check_is_one(
            domain,
            z_challenge,
            alpha2,
            z_poly
        )
        ab = poly_add_poly(a,b)
        abc = poly_add_poly(ab,c)
        return abc
    
    # Computes the following:
    # -(a_eval + beta * sigma_1 + gamma)(b_eval + beta * sigma_2 + gamma)
    # (c_eval + beta * sigma_3 + gamma) * beta *z_eval * alpha^2 * Sigma_4(X)
    def compute_lineariser_copy_range_check(
        self,
        a_eval: fr.Fr, b_eval: fr.Fr, c_eval: fr.Fr,
        z_eval: fr.Fr,
        sigma_1_eval: fr.Fr,
        sigma_2_eval: fr.Fr,
        sigma_3_eval: fr.Fr,
        alpha: fr.Fr, beta: fr.Fr, gamma: fr.Fr,
        fourth_sigma_poly: List[fr.Fr],
    ):
        # a_eval + beta * sigma_1 + gamma
        beta_sigma_1 = beta.mul(sigma_1_eval)
        a_0 = a_eval.add(beta_sigma_1)
        a_0 = a_0.add(gamma)

        # b_eval + beta * sigma_2 + gamma
        beta_sigma_2 = beta.mul(sigma_2_eval)
        a_1 = b_eval.add(beta_sigma_2)
        a_1 = a_1.add(gamma)

        # c_eval + beta * sigma_3 + gamma
        beta_sigma_3 = beta.mul(sigma_3_eval)
        a_2 = c_eval.add(beta_sigma_3)
        a_2 = a_2.add(gamma)

        beta_z_eval = beta.mul(z_eval)
        a = a_0.mul(a_1)
        a = a.mul(a_2)
        a = a.mul(beta_z_eval)
        a = a.mul(alpha)
        neg_a = a.neg()

        res = poly_mul_const(fourth_sigma_poly,neg_a)
        return res
    
    # Computes the following:
    # (a_eval + beta * z_challenge + gamma)(b_eval + beta * K1 * z_challenge +
    # gamma)(c_eval + beta * K2 * z_challenge + gamma) * alpha z(X)
    def compute_lineariser_identity_range_check(
        self,
        a_eval: fr.Fr, b_eval: fr.Fr, c_eval: fr.Fr, d_eval: fr.Fr,
        z_challenge: fr.Fr,
        alpha: fr.Fr, beta: fr.Fr, gamma: fr.Fr,
        z_poly: List[fr.Fr]
    ):
        beta_z = beta.mul(z_challenge)

        # a_eval + beta * z_challenge + gamma
        a_0 = a_eval.add(beta_z)
        a_0 = a_0.add(gamma)

        # b_eval + beta * K1 * z_challenge + gamma
        beta_z_k1 = K1().mul(beta_z)
        a_1 = b_eval.add(beta_z_k1)
        a_1 = a_1.add(gamma)

        # c_eval + beta * K2 * z_challenge + gamma
        beta_z_k2 = K2().mul(beta_z)
        a_2 = c_eval.add(beta_z_k2)
        a_2 = a_2.add(gamma)

        # d_eval + beta * K3 * z_challenge + gamma
        beta_z_k3 = K3().mul(beta_z)
        a_3 = d_eval.add(beta_z_k3)
        a_3 = a_3.add(gamma)

        a = a_0.mul(a_1)
        a = a.mul(a_2)
        a = a.mul(a_3)
        a = a.mul(alpha)
        res = poly_mul_const(z_poly,a)
        return res
    
    def compute_lineariser_check_is_one(
        self, 
        domain: Radix2EvaluationDomain, 
        z_challenge: fr.Fr, 
        alpha_sq: fr.Fr, 
        z_coeffs: List[fr.Fr]):

        lagrange_coefficients = domain.evaluate_all_lagrange_coefficients(z_challenge)
        l_1_z = lagrange_coefficients[0]
        const = l_1_z.mul(alpha_sq)
        res = poly_mul_const(z_coeffs,const)
        return res
