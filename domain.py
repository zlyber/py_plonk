from dataclasses import dataclass
import gmpy2
from bls12_381 import fr
import math

@dataclass
class Radix2EvaluationDomain:
    size: int
    log_size_of_group: int
    size_as_field_element: fr.Fr
    size_inv: fr.Fr
    group_gen: fr.Fr
    group_gen_inv: fr.Fr
    generator_inv: fr.Fr

    @classmethod
    def new(cls, num_coeffs: int, params:fr.Fr):
        # Compute the size of our evaluation domain
        size = num_coeffs if num_coeffs & (num_coeffs - 1) == 0 else 2 ** num_coeffs.bit_length()
        log_size_of_group = size.bit_length()-1
        
        # Check if log_size_of_group exceeds TWO_ADICITY
        if log_size_of_group > params.TWO_ADICITY:
            return None

        # Compute the generator for the multiplicative subgroup.
        # It should be the 2^(log_size_of_group) root of unity.
        group_gen = fr.Fr(value=params.get_root_of_unity(size))
        
        # Check that it is indeed the 2^(log_size_of_group) root of unity.
        group_gen_pow = group_gen.pow(size)
        assert group_gen_pow == params.one()

        size_as_field_element=fr.Fr.from_repr(size)
        size_inv = fr.Fr.inverse(size_as_field_element)
        group_gen_inv = fr.Fr.inverse(group_gen)
        generator_inv = fr.Fr.inverse(params.multiplicative_generator())

        return cls(size, log_size_of_group, size_as_field_element, size_inv, group_gen, group_gen_inv, generator_inv)
    
    # Evaluate all Lagrange polynomials at tau to get the lagrange coefficients.
    # Define the following as
    # - H: The coset we are in, with generator g and offset h
    # - m: The size of the coset H
    # - Z_H: The vanishing polynomial for H. Z_H(x) = prod_{i in m} (x - hg^i) = x^m - h^m
    # - v_i: A sequence of values, where v_0 = 1/(m * h^(m-1)), and v_{i + 1} = g * v_i
    #
    # We then compute L_{i,H}(tau) as `L_{i,H}(tau) = Z_H(tau) * v_i / (tau - h g^i)`
    #
    # However, if tau in H, both the numerator and denominator equal 0
    # when i corresponds to the value tau equals, and the coefficient is 0 everywhere else.
    # We handle this case separately, and we can easily detect by checking if the vanishing poly is 0.
    def evaluate_all_lagrange_coefficients(self, tau: fr.Fr):
        size = self.size
        t_size = tau.pow(size)
        domain_offset = tau.one()
        one = tau.one()
        z_h_at_tau = t_size.sub(domain_offset)
        zero = fr.Fr.zero()

        if z_h_at_tau.value == 0:
            u = [zero for _ in range(size)]
            omega_i = domain_offset
            for i in range(size):
                if omega_i.value == tau.value:
                    u[i] = one
                    break
                omega_i = omega_i.mul(self.group_gen)
            return u
        else:
            # In this case we have to compute `Z_H(tau) * v_i / (tau - h g^i)`
            # for i in 0..size
            # We actually compute this by computing (Z_H(tau) * v_i)^{-1} * (tau - h g^i)
            # and then batch inverting to get the correct lagrange coefficients.
            # We let `l_i = (Z_H(tau) * v_i)^-1` and `r_i = tau - h g^i`
            # Notice that since Z_H(tau) is i-independent,
            # and v_i = g * v_{i-1}, it follows that
            # l_i = g^-1 * l_{i-1}
            # TODO: consider caching the computation of l_i to save N multiplications
            from arithmetic import batch_inversion

            # v_0_inv = m * h^(m-1)
            f_size = fr.Fr.from_repr(size)
            pow_dof = domain_offset.pow(size-1)
            v_0_inv = f_size.mul(pow_dof)
            z_h_at_tau_inv = fr.Fr.inverse(z_h_at_tau)
            l_i = z_h_at_tau_inv.mul(v_0_inv)
            negative_cur_elem = domain_offset.neg()
            lagrange_coefficients_inverse = [zero for _ in range(size)]
            for i in range(size):
                r_i = tau.add(negative_cur_elem)
                lagrange_coefficients_inverse[i] = l_i.mul(r_i)
                # Increment l_i and negative_cur_elem
                l_i = l_i.mul(self.group_gen_inv)
                negative_cur_elem = negative_cur_elem.mul(self.group_gen)
            batch_inversion(lagrange_coefficients_inverse)
            return lagrange_coefficients_inverse
        
    # def evaluate_all_lagrange_coefficients(self, tau: field):
    #     size = self.size
    #     t_size = tau.pow(size)
    #     one = tau.one()
    #     zero = field.zero(tau.params)

    #     if t_size.value == one.value:
    #         u = [zero for _ in range(size)]
    #         omega_i = one
    #         for i in range(size):
    #             if omega_i.value == tau.value:
    #                 u[i] = one
    #                 break
    #             omega_i =omega_i.mul(self.group_gen)
    #         return u
    #     else:
    #         from arithmetic import batch_inversion
    #         t_size_sub_one = t_size.sub(one)
    #         l = t_size_sub_one.mul(self.size_inv)
    #         r = one
    #         u = [zero] * size
    #         ls = [zero] * size
    #         for i in range(size):
    #             u[i] = tau.sub(r)
    #             ls[i] = l
    #             l = l.mul(self.group_gen)
    #             r = r.mul(self.group_gen)

    #         batch_inversion(u)

    #         for i in range(size):
    #             u[i] = ls[i].mul(u[i])

    #         return u
    
    # This evaluates the vanishing polynomial for this domain at tau.
    # For multiplicative subgroups, this polynomial is `z(X) = X^self.size - 1`.
    def evaluate_vanishing_polynomial(self, tau: fr.Fr):
        one = tau.one()
        pow_tau = tau.pow(self.size)
        return pow_tau.sub(one)
    
    # Returns the `i`-th element of the domain, where elements are ordered by
    # their power of the generator which they correspond to.
    # e.g. the `i`-th element is g^i
    def element(self, i):
        # TODO: Consider precomputed exponentiation tables if we need this to be faster.
        res = self.group_gen.pow(i)
        return res
    
        