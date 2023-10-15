from dataclasses import dataclass
from field import field
import gmpy2
from bls12_381 import fq,fr
import math

@dataclass
class Radix2EvaluationDomain:
    size: int
    log_size_of_group: int
    size_as_field_element: field
    size_inv: field
    group_gen: field
    group_gen_inv: field
    generator_inv: field

    @classmethod
    def new(cls, num_coeffs: int, params):
        # Compute the size of our evaluation domain
        size = num_coeffs if num_coeffs & (num_coeffs - 1) == 0 else 2 ** num_coeffs.bit_length()
        log_size_of_group = size.bit_length()-1
        
        # Check if log_size_of_group exceeds TWO_ADICITY
        if log_size_of_group > params.TWO_ADICITY:
            return None

        # Compute the generator for the multiplicative subgroup.
        # It should be the 2^(log_size_of_group) root of unity.
        group_gen = field(params.get_root_of_unity(size),params)
        
        # Check that it is indeed the 2^(log_size_of_group) root of unity.
        group_gen_pow = group_gen.pow(size)
        assert group_gen_pow.value == params.one()

        size_as_field_element=field.from_repr(size,params)
        size_inv = field.inverse(size_as_field_element,params)
        group_gen_inv = field.inverse(group_gen,params)
        generator_inv = field.inverse(params.multiplicative_generator(),params)

        return cls(size, log_size_of_group, size_as_field_element, size_inv, group_gen, group_gen_inv, generator_inv)
    
    def evaluate_all_lagrange_coefficients(self, tau: field):
        size = self.size
        t_size = tau.pow(size)
        one = tau.one()
        zero = field.zero(tau.params)

        if t_size.value == one.value:
            u = [zero] * size
            omega_i = one
            for i in range(size):
                if omega_i.value == tau.value:
                    u[i] = one
                    break
                omega_i =omega_i.mul(self.group_gen)
            return u
        else:
            from arithmetic import batch_inversion
            t_size_sub_one = t_size.sub(one)
            l = t_size_sub_one.mul(self.size_inv)
            r = one
            u = [zero] * size
            ls = [zero] * size
            for i in range(size):
                u[i] = tau.sub(r)
                ls[i] = l
                l = l.mul(self.group_gen)
                r = r.mul(self.group_gen)

            batch_inversion(u)

            for i in range(size):
                u[i] = ls[i].mul(u[i])

            return u
    
    # This evaluates the vanishing polynomial for this domain at tau.
    # For multiplicative subgroups, this polynomial is `z(X) = X^self.size - 1`.
    def evaluate_vanishing_polynomial(self, tau: field):
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
    
        