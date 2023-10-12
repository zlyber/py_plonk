from field import field
from domain import Radix2EvaluationDomain
from arithmetic import INTT,from_coeff_vec

def as_evals(public_inputs,pi_pos,n,params):
    pi = [field.zero(params) for _ in range(n)]
    for pos in pi_pos:
        pi[pos] = public_inputs
    return pi

def into_dense_poly(public_inputs,pi_pos,n,params):
    domain = Radix2EvaluationDomain.new(n,params)
    evals = as_evals(public_inputs,pi_pos,n,params)
    pi_coeffs = INTT(domain,evals)
    pi_poly = from_coeff_vec(pi_coeffs)
    return pi_poly