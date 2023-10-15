from plonk_core.src.permutation import constants
from field import field
from arithmetic import NTT,INTT,from_coeff_vec
import copy
import math

def numerator_irreducible(root:field, w:field, k:field, beta:field, gamma:field):
    mid1 = beta.mul(k)
    mid2 = mid1.mul(root)
    mid3 = w.add(mid2)
    mid4 = mid3.add(gamma)
    return mid4

def denominator_irreducible(w:field, sigma:field, beta:field, gamma:field):
    mid1 = beta.mul(sigma)
    mid2 = w.add(mid1)
    mid3 = mid2.add(gamma)
    return mid3

def lookup_ratio(delta:field, epsilon:field, f:field, t:field, t_next:field,
                h_1:field, h_1_next:field, h_2:field):
    one = delta.one()
    one_plus_delta =delta.add(one)
    epsilon_one_plus_delta = epsilon.mul(one_plus_delta)

    mid1 = epsilon.add(f)
    mid2 = epsilon_one_plus_delta.add(t)
    mid3 = delta.mul(t_next)
    mid4 = mid2.add(mid3)
    mid5 = one_plus_delta.mul(mid1)
    result = mid5.mul(mid4)
    mid6 = h_2.add(delta)
    mid7 = result.add(h_1)
    mid8 = mid7.add(mid6)
    result = result.mul(mid8)
    mid9 = epsilon_one_plus_delta.add(h_2)
    mid10 = h_1_next.mul(delta)
    mid11 = mid9.add(mid10)
    result = result.mul(mid11)

    result_inv = field.inverse(result,result.params)

    return result_inv


def compute_permutation_poly(domain, wires, beta:field, gamma, sigma_polys):
    n = domain.size
    #get Fr
    params = beta.params

    # Constants defining cosets H, k1H, k2H, etc
    ks = [field.zero(params),constants.K1(params),constants.K2(params),constants.K3(params)]
    sigma_mappings = [[],[],[],[]]

    sigma_mappings[0] = NTT(domain,sigma_polys[0],params)
    sigma_mappings[1] = NTT(domain,sigma_polys[1],params)
    sigma_mappings[2] = NTT(domain,sigma_polys[2],params)
    sigma_mappings[3] = NTT(domain,sigma_polys[3],params)

    # Transpose wires and sigma values to get "rows" in the form [wl_i,
    # wr_i, wo_i, ... ] where each row contains the wire and sigma
    # values for a single gate
    gatewise_wires = [
        [w0, w1, w2, w3] for w0, w1, w2, w3 in zip(wires[0],wires[1],wires[2],wires[3])
    ]
    gatewise_sigmas = [
        [s0, s1, s2, s3] for s0, s1, s2, s3 in zip(sigma_mappings[0],sigma_mappings[1],
                                                   sigma_mappings[2],sigma_mappings[3])
    ]

    # Compute all roots, same as calculating twiddles, but doubled in size
    log_size = int(math.log2(n))
    roots = [field.zero(params)] * (1 << log_size )
    roots[0] = beta.one()
    for idx in range(1, len(roots)):
        roots[idx] = roots[idx - 1].mul(domain.group_gen)
    
    # Initialize an empty list for product_argument
    product_argument = []

    # Associate each wire value in a gate with the k defining its coset
    for gate_root, gate_sigmas, gate_wires in zip(roots, gatewise_sigmas, gatewise_wires):
        # Initialize numerator and denominator products
        numerator_product = beta.one()
        denominator_product = beta.one()

        # Now the ith element represents gate i and will have the form:
        # (root_i, ((w0_i, s0_i, k0), (w1_i, s1_i, k1), ..., (wm_i, sm_i,
        # km)))   for m different wires, which is all the
        # information   needed for a single product coefficient
        # for a single gate Multiply up the numerator and
        # denominator irreducibles for each gate and pair the results
        for sigma, wire, k in zip(gate_sigmas, gate_wires, ks):

            # Calculate numerator and denominator for each wire
            numerator_temp = numerator_irreducible(gate_root, wire, k, beta, gamma)
            numerator_product= numerator_product.mul(numerator_temp)
            denominator_temp = denominator_irreducible(wire, sigma, beta, gamma)
            denominator_product = denominator_product.mul(denominator_temp)
        
        # Calculate the product coefficient for the gate
        denominator_product_under = field.inverse(denominator_product,params)
        gate_coefficient = numerator_product.mul(denominator_product_under)
        
        # Append the gate coefficient to the product_argument list
        product_argument.append(gate_coefficient)

    z=[]
    # First element is one
    state = beta.one()
    z.append(state)

    # Accumulate by successively multiplying the scalars        
    for s in product_argument:
        state = state.mul(s)
        z.append(state)

    # Remove the last(n+1'th) element
    z.pop()
    
    #Compute z poly
    INTT(domain,z)
    z_poly = from_coeff_vec(z)
    
    return z_poly

# Define a Python function that mirrors the Rust function
def compute_lookup_permutation_poly(domain, f, t, h_1, h_2, delta, epsilon):
    n = domain.size

    #get Fr
    params = delta.params

    assert len(f) == n
    assert len(t) == n
    assert len(h_1) == n
    assert len(h_2) == n

    t_next = t[1:] + [t[0]]
    h_1_next = h_1[1:] + [h_1[0]]

    product_arguments = []
    for f_val, t_val, t_next_val, h_1_val, h_1_next_val, h_2_val in zip(f, t, t_next, h_1, h_1_next, h_2):
        portion = lookup_ratio(delta, epsilon, f_val, t_val, t_next_val, h_1_val, h_1_next_val, h_2_val)
        product_arguments.append(portion)

    state = delta.one()
    p = [state]
    for s in product_arguments:
        state = state.mul(s)
        p.append(state)
    
    p.pop()
    INTT(domain,p)
    p_poly = from_coeff_vec(p)
    
    return p_poly

