import gmpy2
import copy
from bls12_381 import fr,fq
from domain import Radix2EvaluationDomain
from structure import AffinePointG1
from jacobian import ProjectivePointG1
import math
import random

def reverse_bits(operand, bit_count):
    acc = 0
    for i in range(bit_count):
        acc = (acc << 1) | ((operand >> i) & 1)
    return acc

def derange(xi, log_len):
    for idx in range(1, len(xi) - 1):
        ridx = reverse_bits(idx, log_len)
        if idx < ridx:
            xi[idx], xi[ridx] = xi[ridx], xi[idx]
    return xi

def precompute_twiddles(domain:Radix2EvaluationDomain, root:fr.Fr):
    log_size = int(math.log2(domain.size))
    powers = [root.zero()] * (1 << (log_size - 1))
    powers[0] = root.one()
    for idx in range(1, len(powers)):
        powers[idx] = powers[idx - 1].mul(root)
    return powers

def operator(domain:Radix2EvaluationDomain, xi:list[fr.Fr], root:fr.Fr):
    log_size = int(math.log2(domain.size))
    xi = derange(xi,log_size)
    twiddles=precompute_twiddles(domain,root)
    chunk = 2
    twiddle_chunk = domain.size // 2
    for i in range(log_size):
        for j in range(0, domain.size, chunk):
            t = xi[j + chunk // 2]       # Copy Right[0]
            xi[j + chunk // 2] = xi[j]   # Right[0] = Left[0]
            xi[j] = xi[j].add(t)
            xi[j + chunk // 2] = xi[j + chunk // 2].sub(t)
            for m in range(chunk // 2 - 1):
                twiddle = twiddles[(m + 1) * twiddle_chunk]
                t1 = xi[j + chunk // 2 + m + 1]
                t1 = t1.mul(twiddle)
                xi[j + chunk // 2 + m + 1] = xi[j + m + 1]
                xi[j + m + 1] = xi[j + m + 1].add(t1)  # a + b * w
                # a - b * w
                xi[j + chunk // 2 + m + 1] = xi[j + chunk // 2 + m + 1].sub(t1)  
        chunk *= 2  # Merge up
        twiddle_chunk //= 2
    return xi

def resize(self, target_len, padding):
    res = self[:]
    if len(self) < target_len:
        num_to_pad = target_len - len(self)
        res.extend([padding for _ in range(num_to_pad)])
    return res

def distribute_powers(coeffs:list[fr.Fr], g):
    g_field = fr.Fr(value = g)
    one = g_field.one()
    distribute_powers_and_mul_by_const(coeffs, g_field, one)

# Multiply the `i`-th element of `coeffs` with `c*g^i`.
def distribute_powers_and_mul_by_const(coeffs, g, c):
    pow = c
    for i in range(len(coeffs)):
        coeffs[i] = coeffs[i].mul(pow)
        pow = pow.mul(g)

def degree(poly):
    if len(poly)==0:
        return 0
    else:
        assert not poly[-1]==0
        return len(poly) - 1

def convert_to_bigints(p: list[fr.Fr]):
    coeffs = [fr.Fr(value = s.into_repr()) for s in p]
    return coeffs

def skip_leading_zeros_and_convert_to_bigints(p: list[fr.Fr]):
    num_leading_zeros = 0
    while num_leading_zeros < len(p) and p[num_leading_zeros] == 0:
        num_leading_zeros += 1

    coeffs = convert_to_bigints(p[num_leading_zeros:])
    return num_leading_zeros, coeffs

def NTT(domain,coeffs):
    zero = fr.Fr.zero()
    #add zero to resize
    resize_coeffs = resize(coeffs,domain.size,zero)
    evals = operator(domain,resize_coeffs,domain.group_gen)
    return evals

def INTT(domain,evals):
    #add zero to resize
    zero = fr.Fr.zero()
    resize_evals = resize(evals,domain.size,zero)
    evals = operator(domain,resize_evals,domain.group_gen_inv)
    for i in range(len(evals)):
        evals[i] = evals[i].mul(domain.size_inv)
    return evals

# Compute a NTT over a coset of the domain, modifying the input vector in place.
def coset_NTT(coeffs:list[fr.Fr], domain):
    modified_coeffs = coeffs[:]
    distribute_powers(modified_coeffs, fr.Fr.GENERATOR)
    evals = NTT(domain,modified_coeffs)
    return evals

# Compute a INTT over a coset of the domain, modifying the input vector in place.
def coset_INTT(evals:list[fr.Fr], domain):
    #add zero to resize
    zero = fr.Fr.zero()
    resize_evals = resize(evals,domain.size,zero)
    evals = operator(domain,resize_evals,domain.group_gen_inv)
    distribute_powers_and_mul_by_const(evals, domain.generator_inv,domain.size_inv)
    return evals

def from_coeff_vec(coeffs:list):
    result = coeffs[:]
    while result and result[-1].value == 0:
        result.pop()
    return result

def poly_add_poly(self: list[fr.Fr], other: list[fr.Fr]):
    if len(self) == 0:
        res = other[:]
        return res
    if len(other) == 0:
        res = self[:]
        return res
    elif len(self) >= len(other):
        result = self[:]
        for i in range(len(other)):
            result[i] = result[i].add(other[i])

        result = from_coeff_vec(result)
        return result
    else:
        result = other[:]
        for i in range(len(self)):
            result[i] = result[i].add(self[i])
    
        result = from_coeff_vec(result)
        return result

def poly_mul_const(poly:list[fr.Fr],elem:fr.Fr):
    if len(poly) == 0 or elem.value == 0:
        return poly
    else:
        result = poly[:]
        for i in range(len(result)):
            result[i] = result[i].mul(elem)
        return result
    
def divide_with_q_and_r(self: list[fr.Fr], divisor: list[fr.Fr]):
    if len(self) == 0:
        return self
    elif len(divisor) == 0:
        raise ValueError("Dividing by zero polynomial")
    elif len(self) < len(divisor):
        zero = fr.Fr.zero()
        return [zero]
    else:
        zero = fr.Fr.zero()
        quotient = [zero for _ in range(len(self) - len(divisor) + 1)]
        remainder = self[:]
        divisor_leading = divisor[-1]
        divisor_leading_inv = fr.Fr.inverse(divisor_leading)
        while len(remainder) != 0 and len(remainder) >= len(divisor):
            remainder_leading = remainder[-1]
            cur_q_coeff = remainder_leading.mul(divisor_leading_inv)
            cur_q_degree = len(remainder) - len(divisor)
            quotient[cur_q_degree] = cur_q_coeff
            for i, div_coeff in enumerate(divisor):
                temp = cur_q_coeff.mul(div_coeff)
                remainder[cur_q_degree + i] = remainder[cur_q_degree + i].sub(temp)
            while remainder[-1] and remainder[-1].value == 0:
                remainder.pop()
        res_quotient = from_coeff_vec(quotient)
        return res_quotient, remainder
            
def poly_div_poly(self: list[fr.Fr], divisor: list[fr.Fr]):
        res, remainder = divide_with_q_and_r(self,divisor)
        return res

def rand_poly(d):
    random.seed(42)
    random_coeffs = [fr.Fr.from_repr(random.random) for _ in range(d + 1)]
    return from_coeff_vec(random_coeffs)

def ln_without_floats(a):
    # log2(a) * ln(2)
    return int(math.log2(a) * 69 / 100)

# Evaluates `self` at the given `point` in `Self::Point`.
def evaluate(self, point: fr.Fr):
    zero = fr.Fr.zero()
    if len(self) == 0:
        return zero
    elif point.value == 0:
        return self[0]
    return horner_evaluate(self, point)

# Horner's method for polynomial evaluation
def horner_evaluate(poly_coeffs: list[fr.Fr], point: fr.Fr):
    result = fr.Fr.zero()
    for coeff in reversed(poly_coeffs):
        result = result.mul(point)
        result = result.add(coeff)
    return result

def poly_add_poly_mul_const(self: list[fr.Fr], f: fr.Fr, other: list[fr.Fr]):
    if len(self) == 0:
        self = other[:]
        for i in range(len(self)):
            self[i] = self[i].mul(f)
        return self
    elif len(other) == 0:
        return self
    elif len(self) >= len(other):
        pass
    else:
        zero = fr.Fr.zero()
        self = resize(self,len(other),zero)
    for i in range(len(other)):
        temp = f.mul(other[i])
        self[i] = self[i].add(temp)
    self = from_coeff_vec(self)
    return self

# Given a vector of field elements {v_i}, compute the vector {coeff * v_i^(-1)}
# This method is explicitly single core.
def serial_batch_inversion_and_mul(v: list[fr.Fr], coeff: fr.Fr):
    prod = []
    tmp = coeff.one()

    for f in v:
        if f.value != 0:
            tmp = tmp.mul(f)
            prod.append(tmp)

    # Invert `tmp`.
    tmp = fr.Fr.inverse(tmp)

    # Multiply product by coeff, so all inverses will be scaled by coeff
    tmp = tmp.mul(coeff)
    rev_prod = list(reversed(prod[:-1]))
    rev_prod.append(coeff.one())
    # Backwards, skip last element, fill in one for last term.
    for i,(f, s) in enumerate(zip(reversed(v), rev_prod)):
        if f.value != 0:
            # tmp := tmp * f; f := tmp * s = 1/f
            new_tmp = tmp.mul(f)
            f = tmp.mul(s)
            tmp = new_tmp
            v[len(v) - 1 - i] = f  # Update the value of v with the new result

# Given a vector of field elements {v_i}, compute the vector {coeff * v_i^(-1)}
def batch_inversion_and_mul(v: list[fr.Fr], coeff: fr.Fr):
    serial_batch_inversion_and_mul(v, coeff)

# Given a vector of field elements {v_i}, compute the vector {v_i^(-1)}
def batch_inversion(v: list[fr.Fr]):
    one = v[0].one()
    batch_inversion_and_mul(v, one)

# The first lagrange polynomial has the expression:
# L_0(X) = mul_from_1_to_(n-1) [(X - omega^i) / (1 - omega^i)]
#
# with `omega` being the generator of the domain (the `n`th root of unity).
#
# We use two equalities:
#   1. `mul_from_2_to_(n-1) [1 / (1 - omega^i)] = 1 / n`
#   2. `mul_from_2_to_(n-1) [(X - omega^i)] = (X^n - 1) / (X - 1)`
# to obtain the expression:
# L_0(X) = (X^n - 1) / n * (X - 1)
def compute_first_lagrange_evaluation(
    domain: Radix2EvaluationDomain, z_h_eval: fr.Fr, z_challenge: fr.Fr):
    one = z_challenge.one()
    n_fr = fr.Fr.from_repr(domain.size)  
    z_challenge_sub_one = z_challenge.sub(one)
    denom = n_fr.mul(z_challenge_sub_one)
    denom_in = fr.Fr.inverse(denom)
    res = z_h_eval.mul(denom_in)
    return res  

def MSM(bases:list[AffinePointG1], scalars:list[fr.Fr], params):
    size = min(len(bases), len(scalars))
    fq_elem = bases[0].x
    scalars = scalars[:size]
    bases = bases[:size]
    scalars_and_bases_iter = [(s, b) for s, b in zip(scalars, bases) if not s==0]

    c = 3 if size < 32 else ln_without_floats(size) + 2
    num_bits = fr.Fr.MODULUS_BITS
    fr_one = params.one()
    fr_one = fr_one.into_repr()

    zero:ProjectivePointG1 = ProjectivePointG1.zero(fq_elem)
    window_starts = list(range(0, num_bits, c))

    window_sums = []
    for w_start in window_starts:
        res = zero
        buckets = [zero for _ in range((1 << c) - 1)]

        for org_scalar, org_base in scalars_and_bases_iter:
            scalar = copy.copy(org_scalar)
            base = copy.copy(org_base)
            if scalar.value == fr_one:
                if w_start == 0:
                    res = res.add_assign_mixed(base)
            else:
                # We right-shift by w_start, thus getting rid of the lower bits
                scalar.value >>= w_start
                # We mod the remaining bits by 2^{window size}, thus taking `c` bits.
                scalar.value %= (1 << c)
                if scalar.value != 0:
                    buckets[scalar.value - 1] = buckets[scalar.value - 1].add_assign_mixed(base)

        running_sum:ProjectivePointG1 = ProjectivePointG1.zero(fq_elem)
        for b in reversed(buckets):
            running_sum = running_sum.add_assign(b)
            res = res.add_assign(running_sum)

        window_sums.append(res)

    lowest = window_sums[0]

    total:ProjectivePointG1 = zero
    for sum_i in reversed(window_sums[1:]):
        total = total.add_assign(sum_i)
        for _ in range(c):
            total = total.double()
    total = lowest.add_assign(total)
    return total
