from dataclasses import dataclass
from structure import UniversalParams,OpenProof
from jacobian import ProjectivePointG1
from field import field
from bls12_381 import fr,fq
from typing import List
from arithmetic import MSM,skip_leading_zeros_and_convert_to_bigints,convert_to_bigints,rand_poly,poly_add_poly_mul_const,evaluate,from_coeff_vec,poly_div_poly
from plonk_core.src.proof_system.linearisation_poly import ProofEvaluations
import random

class Randomness:
    def __init__(self, blind_poly: List[fr.Fr]):
        self.blind_poly = blind_poly

    @classmethod
    def empty(cls):
        return cls([])

    @classmethod
    def calculate_hiding_polynomial_degree(cls, hiding_bound):
        return hiding_bound + 1

    def push(self, a):
        self.blind_poly.append(a)

    @classmethod
    def rand(cls, hiding_bound):
        hiding_poly_degree = cls.calculate_hiding_polynomial_degree(hiding_bound)
        return cls(blind_poly = rand_poly(hiding_poly_degree))
    
    def add_assign(self, f:field, other: 'Randomness'):
        self.blind_poly = poly_add_poly_mul_const(self.blind_poly, f, other.blind_poly)

class Commitment:
    def __init__(self,value):
        self.value = value
    @classmethod
    def commit(cls,powers,polynomial:list[fr.Fr],hiding_bound,params):
        num_leading_zeros, plain_coeffs = skip_leading_zeros_and_convert_to_bigints(polynomial)
        commitment:ProjectivePointG1 = MSM(
            powers[0][num_leading_zeros:],
            plain_coeffs,
            params
        )
        randomness = Randomness.empty()
        if hiding_bound:
            randomness = Randomness.rand(hiding_bound)

        random_ints = convert_to_bigints(randomness.blind_poly)
        random_commitment:ProjectivePointG1 = MSM(powers[1],random_ints,params)
        random_commitment_affine = random_commitment.to_affine()
        commitment = commitment.add_assign_mixed(random_commitment_affine)
        commitment_affine = commitment.to_affine()
        return Commitment(value=commitment_affine),randomness
    
# On input a list of labeled polynomials and a query point, `open` outputs a proof of evaluation
# of the polynomials at the query point.
def open(
    ck: UniversalParams,
    labeled_polynomials: 'LabeledPoly',
    _commitments: 'LabeledCommitment',
    point,
    opening_challenge: field,
    rands,
    _rng=None
):
    combined_polynomial = []
    combined_rand = Randomness.empty()

    opening_challenge_counter = 0

    curr_challenge = opening_challenges(opening_challenge, opening_challenge_counter)
    opening_challenge_counter += 1

    for polynomial, rand in zip(labeled_polynomials, rands):

        combined_polynomial = poly_add_poly_mul_const(combined_polynomial,curr_challenge, polynomial.poly)
        combined_rand.add_assign(curr_challenge, rand)
        curr_challenge = opening_challenges(opening_challenge, opening_challenge_counter)
        opening_challenge_counter += 1

    powers = [ck.powers_of_g,ck.powers_of_gamma_g]
    proof = open_proof(powers, combined_polynomial, point, combined_rand)
    return proof

dataclass
class LabeledCommitment:
    def __init__(self,label,commitment):
        self.label = label
        self.commitment =commitment

    @classmethod
    def new(cls,label,commitment):
        return cls(label = label,commitment = commitment)

class LabeledPoly:
    def __init__(self, label, hiding_bound, poly):
        self.label = label
        self.hiding_bound = hiding_bound
        self.poly = poly

    @classmethod
    def new(cls, label, hiding_bound, poly):
        return cls(label=label, hiding_bound=hiding_bound, poly=poly)


def commit_poly(ck:UniversalParams,polys,params):
    random.seed(42)
    randomness = []
    labeled_comm = []
    for labeled_poly in polys:
        polynomial = labeled_poly.poly
        hiding_bound = labeled_poly.hiding_bound
        label = labeled_poly.label

        powers = [ck.powers_of_g,ck.powers_of_gamma_g]

        comm,rand = Commitment.commit(powers,polynomial,hiding_bound,params)
        labeled_comm.append(LabeledCommitment.new(label,comm))
        randomness.append(rand)
    return labeled_comm,randomness

def opening_challenges(opening_challenge: fr.Fr, pow):
    return opening_challenge.pow(pow)

# Compute witness polynomial.
#
# The witness polynomial w(x) the quotient of the division (p(x) - p(z)) / (x - z)
# Observe that this quotient does not change with z because
# p(z) is the remainder term. We can therefore omit p(z) when computing the quotient.
def compute_witness_polynomial(p: List[fr.Fr], point: fr.Fr, randomness: Randomness):
    neg_p = point.neg()
    one = point.one()
    divisor = from_coeff_vec([neg_p,one])
    witness_polynomial = p[:]
    if len(p) != 0:
        witness_polynomial = poly_div_poly(p, divisor)
    random_witness_polynomial = None
    if len(randomness.blind_poly) != 0:
        random_p = randomness.blind_poly
        random_witness_polynomial = poly_div_poly(random_p, divisor)
    return witness_polynomial, random_witness_polynomial

def open_with_witness_polynomial(
    powers, 
    point: fr.Fr, 
    randomness: Randomness,
    witness_polynomial, 
    hiding_witness_polynomial):

    num_leading_zeros, witness_coeffs =skip_leading_zeros_and_convert_to_bigints(witness_polynomial)
    w = MSM(powers[0][num_leading_zeros:],witness_coeffs,point)
    random_v = None
    if hiding_witness_polynomial is not None:
        blinding_p = randomness.blind_poly
        blinding_evaluation = evaluate(blinding_p, point)
        random_witness_coeffs = convert_to_bigints(hiding_witness_polynomial)
        random_commit = MSM(powers[1],random_witness_coeffs,point)
        w = w.add_assign(random_commit)
        random_v = blinding_evaluation
    
    return OpenProof(w.to_affine(), random_v)

# On input a polynomial `p` and a point `point`, outputs a proof for the same.
def open_proof(powers, p: List[field], point: field, rand: Randomness):
    witness_poly, hiding_witness_poly = compute_witness_polynomial(p, point, rand)
    proof = open_with_witness_polynomial(
            powers,
            point,
            rand,
            witness_poly,
            hiding_witness_poly,
        )
    return proof
@dataclass
class Proof:
    # Commitment to the witness polynomial for the left wires.
    a_comm: Commitment

    # Commitment to the witness polynomial for the right wires.
    b_comm: Commitment

    # Commitment to the witness polynomial for the output wires.
    c_comm: Commitment

    # Commitment to the witness polynomial for the fourth wires.
    d_comm: Commitment

    # Commitment to the permutation polynomial.
    z_comm: Commitment

    # Commitment to the lookup query polynomial.
    f_comm: Commitment

    # Commitment to first half of sorted polynomial
    h_1_comm: Commitment

    # Commitment to second half of sorted polynomial
    h_2_comm: Commitment

    # Commitment to the lookup permutation polynomial.
    z_2_comm: Commitment

    # Commitment to the quotient polynomial.
    t_1_comm: Commitment

    # Commitment to the quotient polynomial.
    t_2_comm: Commitment

    # Commitment to the quotient polynomial.
    t_3_comm: Commitment

    # Commitment to the quotient polynomial.
    t_4_comm: Commitment

    # Commitment to the quotient polynomial.
    t_5_comm: Commitment

    # Commitment to the quotient polynomial.
    t_6_comm: Commitment

    # Commitment to the quotient polynomial.
    t_7_comm: Commitment

    # Commitment to the quotient polynomial.
    t_8_comm: Commitment

    # Batch opening proof of the aggregated witnesses
    aw_opening: OpenProof

    # Batch opening proof of the shifted aggregated witnesses
    saw_opening: OpenProof

    # Subset of all of the evaluations added to the proof.
    evaluations: ProofEvaluations
