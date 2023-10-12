import gmpy2
import math
from dataclasses import dataclass
from bls12_381.fr import FrParameters as Fq
from bls12_381.edwards_fr import FrParameters as Fr
from field import field

fr = Fr()
fq = Fq()

GENERATOR_X = field(gmpy2.mpz(8076246640662884909881801758704306714034609987455869804520522091855516602923),
                    fq)
GENERATOR_Y = field(gmpy2.mpz(13262374693698910701929044844600465831413122818447359594527400194675274060458),
                    fq)
@dataclass
class EdwardsParameters:
    COEFF_A = field(gmpy2.mpz(41515536288062376014772236515869989659801672835942349428353392091902613913603),
                    fq)

    COEFF_D = field(gmpy2.mpz(39791098284436708367363857153769964807626706786434173600033581261390821521072),
                    fq)
    
    COFACTOR = [8]

    COFACTOR_INV = field(gmpy2.mpz(819310549611346726241370945440405716213240158234039660170669895299022906775),
                         fr)
    
    # AFFINE_GENERATOR_COEFFS = (GENERATOR_X, GENERATOR_Y)
    AFFINE_GENERATOR_COEFFS = (GENERATOR_X, GENERATOR_Y)
