from dataclasses import dataclass
import gmpy2
import math

@dataclass
class FftParameters:
    TWO_ADICITY: int = 32

    TWO_ADIC_ROOT_OF_UNITY: gmpy2.mpz = gmpy2.mpz(
        "0x5bf3adda19e9b27b0af53ae352a31e645b1b4c801819d7ecb9b58d8c5f0e466a", 16

    )

@dataclass
class FpParameters:

    MODULUS: gmpy2.mpz = gmpy2.mpz(
        52435875175126190479447740508185965837690552500527637822603658699938581184513
    )

    MODULUS_BITS: int = 255

    CAPACITY: int = MODULUS_BITS - 1

    REPR_SHAVE_BITS: int = 1

    R:gmpy2.mpz = gmpy2.mpz(
        10920338887063814464675503992315976177888879664585288394250266608035967270910
    )

    R2: gmpy2.mpz = gmpy2.mpz(
        "0x748d9d99f59ff1105d314967254398f02b6cedcb87925c23c999e990f3f29c6d", 16
    )

    R_INV: gmpy2.mpz = gmpy2.mpz(
        12549076656233958353659347336803947287922716146853412054870763148006372261952
    )

    GENERATOR: gmpy2.mpz = gmpy2.mpz(
        24006497034320510773280787438025867407531605151569380937148207556313189711857
    )

    MODULUS_MINUS_ONE_DIV_TWO: gmpy2.mpz = gmpy2.mpz(
        "0x39f6d3a994cebea4199cec0404d0ec02a9ded2017fff2dff7fffffff80000000", 16
    )

    T: gmpy2.mpz = gmpy2.mpz(
        12208678567578594777604504606729831043093128246378069236549469339647
    )

    T_MINUS_ONE_DIV_TWO: gmpy2.mpz = gmpy2.mpz(
        6104339283789297388802252303364915521546564123189034618274734669823
    )

@dataclass
class FrParameters(FftParameters,FpParameters):
    #256bits
    BYTE_SIZE:int = 32
    # Return the Multiplicative identity
    def one(cls):
        return cls.R
    
    # Returns the 2^s root of unity.
    def two_adic_root_of_unity(self):
        return self.TWO_ADIC_ROOT_OF_UNITY 

    # Returns the 2^s * small_subgroup_base^small_subgroup_base_adicity root of unity
    # if a small subgroup is defined.
    def large_subgroup_root_of_unity():
        pass

    # Returns the multiplicative generator of `char()` - 1 order.
    def multiplicative_generator(cls):
        return cls.GENERATOR

    # Returns the root of unity of order n, if one exists.
    # If no small multiplicative subgroup is defined, this is the 2-adic root of unity of order n
    # (for n a power of 2).
    def get_root_of_unity(self,n):
        size = 2 ** (n.bit_length()-1)
        log_size_of_group = int(math.log2(size))

        if n != size or log_size_of_group > self.TWO_ADICITY:
            return None

        # Compute the generator for the multiplicative subgroup.
        # It should be 2^(log_size_of_group) root of unity.
        omega = self.two_adic_root_of_unity()
        R_inv=gmpy2.invert(self.R,self.MODULUS)
        for _ in range(log_size_of_group, self.TWO_ADICITY):
            #modsquare
            omega *=omega
            omega *=R_inv
            omega %=self.MODULUS
        return omega



