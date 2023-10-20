from dataclasses import dataclass
import gmpy2
import math
from field import field


    

@dataclass
class Fr(field):
    TWO_ADICITY: int = 32

    TWO_ADIC_ROOT_OF_UNITY: gmpy2.mpz = gmpy2.mpz(
        "0x4d6b87b1da259e207342261215ac260bb3524a6466112932aa9f02ab1d6124de", 16

    )

    MODULUS: gmpy2.mpz = gmpy2.mpz(
        6554484396890773809930967563523245729705921265872317281365359162392183254199
    )

    MODULUS_BITS: int = 252

    CAPACITY: int = MODULUS_BITS - 1

    REPR_SHAVE_BITS: int = 4

    R:gmpy2.mpz = gmpy2.mpz(
        "0x9a6fc6f479155c60932514eeeb8814f4f315d62f66b6e75025f80bb3b99607d9", 16
    )

    R2: gmpy2.mpz = gmpy2.mpz(
        "0x4f6547b8d127688069dab7fac026e9a551b0cef09ce3fc2667719aa495e57731", 16
    )
    # TODO:calculate R_INV
    R_INV: gmpy2.mpz = gmpy2.mpz(
        12549076656233958353659347336803947287922716146853412054870763148006372261952
    )

    GENERATOR: gmpy2.mpz = gmpy2.mpz(
        "0xe70cbdc7dccf3ac05fa8cc968193ccbbbf4aa36101f13a58720b1b19d49ea8f1"
    )

    MODULUS_MINUS_ONE_DIV_TWO: gmpy2.mpz = gmpy2.mpz(
        3277242198445386904965483781761622864852960632936158640682679581196091627099
    )

    T: gmpy2.mpz = MODULUS_MINUS_ONE_DIV_TWO

    T_MINUS_ONE_DIV_TWO: gmpy2.mpz = gmpy2.mpz(
        1638621099222693452482741890880811432426480316468079320341339790598045813549
    )

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



