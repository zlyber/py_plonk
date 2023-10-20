import gmpy2
import math
from dataclasses import dataclass
from field import field


class Fq(field):
    def __init__(self, value:gmpy2.mpz):
        self.value = value

    TWO_ADICITY: int = 1

    TWO_ADIC_ROOT_OF_UNITY: gmpy2.mpz = gmpy2.mpz(
        "0x40ab3263eff02060ef148d1ea0f4c069eca8f3318332bb7a7e83a49a2e99d69032b7fff2ed47fffd43f5fffffffcaaae",
        16
    )

    MODULUS: gmpy2.mpz = gmpy2.mpz(
        4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
    )

    MODULUS_BITS: int = 381

    CAPACITY: int = MODULUS_BITS - 1

    REPR_SHAVE_BITS: int = 3

    R: gmpy2.mpz = gmpy2.mpz(
        3380320199399472671518931668520476396067793891014375699959770179129436917079669831430077592723774664465579537268733
    )

    R2: gmpy2.mpz = gmpy2.mpz(
        2708263910654730174793787626328176511836455197166317677006154293982164122222515399004018013397331347120527951271750
    )

    R_INV: gmpy2.mpz = gmpy2.mpz(
        3231460744492646417066832100176244795738767926513225105051837195607029917124509527734802654356338138714468589979680
    )
    GENERATOR: gmpy2.mpz = gmpy2.mpz(
        2758230843577277949620073511305048635578704962089743514587482222134842183668501798417467556318533664893264801977679
    )

    MODULUS_MINUS_ONE_DIV_TWO: gmpy2.mpz = gmpy2.mpz(
        "0xd0088f51cbff34d0258dd3db21a5d66bb23ba5c279c2895fb39869507b587b12f55ffff58a9ffff0dcff7fffffffd555",
        16
    )

    T: gmpy2.mpz = gmpy2.mpz(
        "0xd0088f51cbff34d0258dd3db21a5d66bb23ba5c279c2895fb39869507b587b12f55ffff58a9ffff0dcff7fffffffd555",
        16
    )

    T_MINUS_ONE_DIV_TWO: gmpy2.mpz = gmpy2.mpz(
        "0x680447a8e5ff9a6092c6e9ed90d2eb35d91dd2e13ce144afd9cc34a83dac3d897aaffffac54ffff0ee7fbfffffffeaaa",
        16
    )

    #384bits
    BYTE_SIZE: int = 48
    # # Return the Multiplicative identity
    # def one(cls):
    #     return cls(cls.R)
    
    # # Returns the 2^s root of unity.
    # def two_adic_root_of_unity(self):
    #     return self.TWO_ADIC_ROOT_OF_UNITY 

    # # Returns the 2^s * small_subgroup_base^small_subgroup_base_adicity root of unity
    # # if a small subgroup is defined.
    # def large_subgroup_root_of_unity():
    #     pass

    # # Returns the multiplicative generator of `char()` - 1 order.
    # def multiplicative_generator(cls):
    #     return cls.GENERATOR

    # # Returns the root of unity of order n, if one exists.
    # # If no small multiplicative subgroup is defined, this is the 2-adic root of unity of order n
    # # (for n a power of 2).
    # def get_root_of_unity(self,n):
    #     size = 2 ** (n.bit_length()-1)
    #     log_size_of_group = int(math.log2(size))

    #     if n != size or log_size_of_group > self.TWO_ADICITY:
    #         return None

    #     # Compute the generator for the multiplicative subgroup.
    #     # It should be 2^(log_size_of_group) root of unity.
    #     omega = self.two_adic_root_of_unity()
    #     R_inv=gmpy2.invert(self.R,self.MODULUS)
    #     for _ in range(log_size_of_group, self.TWO_ADICITY):
    #         #modsquare
    #         omega *=omega
    #         omega *=R_inv
    #         omega %=self.MODULUS
    #     return omega

FQ_ONE = gmpy2.mpz(1)
FQ_ZERO = gmpy2.mpz(0)
COEFF_A = gmpy2.mpz(0)
COEFF_B = gmpy2.mpz(4)