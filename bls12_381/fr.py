import gmpy2
import math
from field import field
from transcript import flags
from serialize import buffer_byte_size
from bytes import read
class Fr(field):
    def __init__(self, value:gmpy2.mpz):
        self.value = value

    TWO_ADICITY: int = 32

    TWO_ADIC_ROOT_OF_UNITY: gmpy2.mpz = gmpy2.mpz(
        "0x5bf3adda19e9b27b0af53ae352a31e645b1b4c801819d7ecb9b58d8c5f0e466a", 16

    )

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
        3294906474794265442129797520630710739278575682199800681788903916070560242797
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

    #256bits
    BYTE_SIZE:int = 32
    
def deserialize(params,reader):
    output_byte_size = buffer_byte_size(params.MODULUS_BITS + flags.EmptyFlags.BIT_SIZE)

    masked_bytes = bytearray([0] * (params.BYTE_SIZE + 1))
    masked_bytes[:output_byte_size] = reader[:output_byte_size]

    flag = flags.EmptyFlags.from_u8_remove_flags(masked_bytes[output_byte_size - 1])

    element = read(masked_bytes,params)
    field_element = Fr.from_repr(element)
    return field_element, flag

    
def from_random_bytes(params, bytes: bytes):
    limbs = (len(bytes) + 1) // 8
    if flags.EmptyFlags.BIT_SIZE > 8:
        return None

    result_bytes = bytearray([0] * (limbs * 8 + 1))
    result_bytes[:len(bytes)] = bytes
    last_bytes_mask = bytearray(9)
    last_limb_mask = ((2 ** 64 - 1)>>params.REPR_SHAVE_BITS).to_bytes(8, byteorder='little')
    last_bytes_mask[:8] = last_limb_mask[:]
    output_byte_size = buffer_byte_size(params.MODULUS_BITS + flags.EmptyFlags.BIT_SIZE)
    flag_location = output_byte_size - 1
    flag_location_in_last_limb = flag_location - (8 * (limbs - 1))

    last_bytes = result_bytes[8 * (limbs - 1):]

    flags_mask = 0xFF >> (8 - flags.EmptyFlags.BIT_SIZE)
    flag = 0
    for i, (b, m) in enumerate(zip(last_bytes, last_bytes_mask)):
        if i == flag_location_in_last_limb:
            flag = b & flags_mask
        b &= m

    field_element,flag = deserialize(params,result_bytes[:limbs * 8])
    #flags_obj = flags.SWFlags.from_u8(flag)

    return field_element


