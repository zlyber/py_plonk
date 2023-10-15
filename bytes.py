from io import IOBase
import gmpy2
import struct

def write(content:gmpy2.mpz,writer:list):
    limbs = len(writer)//8
    content_part = [[] for _ in range(limbs)]
    bytes_content = [[] for _ in range(limbs)]
    
    for i in range(limbs):
        content_part[i] = int(content & 0xFFFFFFFFFFFFFFFF)
        if content_part[i] == 0:
            bytes_content[i] = b'\x00' * 8  # 以8个字节的空字节填充
        else:
            bytes_content[i] = content_part[i].to_bytes(8, byteorder='little')
        content = content >> 64

    little_endian_bytes = [item for sublist in bytes_content for item in sublist]
    byte_list = [byte for byte in little_endian_bytes]

    for i in range(len(writer)):
        writer[i] = byte_list[i]
    return writer

def read(reader,params):
    limbs = (params.REPR_SHAVE_BITS+params.MODULUS_BITS)//64
    format_string = "<" + "Q" * limbs
    integer_data = struct.unpack_from(format_string, reader)
    num = 0
    for data in reversed(integer_data):
        num = num << 64 | data
    num = gmpy2.mpz(num)
    return num