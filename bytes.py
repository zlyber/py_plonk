from io import IOBase
import gmpy2
import struct

def write(content:gmpy2.mpz,writer:list):
    bytes_content = gmpy2.to_binary(content)
    little_endian_bytes = bytes_content[::-1]
    byte_list = [byte for byte in little_endian_bytes]
    writer.extend(byte_list)
    return writer

def read(reader):
    integer_data = gmpy2.mpz(struct.unpack_from('<Q', reader)[0])
    return integer_data