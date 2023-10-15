from dataclasses import dataclass
from transcript import strobe
from structure import AffinePointG1
from field import from_random_bytes
import struct

MERLIN_PROTOCOL_LABEL = b"Merlin v1.0"

def write_u32(buf, n):
    # little endian
    packed_data = struct.pack("<I", n)
    buf[:4] = packed_data
    return buf

def encode_usize_as_u32(x):
    if x > 0xFFFFFFFF:
        raise ValueError("Value too large to fit in u32")
    buf=[0,0,0,0]
    return write_u32(buf,x)

@dataclass
class Transcript:
    strobe:strobe.Strobe128

    @classmethod
    def new(cls,label):
        transcript = cls(strobe.Strobe128.new(MERLIN_PROTOCOL_LABEL))
        transcript.append_message(b"dom-sep", label)
        return transcript

    def append_message(self, label, message):
        data_len = encode_usize_as_u32(len(message))
        self.strobe.meta_ad(label, False)
        self.strobe.meta_ad(data_len, True)
        self.strobe.ad(message, False)
    
    def append_pi(self, label):
        pi_bytes = [1, 0, 0, 0, 0, 0, 0, 0, 71, 2, 0, 0, 0, 0, 0, 0, 175, 183, 188, 23, 25,
                    155, 195, 172, 246, 232, 114, 89, 255, 117, 123, 226, 161, 86, 112, 46,
                    124, 69, 74, 38, 91, 197, 152, 114, 132, 231, 67, 72]
        self.append_message(label,pi_bytes)
    
    def append(self,label,item):
        bytes = []
        bytes = item.serialize(bytes)
        self.append_message(label,bytes)
    
    def challenge_bytes(self, label, dest):
        data_len = encode_usize_as_u32(len(dest))
        self.strobe.meta_ad(label, False)
        self.strobe.meta_ad(data_len, True)
        modified_dest = self.strobe.prf(dest, False)
        return modified_dest
    
    def challenge_scalar(self, label: bytes, params):
        size = params.MODULUS_BITS // 8
        buf = bytes([0] * size)
        modified_buf = self.challenge_bytes(label, buf)
        c_s = from_random_bytes(params,modified_buf)
        return c_s
