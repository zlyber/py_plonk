import gmpy2
import math
from typing import List, TypeVar,Generic,get_args 
from dataclasses import dataclass
from transcript import flags
from bls12_381 import fq,fr
from serialize import buffer_byte_size
from bytes import write,read
@dataclass
class field:

    value:gmpy2.mpz
    params:any

    @classmethod
    def zero(cls,params):
        return cls(gmpy2.mpz(0),params)

    def one(self):
        return field(self.params.one(),self.params)
    
    def add(self,b):
        res = self.value + b.value
        res %= self.params.MODULUS
        return field(res,self.params)
    
    def sub(self,b):
        res = gmpy2.f_mod(self.value - b.value, self.params.MODULUS)
        return field(res,self.params)
    
    def neg(self):
        if self.value!=0:
            temp = field(self.params.MODULUS,self.params)
            res = temp.sub(self)
            return res
        else:
            return self
        
    def double(self):
        res = self.add(self)
        return res

    #mongomery mul
    def mul(self, b):
        res = self.value * b.value
        #keep the result in mongomery form, same as mongomery reduce
        res = res * self.params.R_INV
        res %= self.params.MODULUS
        res = field(res,self.params)
        return res
    
    def square(self):
        self = self.mul(self)
        return self
    
    def pow(self, exp):
        res = field(self.params.one(),self.params)
        for i in range(63,-1,-1):
            #modsquare
            res = res.square()
            if ((exp >> i) & 1) == 1:
                res = res.mul(self)
        return res
    
    
    #new
    @classmethod
    def from_repr(cls, r, params):
        r=field(r,params)
        if r == 0:
            return r
        else:  
            R2 = field(params.R2,params)
            r = r.mul(R2)
            return r
        
    #Montgomery Reduction
    def into_repr(self):
        if self.value == 0:
            return self.value
        else:
            R_INV = field(self.params.R_INV,self.params)
            self.mul(R_INV)
            return self.value

    @classmethod
    def inverse(cls,self,params):
        if self==0:
            print("cannot invert 0!\n")
            return  None
        #judge if self is already in field
        elif type(self)==field:
            x = field(self.value,self.params) 
            x.value = gmpy2.invert(x.value, params.MODULUS)
            return x
        else:
            x=field.from_repr(self,params)
            return field(gmpy2.invert(x.value, params.MODULUS),params)
    
    def write(self,writer):
        content = self.into_repr()
        new_writer = write(content,writer)
        return new_writer
    
    def serialize_with_flags(self, writer:list, flag):
        if flag.BIT_SIZE > 8:
            print("Not enough space")
            return

        output_byte_size = buffer_byte_size(self.params.MODULUS_BITS + flag.BIT_SIZE)

        bytes = bytearray(self.params.BYTE_SIZE)
        modified_bytes = self.write(bytes[:self.params.BYTE_SIZE])

        modified_bytes[output_byte_size - 1] |= flag.u8_bitmask()

        writer.extend(modified_bytes[:output_byte_size])
        return writer
    
    def serialize(self,writer):
        writer = self.serialize_with_flags(writer,flags.EmptyFlags)
        return writer
    
def deserialize(params,reader):
    output_byte_size = buffer_byte_size(params.MODULUS_BITS + flags.EmptyFlags.BIT_SIZE)

    masked_bytes = bytearray([0] * (params.BYTE_SIZE + 1))
    masked_bytes[:output_byte_size] = reader[:output_byte_size]

    flag = flags.EmptyFlags.from_u8_remove_flags(masked_bytes[output_byte_size - 1])

    element = read(masked_bytes)
    field_element = field.from_repr(element,params)
    return field_element, flag

    
def from_random_bytes(params, bytes: bytes):
    limbs = (len(bytes) + 1) // 8
    if flags.EmptyFlags.BIT_SIZE > 8:
        return None

    result_bytes = bytearray([0] * (limbs * 8 + 1))
    result_bytes[:len(bytes)] = bytes

    last_limb_mask = ((2 ** 64 - 1)>>params.REPR_SHAVE_BITS).to_bytes(8, byteorder='little')

    output_byte_size = buffer_byte_size(params.MODULUS_BITS + flags.EmptyFlags.BIT_SIZE)
    flag_location = output_byte_size - 1
    flag_location_in_last_limb = flag_location - (8 * (limbs - 1))

    last_bytes = result_bytes[8 * (limbs - 1):]

    flags_mask = 0xFF >> (8 - flags.EmptyFlags.BIT_SIZE)
    flag = 0

    for i, (b, m) in enumerate(zip(last_bytes, last_limb_mask)):
        if i == flag_location_in_last_limb:
            flag = b & flags_mask
        b &= m

    field_element,flag = deserialize(params,result_bytes[:limbs * 8])
    #flags_obj = flags.SWFlags.from_u8(flag)

    return field_element

        