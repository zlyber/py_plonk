from dataclasses import dataclass
import struct
# Strobe R value; security level 128 is hardcoded
STROBE_R = 166

FLAG_I = 1
FLAG_A = 1 << 1
FLAG_C = 1 << 2
FLAG_T = 1 << 3
FLAG_M = 1 << 4
FLAG_K = 1 << 5


KECCAK_F_ROUND_COUNT = 24

PLEN = 25

RHO = [1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44]

PI = [10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1]

RC = [
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808a,
    0x8000000080008000,
    0x000000000000808b,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008a,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000a,
    0x000000008000808b,
    0x800000000000008b,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800a,
    0x800000008000000a,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008
]


def transmute_state(st):
    # Assuming st is a bytearray with the same layout as [u64; 25] in Rust
    u64_values = []
    
    for i in range(0, len(st), 8):
        u64_values.append(struct.unpack("<Q", st[i:i+8])[0])
    
    return u64_values

def transmute_inverse(st):
    u8_array = bytearray()
    for value in st:
        #packing
        u8_values = struct.pack("<Q", value)  
        u8_array.extend(u8_values)  
    return u8_array

#we don't impl SIMD
def truncate_rc(rc):
    return rc


def rotate_left(value, n):
    return ((value << n) | (value >> (64 - n))) & 0xFFFFFFFFFFFFFFFF


def keccak_p(state, round_count):
    if round_count > KECCAK_F_ROUND_COUNT:
        raise ValueError("A round_count greater than KECCAK_F_ROUND_COUNT is not supported!")

    round_consts = RC[(KECCAK_F_ROUND_COUNT - round_count): KECCAK_F_ROUND_COUNT]

    for rc in round_consts:
        array = [0] * 5

        # Theta
        for x in range(5):
            for y in range(5):
                array[x] ^= state[5 * y + x]

        for x in range(5):
            for y in range(5):
                t1 = array[(x + 4) % 5]
                t2 = rotate_left(array[(x + 1) % 5],1)
                state[5 * y + x] ^= t1 ^ t2

        # Rho and Pi
        last = state[1]
        for x in range(24):
            array[0] = state[PI[x]]
            state[PI[x]] = rotate_left(last,RHO[x])
            last = array[0]

        # Chi
        for y_step in range(5):
            y = 5 * y_step
            for x in range(5):
                array[x] = state[y + x]

            for x in range(5):
                t1 = ~array[(x + 1) % 5]
                t2 = array[(x + 2) % 5]
                state[y + x] = array[x] ^ (t1 & t2)

        # Iota
        state[0] ^= truncate_rc(rc)

    return state


@dataclass
class Strobe128:
    state: bytearray
    pos: int
    pos_begin: int
    cur_flags: int

    @classmethod
    def new(cls,protocol_label):
        st = bytearray([0] * 200)
        st[0:6] = [1, STROBE_R + 2, 1, 0, 1, 96]
        st[6:18] = b"STROBEv1.0.2"
        # Simulate the keccak::f1600 function
        st=transmute_state(st)
        st=keccak_p(st,KECCAK_F_ROUND_COUNT)
        st=transmute_inverse(st)
        strobe = cls(st, 0, 0, 0)
        strobe.meta_ad(protocol_label, False)

        return strobe

    def run_f(self):
        self.state[self.pos] ^= self.pos_begin
        self.state[self.pos + 1] ^= 0x04
        self.state[STROBE_R + 1] ^= 0x80
        self.state = transmute_state(self.state)
        self.state = keccak_p(self.state,KECCAK_F_ROUND_COUNT)
        self.state = transmute_inverse(self.state)
        self.pos = 0
        self.pos_begin = 0

    def absorb(self, data):
        for byte in data:
            self.state[self.pos] ^= byte
            self.pos += 1
            if self.pos == STROBE_R:
                self.run_f()
    
    def squeeze(self, data):
        data_array = bytearray(data)
        for i in range(len(data)):
            data_array[i] = self.state[self.pos]
            self.state[self.pos] = 0
            self.pos += 1
            if self.pos == STROBE_R:
                self.run_f()
        modified_data = bytes(data_array)
        return modified_data

    def begin_op(self, flags, more):
        # judge whether go on
        if more:
            assert (
                self.cur_flags == flags
            ), f"You tried to continue op {bin(self.cur_flags)} but changed flags to {bin(flags)}"
            return

        # Skip adjusting direction information (we just use AD, PRF)
        assert (
            (flags & FLAG_T) == 0
        ), "You used the T flag, which this implementation doesn't support"

        old_begin = self.pos_begin
        self.pos_begin = self.pos + 1
        self.cur_flags = flags

        self.absorb([old_begin, flags])
        # Force running F if C or K is set
        force_f = (flags & (FLAG_C | FLAG_K)) != 0

        if force_f and self.pos != 0:
            self.run_f()

    def meta_ad(self, data, more):
        self.begin_op(FLAG_M | FLAG_A, more)
        self.absorb(data)

    def ad(self,data,more):
        self.begin_op(FLAG_A, more)
        self.absorb(data)
    
    def prf(self, data, more):
        flags = FLAG_I | FLAG_A | FLAG_C
        self.begin_op(flags, more)
        modified_data = self.squeeze(data)
        return modified_data
