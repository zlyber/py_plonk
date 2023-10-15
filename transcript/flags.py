from abc import ABC, abstractmethod

Infinity : int = 0
PositiveY: int = 1
NegativeY: int = 2

class Flags(ABC):
    @property
    @classmethod
    @abstractmethod
    def BIT_SIZE(cls):
        pass

    @abstractmethod
    def u8_bitmask(self):
        pass

    @classmethod
    @abstractmethod
    def from_u8(cls, value):
        pass

    @classmethod
    def from_u8_remove_flags(self,value):
        flags = self.from_u8(value)
        if flags:
            bitmask = flags.u8_bitmask()
            value &= ~bitmask
        return flags

class SWFlags(Flags):
    BIT_SIZE = 2
    def __init__(self,flag):
        self.flag = flag

    @classmethod
    def infinity(cls):
        return SWFlags(Infinity)
    @classmethod
    def from_y_sign(cls,is_positive:bool):
        if is_positive:
            return cls(flag = PositiveY)
        else:
            return cls(flag = NegativeY)
        
    def u8_bitmask(self):
        mask = 0
        #infinity
        if self.flag == Infinity:
            mask |= 1 << 6
        #positive
        elif self.flag == PositiveY:
            mask |= 1 << 7
        return mask

    @classmethod
    def from_u8(cls, value):
        x_sign = (value >> 7) & 1 == 1
        is_infinity = (value >> 6) & 1 == 1
        if x_sign and is_infinity:
            return None
        elif not x_sign and is_infinity:
            return cls(Infinity)
        elif x_sign and not is_infinity:
            return cls(PositiveY)
        else:
            return cls(NegativeY)  
    

class EmptyFlags(Flags):
    BIT_SIZE = 0
    
    @classmethod
    def u8_bitmask(cls):
        return 0

    @classmethod
    def from_u8(cls, value):
        if (value >> 7) == 0:
            return cls()
        else:
            return None


