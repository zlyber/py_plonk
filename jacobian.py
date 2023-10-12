from dataclasses import dataclass
from typing import List
from field import field
from structure import AffinePointG1
from transcript import transcript
import gmpy2
from bls12_381 import fq


@dataclass
class ProjectivePointG1:
    x: field
    y: field
    z: field

    @classmethod
    #Returns the point at infinity, which always has Z = 0.
    def zero(cls,field_elemnt:field):
        x=field_elemnt.one()
        y=field_elemnt.one()
        z=field.zero(field_elemnt)
        return cls(x,y,z)
    
    def is_zero(self):
        return self.z.value == 0
    
    def double(self):
        if self.is_zero():
            return self

        # if fq.COEFF_A == 0:
        # A = X1^2
        a = self.x.square()

        # B = Y1^2
        b = self.y.square()

        # C = B^2
        c = b.square()

        # D = 2*((X1+B)^2-A-C)
        mid1 = self.x.add(b)
        mid1 = mid1.square()
        mid2 = mid1.sub(a)
        mid2 = mid2.sub(c)
        d = mid2.double()

        # E = 3*A
        mid1 = a.double()
        mid2 = mid1.double()
        e = a.add(mid2)

        # F = E^2
        f = e.square()

        # Z3 = 2*Y1*Z1
        mid1 = self.z.mul(self.y)
        self.z = mid1.double()

        # X3 = F-2*D
        mid1 = f.sub(d)
        self.x = mid1.sub(d)

        # Y3 = E*(D-X3)-8*C
        mid1 = d.sub(self.x)
        mid2 = c.double()
        mid2 = mid2.double()
        mid2 = mid2.double()
        mid3 = mid1.mul(e)
        self.y = mid3.sub(mid2)
        
        return self
        # else:
        #     # http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
        #     # XX = X1^2
        #     xx = self.x.square()

        #     # YY = Y1^2
        #     yy = self.y.square()

        #     # YYYY = YY^2
        #     yyyy = yy.square()

        #     # ZZ = Z1^2
        #     zz = self.z.square()

        #     # S = 2*((X1+YY)^2-XX-YYYY)
        #     mid1 = self.x.add(yy)
        #     mid1 = mid1.square()
        #     mid2 = mid1.sub(xx)
        #     mid1 = mid2.sub(yyyy)
        #     s = mid1.double()

        #     # M = 3*XX+a*ZZ^2
        #     mid1 = xx.double()
        #     mid1 = mid1.add(xx)
        #     mid2 = zz.square()

        #     m = xx + xx + xx + P.mul_by_a(zz.square())

        #     # T = M^2-2*S
        #     t = m.square() - s.double()

        #     # X3 = T
        #     self.x = t
        #     # Y3 = M*(S-T)-8*YYYY
        #     old_y = self.y
        #     self.y = m * (s - t) - yyyy.double_in_place().double_in_place().double_in_place()
        #     # Z3 = (Y1+Z1)^2-YY-ZZ
        #     self.z = (old_y + self.z).square() - yy - zz
        #     return self

    def to_affine(p:'ProjectivePointG1'):
        if p.is_zero():
            return AffinePointG1.zero(p.x.params)
        elif p.z.value == p.z.params.one():
            return AffinePointG1.new(p.x,p.y)
        else:
            zinv = field.inverse(p.z,p.z.params)
            zinv_squared = zinv.square()

            x = p.x.mul(zinv_squared)
            mid1 = zinv_squared.mul(zinv)
            y = p.y.mul(mid1)
            return AffinePointG1.new(x,y)
        
    def add_assign(self, other:'ProjectivePointG1'):
        if self.is_zero():
            self.x, self.y, self.z = other.x, other.y, other.z
            return

        if other.is_zero():
            return

        # http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
        # Works for all curves.

        # Z1Z1 = Z1^2
        z1z1 = self.z.square()

        # Z2Z2 = Z2^2
        z2z2 = other.z.square()

        # U1 = X1*Z2Z2
        u1 = self.x.mul(z2z2)

        # U2 = X2*Z1Z1
        u2 = other.x.mul(z1z1)

        # S1 = Y1*Z2*Z2Z2
        s1 = self.y.mul(other.z)
        s1 = s1.mul(z2z2)
        # S2 = Y2*Z1*Z1Z1
        s2 = other.y.mul(self.z)
        s2 = s2.mul(z1z1)

        if u1 == u2 and s1 == s2:
            # The two points are equal, so we double.
            self.double()
        else:
            # If we're adding -a and a together, self.z becomes zero as H becomes zero.

            # H = U2-U1
            h = u2.sub(u1)

            # I = (2*H)^2
            mid1 = h.double()
            i = mid1.square()

            # J = H*I
            j = h.mul(i)

            # r = 2*(S2-S1)
            r = s2.sub(s1)
            r = r.double()

            # V = U1*I
            v = u1.mul(i)

            # X3 = r^2 - J - 2*V
            mid1 = r.square()
            mid2 = mid1.sub(j)
            mid1 = v.double()
            self.x = mid2.sub(mid1)

            # Y3 = r*(V - X3) - 2*S1*J
            mid1 = v.sub(self.x)
            mid2 = r.mul(mid1)
            mid3 = s1.mul(j)
            mid4 = mid3.double()
            self.y = mid2.sub(mid4)


            # Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
            mid1 = self.z.add(other.z)
            mid2 = mid1.square()
            mid3 = mid2.sub(z1z1)
            mid4 = mid3.sub(z2z2)
            self.z = mid4.mul(h)

    def add_assign_mixed(self, other:'AffinePointG1'):
        if other.is_zero():
            return

        elif self.is_zero():
            self.x = other.x
            self.y = other.y
            self.z = self.x.one()

        else:
            # Z1Z1 = Z1^2
            z1z1 = self.z.square()

            # U2 = X2*Z1Z1
            u2 = other.x.mul(z1z1)

            # S2 = Y2*Z1*Z1Z1
            mid1 = other.y.mul(self.z)
            s2 = mid1.mul(z1z1)
            if self.x == u2 and self.y == s2:
                # The two points are equal, so we double.
                self.double()
            else:
                # If we're adding -a and a together, self.z becomes zero as H becomes zero.

                # H = U2-X1
                h = u2.sub(self.x)

                # HH = H^2
                hh = h.square()

                # I = 4*HH
                mid1 = hh.double()
                i = mid1.double()

                # J = H*I
                j = h.mul(i)

                # r = 2*(S2-Y1)
                mid1 = s2.sub(self.y)
                r = mid1.double()

                # V = X1*I
                v = self.x.mul(i)

                # X3 = r^2 - J - 2*V
                mid1 = r.square()
                mid2 = mid1.sub(j)
                mid3 = v.double()
                self.x = mid2.sub(mid3)

                # Y3 = r*(V-X3)-2*Y1*J
                mid1= j.mul(self.y)
                mid2 = mid1.double()
                mid3 = v.sub(self.x)
                mid4 = r.mul(mid3)
                self.y = mid4.sub(mid2)

                # Z3 = (Z1+H)^2-Z1Z1-HH
                mid1 = self.z.add(h)
                mid2 = mid1.square()
                mid3 = mid2.sub(z1z1)
                self.z = mid3.sub(hh)

    

