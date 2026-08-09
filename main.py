"""
Reed-Solomon (RS) Code Implementation & Finite Field Arithmetic

Author: Nathan Kimutai(tekindet)

This module provides a standalone, zero-dependency implementation of Reed-Solomon
error-correcting codes over Galois Field GF(2^8) using polynomial arithmetic
"""
class RSCodec:

    def __init__(self,n,k):
        self.nsym = 4
        self.n = n

    def encode(self,msg):

        msg_bytes = [ord(c) for c in msg]

        M_x = GF256Poly(msg_bytes)
        G_x = build_generator_poly(self.nsym)

        # x^4 for parity symbols
        shift_factor = GF256Poly([1] + [0] * self.nsym)

        M_shifted = M_x * shift_factor

        Q_x,P_x = M_shifted.divmod(G_x)

        # this is the codeword, this is what we are sending
        # to the receiver
        C_x = M_shifted + P_x

        # todo(test): introduce noise in a channel, also 
        # refactor it so that it is easy to use

        #print(f"Final Codeword C(x):\n{C_x}")

        #Q_verify,R_verify = C_x.divmod(G_x)

        return C_x

    def decode(self,msg):

        R_x = GF256Poly(msg)

        rems = []

        for i in range(self.nsym):

            alpha_i = GF256.pow(2,i)

            factor = GF256Poly([1,alpha_i])

            Q_x,P_x = R_x.divmod(factor)

            rems.append(P_x)


        return rems

    def berlekamp_massey(synds):

        # Step 1 : Find the error locator polynomial
        E_x = GF256Poly([])

        v = self.nsym / 2
        for j in range(v):
            


        # Step 2 : Find the error locations(Chien search)

        # Step 3 : Calculate the error magnitudes(Forney Algo)

        return E_x


def build_generator_poly(nsym):
    G_x = GF256Poly([1])

    for i in range(nsym):

        alpha_i = GF256.pow(2,i)

        factor = GF256Poly([1,alpha_i])

        G_x *= factor

    return G_x

class GF256:
    PRIM = 0x11D

    log = [0] * 256
    exp = [0] * 512

    x = 1

    for i in range(255):
        log[x] = i
        exp[i] = x

        x <<= 1

        if x & 0x100:
            x ^= PRIM

    for i in range(255,512):
        exp[i] = exp[i - 255]

    @classmethod
    def mul(cls,a,b):
        if a == 0 or b == 0: return 0
        return cls.exp[cls.log[a] + cls.log[b]]

    @classmethod
    def add(cls,a,b):
        return a ^ b

    @classmethod
    def sub(cls,a,b):
        return a ^ b

    @classmethod
    def div(cls,a,b):
        if b == 0: raise ZeroDivisionError("division by zero in galois field")
        if a == 0: return 0
        return cls.exp[(cls.log[a] + 255 - cls.log[b]) % 255]

    @classmethod
    def pow(cls,a,power):
        return cls.exp[(cls.log[a] * power) % 255]


class GF256Poly:
    def __init__(self,coeffs):

        i = 0
        while i < len(coeffs) - 1 and coeffs[i] == 0:
            i = i + 1
        self.coeffs = coeffs[i:]

    def __add__(self,other):
        len_diff = abs(len(self.coeffs) - len(other.coeffs))

        if len(self.coeffs) < len(other.coeffs):
            c1 = [0] * len_diff + self.coeffs
            c2 = other.coeffs
        else:
            c1 = self.coeffs
            c2 = [0] * len_diff + other.coeffs

        res = [GF256.add(a,b) for a,b in zip(c1,c2)]

        return GF256Poly(res)

    def __mul__(self,other):
        # if you multiply two polynomials a and b of 
        # degree N-1 then the result will be of degree
        # 2N - 2 = 2(N-1)
        res = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i,a in enumerate(self.coeffs):
            for j,b in enumerate(other.coeffs):
                res[i + j] = GF256.add(res[i + j],GF256.mul(a,b)) 

        return GF256Poly(res)

    def divmod(self,divisor):
        # return a quotient and a reminder
        out = list(self.coeffs)
        # leading term of divisor
        norm = divisor.coeffs[0]

        for i in range(len(self.coeffs) - len(divisor.coeffs) + 1):
            out[i] = GF256.div(out[i],norm)
            coeff = out[i]

            if coeff != 0:
                for j in range(1,len(divisor.coeffs)):

                    out[i + j] = GF256.add(
                            out[i + j],
                            GF256.mul(divisor.coeffs[j],coeff))

        separator = len(self.coeffs) - len(divisor.coeffs) + 1

        return GF256Poly(out[:separator]),GF256Poly(out[separator:])

        

    def __repr__(self):

        if self.coeffs == [0]:
            return "0"

        terms = []

        for idx,char in enumerate(self.coeffs):
            if char == 0:
                continue
            deg = len(self.coeffs) - 1

            power = deg - idx

            coeff_str = f"{char}" if(char != 1 or power == 0) else ""

            if power == 0:
                terms.append(f"{char}")
            elif power == 1:
                terms.append(f"{coeff_str}x")
            else:
                terms.append(f"{coeff_str}x^{power}")

        return " + ".join(terms) if terms else "0"

if __name__ == "__main__":
    msg = "this is a test message"

    kodek = RSCodec(8,4)

    C_x = kodek.encode(msg)

    # todo : need to corrupt the message to see if it works

    encoded_msg = bytes(C_x.coeffs)

    corrupted_msg = bytearray(encoded_msg)
    #corrupted_msg[0] ^= 0xFF

    rems = kodek.decode(corrupted_msg)
    print(rems)





