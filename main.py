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


    def berlekamp_massey_second(self, syndromes):
        # C(x) starts as 1
        C = GF256Poly([1])
        B = GF256Poly([1])
        L = 0
        m = 1

        for n in range(len(syndromes)):
            # Compute discrepancy d
            d = syndromes[n]
            for i in range(1, L + 1):
                if i <= len(C.coeffs) - 1:
                    coef = C.coeffs[-(i + 1)]
                    d = GF256.add(d, GF256.mul(coef, syndromes[n - i]))

            if d == 0:
                m += 1
            else:
                T = C
                # C(x) = C(x) + d * B(x) * x^m
                # Construct (d * B(x) * x^m)
                scaled_B_coeffs = [GF256.mul(d, c) for c in B.coeffs] + [0] * m
                scaled_B = GF256Poly(scaled_B_coeffs)
                C = C + scaled_B  # Polynomial addition in GF(2^8) is XOR

                if 2 * L <= n:
                    L = n + 1 - L
                    # B = T / d
                    d_inv = GF256.inv(d)
                    B = GF256Poly([GF256.mul(c, d_inv) for c in T.coeffs])
                    m = 1
                else:
                    m += 1

        return C

    def berlekamp_massey(self,synds):

        C = [1]
        B = [1]

        v = len(synds)

        for N in range(v):
            d = synds[N]

            for i in range(1,len(C)):
                d = GF256.add(d, GF256.mul(C[i], syndromes[N - i]))

                #t_poly = GF256Poly([C[i]]) * syndromes[N - i]

                #d = d + t_poly


            B = [0] + B

            if d != 0:
                T = list(C)

                len_diff = abs(len(B) - len(C))

                c1 = C + [0] * (len(B) - len(C)) if len(C) < len(B) else C
                c2 = B + [0] * (len(C) - len(B)) if len(B) < len(C) else B

                C = [GF256.add(a, GF256.mul(1, b)) for a, b in zip(c1, c2)]

                if 2 * (len(T) - 1) <= N:

                    #u,r = GF256Poly(T).divmod(d)
                    #B = u.coeffs

                    #quotient,remainder = d.divmod(GF256Poly(T))
                    B = [GF256.div(x, d) for x in T]
                    #B = quotient.coeffs
                    #B = [GF256Poly(T).divmod(d)[0].coeffs[0] for x in T]

                    #print(B)

        return GF256Poly(C)


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

        syndromes = []

        for i in range(self.nsym):

            alpha_i = GF256.pow(2,i)

            factor = GF256Poly([1,alpha_i])

            Q_x,P_x = R_x.divmod(factor)

            rem_val = P_x.coeffs[0] if P_x.coeffs else 0
            syndromes.append(rem_val)


        return syndromes

    def correct_errors(self,coeffs):
        """
        Can only correct upto t errors where t = (n - k) / 2
        """
        synds = self.decode(coeffs)

        largest = 0

        R_x = GF256Poly(synds)

        #largest = max(largest,c.coeffs[0]) for c in R_x.coeffs
        largest = [max()for c in R_x.coeffs]

        """
        for i in R_x.coeffs:
            c = i.coeffs[0]
            largest = max(largest,c)

        """


        if largest == 0:
            print("Syndromes all ZERO -> No errors detected")
            return coeffs
        else:
            print("There are some errors you need to deal with")


    def find_error_locations(self,synds):
        """
        Can only correct upto t errors where t = (n - k) / 2
        """

        # Step 1 : Find the error locator polynomial
        A_x = GF256Poly([1])
        v = len(synds)

        t = len(synds)

        for i in range(t - 1):
            alpha_i = GF256.pow(2,i)

            factor = GF256Poly([1,alpha_i])

            #print(factor * synds[i])

            A_x *= factor

        """
        for i in range(v):
            d = synds[i]
            X_i = GF256.pow(2,i)

            factor = GF256Poly([1,X_i])


            A_x *= factor

        # Step 2 : Find the error locations(Chien search)
        for i in range(v):
            continue

        # Step 3 : Calculate the error magnitudes(Forney Algo)
        for i in range(v):
            continue

        """

        return A_x


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
    def inv(cls, a):
        if a == 0: raise ZeroDivisionError("multiplicative inverse of zero is undefined in galois field")
        return cls.exp[255 - cls.log[a]]

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

    def eval(self,x_val):

        res = 0

        for c in self.coeffs:

            res = GF256.add(GF256.mul(res,x_val),c)

        return res


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
    msg = "This is a test message"

    kodek = RSCodec(8,4)

    C_x = kodek.encode(msg)

    #print("C_x",C_x)
    encoded_msg = bytes(C_x.coeffs)

    corrupted_msg = bytearray(encoded_msg)
    corrupted_msg[0] ^= 0xFF

    syndromes = kodek.decode(corrupted_msg)

    res1 = kodek.berlekamp_massey(syndromes)
    print("res1",res1)

    res2 = kodek.berlekamp_massey_second(syndromes)
    print("res2",res2)












