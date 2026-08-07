def build_generator_poly(nsym):
    pass

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
        pass

    def __mul__(self,other):
        pass

    def divmod(self,divisor):
        pass

    def __repr__(self):
        # todo : implement this next so we can start 
        # testing it immediately

        terms = []
        # idx = 0 char = 1

        for idx,char in enumerate(self.coeffs):
            # degree here is zero
            deg = len(self.coeffs) - 1

            # power here is then zero
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

    msg_bytes = [ord(c) for c in msg]

    M_x = GF256Poly(msg_bytes)
    print(M_x)

