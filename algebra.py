from friC import friC
from ctypes import c_uint64, c_uint32, c_void_p, c_double, POINTER, c_bool, cast

class FP16:
    def __init__(self) -> None:
        self.lib = friC.lib
        self.lib.fp_ini.argtypes = (c_uint32,)

        self.lib.fp_ini(2**16 + 1)

        self.fp_t = c_uint32
        self.fp2_t = self.fp_t * 2
        self.fp4_t = self.fp2_t * 2
        self.fp8_t = self.fp4_t * 2
        self.fp16_t = self.fp8_t * 2

        self.zero = self.new(0)
        
        # arith
        self.lib.fp16_mul.argtypes = (self.fp16_t, self.fp16_t, self.fp16_t)
        self.lib.fp16_add.argtypes = (self.fp16_t, self.fp16_t, self.fp16_t)
        self.lib.fp16_sub.argtypes = (self.fp16_t, self.fp16_t, self.fp16_t)
        self.lib.fp16_inv.argtypes = (self.fp16_t, self.fp16_t)
        # cmp
        self.lib.fp16_cmp.argtypes = (self.fp16_t, self.fp16_t)
        self.lib.fp16_cmp.restype = c_uint32

        # io
        self.lib.fp16_in.argtypes = (self.fp16_t, POINTER(c_uint32))
        self.lib.fp16_out.argtypes = (POINTER(c_uint32), self.fp16_t)
        self.lib.fp16_get_powers.argtypes = (POINTER(self.fp16_t), self.fp16_t, c_uint32)
    
    def new(self, v):
        if(type(v) == int):
            v = [v] + 15*[0] 
        res = self.fp16_t()
        res[0][0][0][0] = v[0]
        res[0][0][0][1] = v[1]
        res[0][0][1][0] = v[2]
        res[0][0][1][1] = v[3]
        res[0][1][0][0] = v[4]
        res[0][1][0][1] = v[5]
        res[0][1][1][0] = v[6]
        res[0][1][1][1] = v[7]
        res[1][0][0][0] = v[8]
        res[1][0][0][1] = v[9]
        res[1][0][1][0] = v[10]
        res[1][0][1][1] = v[11]
        res[1][1][0][0] = v[12]
        res[1][1][0][1] = v[13]
        res[1][1][1][0] = v[14]
        res[1][1][1][1] = v[15]
        return res
    
    def get(self, x):
        v = [0]*16
        v[0] = x[0][0][0][0]
        v[1] = x[0][0][0][1]
        v[2] = x[0][0][1][0]
        v[3] = x[0][0][1][1]
        v[4] = x[0][1][0][0]
        v[5] = x[0][1][0][1]
        v[6] = x[0][1][1][0]
        v[7] = x[0][1][1][1]
        v[8] = x[1][0][0][0]
        v[9] = x[1][0][0][1]
        v[10] = x[1][0][1][0]
        v[11] = x[1][0][1][1]
        v[12] = x[1][1][0][0]
        v[13] = x[1][1][0][1]
        v[14] = x[1][1][1][0]
        v[15] = x[1][1][1][1]
        return v
    
fp16 = FP16()

class FieldElement:
    def __init__( self, value=None):
        if(not value):
            self.value = fp16.fp16_t()
        elif type(value) == fp16.fp16_t:
            self.value = value
        else:
            self.value = fp16.new(value)

    def __getstate__(self):
        return fp16.get(self.value)

    def __setstate__(self, state):
        self.value = state

    def __add__( self, right ): 
        res = FieldElement()
        fp16.lib.fp16_add(res.value, self.value, right.value)
        return res

    def __mul__( self, right ):
        res = FieldElement()
        fp16.lib.fp16_mul(res.value, self.value, right.value)
        return res

    def __sub__( self, right ):
        res = FieldElement()
        fp16.lib.fp16_sub(res.value, self.value, right.value)
        return res

    def __truediv__( self, right ):
        assert(False), "not implemented"
        assert(not right.is_zero()), "divide by zero"
        return FieldElement(self.value/right.value, self.field)

    def __neg__( self ):
        res = FieldElement()
        fp16.lib.fp16_sub(res.value, fp16.zero, self.value)
        return res

    def inverse( self ):
        res = FieldElement()
        fp16.lib.fp16_inv(res.value, self.value)
        return res

    # modular exponentiation -- be sure to encapsulate in parentheses!
    def __xor__( self, exponent ):
        assert(False), "use **"

    def __pow__(self, exponent):
        if(exponent == -1):
            return self.inverse()
        if(exponent == 2):
            res = FieldElement()
            fp16.lib.fp16_mul(res.value, self.value, self.value)
            return res
        assert(False), "pow is not implemented"

    def __eq__( self, other ):
        return fp16.get(self.value) == fp16.get(other.value)

    def __neq__( self, other ):
        return fp16.get(self.value) != fp16.get(other.value)

    def __str__( self ):
        return str(fp16.get(self.value))

    def __repr__(self):
        return str(fp16.get(self.value))

    def __bytes__( self ):
        assert(False), "not implemented"
        fp16.get(self.value)
        return int(self.value.polynomial().change_ring(ZZ)(self.field.f.base().order())).to_bytes(self.field.byte_len, "little")

    def is_zero( self ):
        return fp16.get(self.value) == [0]*16
        if self.value == 0:
            return True
        else:
            return False

import math

class Field:
    def __init__( self, p, D):
        self.p = p
        self.D = D
        self.byte_len = int(math.ceil(math.log2(self.p)/8)*self.D)
        # print("Setting up field with p =", self.p, "D =", self.D, "hash byte length =", self.byte_len)

    def get_powers(self, x, size):
        res = (fp16.fp16_t*(size))()
        fp16.lib.fp16_get_powers(res, x.value, size)
        return [FieldElement(i) for i in res]

    def __call__(self, x):
        return FieldElement(x)

    def to_list(self, x):
        return fp16.get(x.value)
        
    def zero( self ):
        return FieldElement(0)

    def one( self ):
        return FieldElement(1)

    def optFP16():
        p = 2**16 + 1
        D = 16
        return Field(p, D)


    # def main(p = 1 + 407 * ( 1 << 119 )):
    #     p = 2**16 + 1
    #     f = GF(p**16)
    #     return Field(f)
    
    # def towerFp16():
    #     def conv_fp16_tower(vec):
    #         pos = [0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15]
    #         return [vec[i] for i in pos]

    #     p = 2**16+1
    #     Fp = GF(p)
    #     R = PolynomialRing(Fp, "x")
    #     x = R.gen()
    #     Fp2 = Fp.extension(x**2 - 3, "x") 
    #     i = Fp2.gen()
    #     S = PolynomialRing(Fp2, "y")
    #     y = S.gen()
    #     Fp4 = Fp2.extension(y**2 - i, "y") 
    #     j = Fp4.gen()
    #     T = PolynomialRing(Fp4, "z")
    #     z = T.gen()
    #     Fp8 = Fp4.extension(z**2 - j, "z")
    #     k = Fp8.gen()
    #     U = PolynomialRing(Fp8, "w")
    #     w = U.gen()
    #     Fp16 = Fp8.extension(w**2 - k, "k")
    #     return Field(Fp16, conv=conv_fp16_tower)

    def primitive_nth_root( self, n ):
        rous = [1, 65536, 65281, 4096, 64, 65529, 8224, 13987, 282, 15028, 19139, 61869, 54449, 6561, 81, 9, 3]
        return rous[int(math.log2(n))]
            
    def sample( self, byte_array ):
        pointer = friC.sample_field_array(byte_array, self.D, self.p)
        val = cast(pointer, POINTER(c_uint64*16))
        return FieldElement(list(val.contents))
        # requires more than 256 bits because p^d > 2^256
        a = [byte_array[i]|(byte_array[i+1]<<8)|((byte_array[i + 2]&1)<<16) for i in range(0,16*3,3)]
        a = [int((i/(2**17))*self.p)%self.p for i in a]
        return FieldElement(a)




