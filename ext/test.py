import ctypes
lib = ctypes.cdll.LoadLibrary('./librs.so')

class polynomial_t(ctypes.Structure):
    _fields_ = [
        ('_length', ctypes.c_uint),
        ('_coeffs', ctypes.c_uint * lib.field_size()),
    ]

    def __repr__(self):
        return repr(self.coeffs)

    @property
    def coeffs(self):
        return [int(c) for c in self._coeffs[:self._length]]

    @coeffs.setter
    def coeffs(self, values):
        self._length = len(values)
        self._coeffs[:self._length] = [int(c) for c in values]


lib.polynomial_alloc.restype = ctypes.POINTER(polynomial_t)

F = lib.field_alloc()
assert F
lib.rs_field(F, 0b100011101)

def test():
    gen = lib.polynomial_alloc(F)
    assert gen

    lib.rs_generator(F, 4, gen)
    assert gen.contents.coeffs == [0x40, 0x78, 0x36, 0x0f, 1]

    msg = [0x56, 0x34, 0x12]
    p = lib.polynomial_alloc(F)
    p.contents.coeffs = msg

    lib.rs_encode(gen, p)
    assert p.contents.coeffs == [0xd9, 0x78, 0xe6, 0x37] + msg

    d = 10
    lib.rs_generator(F, d, gen)
    msg = [236, 112, 150, 198, 198, 150, 38, 39, 6, 50, 23, 118, 71, 117, 210, 64]
    p.contents.coeffs = msg
    lib.rs_encode(gen, p)
    assert p.contents.coeffs == [224, 75, 253, 239, 175, 107, 19, 144, 42, 188] + msg
