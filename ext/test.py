import ctypes
lib = ctypes.cdll.LoadLibrary('./librs.so')

class polynomial_t(ctypes.Structure):
    _fields_ = [
        ('_length', ctypes.c_uint),
        ('_coeffs', ctypes.c_uint * lib.max_size()),
    ]

    @property
    def coeffs(self):
        return [int(c) for c in self._coeffs[:self._length]]

lib.polynomial_alloc.restype = ctypes.POINTER(polynomial_t)

F = lib.field_alloc()
gen = lib.polynomial_alloc(F)

assert F
assert gen

def test_gen():
    lib.rs_init_field(F, 0x11d)
    lib.rs_generator(F, 4, gen)
    assert gen.contents.coeffs == [0x40, 0x78, 0x36, 0x0f, 1]
