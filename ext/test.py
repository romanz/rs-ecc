from ecc import codec

def test_py():
    c = codec.Codec(10)
    x = [0] * 210
    errors = {137: 211L, 155: 255L}
    for k, v in errors.items():
        x[k] = v
    assert c.decode(x) == bytearray([0] * 200)


import ctypes

class Codec(object):

    def __init__(self, distance, libname='./librs.so'):
        self.lib = ctypes.cdll.LoadLibrary(libname)
        self.field = self.lib.field_alloc()
        if not self.field:
            raise MemoryError()

        assert self.lib.rs_field(self.field, 0b100011101) == 0

        self.max_msg_len = self.lib.field_size() - 1 - distance
        self.max_encoded = self.lib.field_size() - 1
        self.distance = distance

        class polynomial_t(ctypes.Structure):
            _fields_ = [
                ('_length', ctypes.c_uint),
                ('_coeffs', ctypes.c_uint * self.lib.field_size()),
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

        self.lib.polynomial_alloc.restype = ctypes.POINTER(polynomial_t)
        self.generator = self.lib.polynomial_alloc()
        if not self.generator:
            raise MemoryError()

        assert self.lib.rs_generator(self.field, distance, self.generator) == 0

    def __del__(self):
        self.lib.polynomial_free(self.generator)
        self.generator = None
        self.lib.field_free(self.field)
        self.field = None
        self.lib = None

    def encode(self, msg):
        if len(msg) > self.max_msg_len:
            raise ValueError('too long message to encode')

        p = self.lib.polynomial_alloc(self.field)
        if not p:
            raise MemoryError()

        p.contents.coeffs = msg
        assert self.lib.rs_encode(self.generator, p) == 0
        encoded = bytearray(p.contents.coeffs)
        self.lib.polynomial_free(p)
        return encoded

    def decode(self, msg):
        if len(msg) > self.max_encoded:
            raise ValueError('too long message to decode')

        p = self.lib.polynomial_alloc(self.field)
        if not p:
            raise MemoryError()

        p.contents.coeffs = msg
        errors = self.lib.rs_decode(self.generator, p)
        if errors < 0:
            raise ValueError('too many errors')
        decoded = bytearray(p.contents.coeffs)
        self.lib.polynomial_free(p)
        return decoded

import random

def test_c():
    c = Codec(64)
    r = random.Random(0)
    for i in range(1000):
        x = bytearray([r.getrandbits(8) for _ in range(c.max_msg_len)])
        y = c.encode(x)
        errors = {}
        for _ in range(c.distance // 2):
            i = r.randrange(len(y))
            errors[i] = r.getrandbits(8)
        for k, v in errors.items():
            y[k] = v
        z = c.decode(y)
        assert list(z) == list(x)
