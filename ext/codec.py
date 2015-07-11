import ctypes
import logging

log = logging.getLogger(__name__)

class Codec(object):

    def __init__(self, distance, libname):
        self.lib = ctypes.cdll.LoadLibrary(libname)
        self.lib.field_alloc.restype = ctypes.c_void_p
        self.field = ctypes.c_void_p(self.lib.field_alloc())
        if not self.field:
            raise MemoryError()

        g = sum(1 << i for i in [16, 14, 12, 1, 0])  # x^16 + x^14 + x^12 + x + 1
        assert self.lib.rs_field(self.field, g) == 0

        self.max_encoded = self.lib.field_size() - 1
        self.max_msg_len = self.max_encoded - distance
        self.distance = distance

        class polynomial_t(ctypes.Structure):
            _fields_ = [
                ('_length', ctypes.c_uint),
                ('_coeffs', ctypes.c_uint * self.lib.field_size()),
                ('_field', ctypes.c_void_p),
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

        self.lib.poly_alloc.restype = ctypes.POINTER(polynomial_t)
        self.generator = self.lib.poly_alloc()
        log.info('generator: %s', self.generator.contents)
        if not self.generator:
            raise MemoryError()

        log.debug('field: %x', self.field.value)
        log.debug('generator: %x', ctypes.addressof(self.generator.contents))
        assert self.lib.rs_generator(self.field, distance, self.generator) == 0
        log.info('generator: %s', self.generator.contents)

    def __del__(self):
        self.lib.poly_free(self.generator)
        self.generator = None
        self.lib.field_free(self.field)
        self.field = None
        self.lib = None

    def encode(self, msg):
        if len(msg) > self.max_msg_len:
            raise ValueError('too long message to encode')

        p = self.lib.poly_alloc(self.field)
        if not p:
            raise MemoryError()

        p.contents.coeffs = msg
        assert self.lib.rs_encode(self.generator, p) == 0
        encoded = list(p.contents.coeffs)
        self.lib.poly_free(p)
        return encoded

    def decode(self, msg):
        if len(msg) > self.max_encoded:
            raise ValueError('too long message to decode')

        p = self.lib.poly_alloc(self.field)
        if not p:
            raise MemoryError()

        p.contents.coeffs = msg
        errors = self.lib.rs_decode(self.generator, p)
        if errors < 0:
            raise ValueError('too many errors')
        decoded = list(p.contents.coeffs)
        self.lib.poly_free(p)
        return decoded
