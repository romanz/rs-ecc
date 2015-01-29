from . import galois
from . import reed_solomon

class Codec(object):

    def __init__(self, distance, field=None):
        if field is None:
            indices = [16, 14, 12, 1, 0]  # x^16 + x^14 + x^12 + x + 1
            coeffs = [0] * (max(indices) + 1)
            for i in indices:
                coeffs[i] = 1

            field = galois.Field(reversed(coeffs), p=2)

        self.max_encoded = field.q - 1
        self.max_msg_len = self.max_encoded - distance
        self.distance = distance
        self.generator = reed_solomon.generator(field=field, n=distance)

    def encode(self, msg):
        msg = list(msg)
        if len(msg) > self.max_msg_len:
            raise ValueError('too long message to encode')
        encoded = reed_solomon.encode(msg=msg, gen=self.generator)
        return encoded

    def decode(self, msg):
        msg = list(msg)
        if len(msg) > self.max_encoded:
            raise ValueError('too long message to decode')
        decoded = reed_solomon.decode(msg=msg, gen=self.generator)
        return decoded
