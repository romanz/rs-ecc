from . import galois
from . import reed_solomon


class Codec(object):

    def __init__(self, distance, field=None):
        if field is None:
            field = galois.Field(reversed([1, 0, 0, 0, 1, 1, 1, 0, 1]), p=2)

        self.max_encoded = field.q - 1
        self.max_msg_len = self.max_encoded - distance
        self.distance = distance
        self.generator = reed_solomon.generator(field=field, n=distance)

    def encode(self, msg):
        msg = list(msg)
        if len(msg) > self.max_msg_len:
            raise ValueError('too long message to encode')
        encoded = reed_solomon.encode(msg=msg, gen=self.generator)
        return bytearray(encoded)

    def decode(self, msg):
        msg = list(msg)
        if len(msg) > self.max_encoded:
            raise ValueError('too long message to decode')
        decoded = reed_solomon.decode(msg=msg, gen=self.generator)
        return bytearray(decoded)
