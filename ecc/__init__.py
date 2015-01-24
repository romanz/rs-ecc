import galois
import reed_solomon


class Codec(object):

    def __init__(self, distance, field=None):
        if field is None:
            field = galois.Field(reversed([1, 0, 0, 0, 1, 1, 1, 0, 1]), p=2)

        self.max_msg_len = field.q - 1
        self.distance = distance
        self.generator = reed_solomon.generator(field=field, n=distance)

    def encode(self, msg):
        encoded = reed_solomon.encode(msg=bytearray(msg), gen=self.generator)
        return bytearray(encoded)

    def decode(self, msg):
        decoded = reed_solomon.decode(msg=bytearray(msg), gen=self.generator)
        return bytearray(decoded)
