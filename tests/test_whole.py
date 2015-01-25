from ecc.codec import Codec
import random


def test_codec(distance=4, iters=100, seed=0):
    c = Codec(distance)
    max_errors = distance // 2

    r = random.Random(seed)

    for i in range(iters):
        msg = bytearray(r.getrandbits(8) for _ in range(c.max_msg_len))
        enc_msg = c.encode(msg)
        assert c.decode(enc_msg) == msg

        for j in range(max_errors):
            value = r.getrandbits(8)
            index = r.randrange(len(enc_msg))
            enc_msg[index] = value
            assert c.decode(enc_msg) == msg
