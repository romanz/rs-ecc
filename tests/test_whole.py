from ecc.codec import Codec
import random


def test_codec(distance=4, iters=10, seed=0):
    c = Codec(distance)
    max_errors = distance // 2

    r = random.Random(seed)

    width = 16
    for i in range(iters):
        msg = [r.getrandbits(width) for _ in range(2048)]
        enc_msg = c.encode(msg)
        assert c.decode(enc_msg) == msg

        for j in range(1, max_errors + 1):
            value = r.getrandbits(width)
            index = r.randrange(len(enc_msg))
            enc_msg[index] = value
            assert c.decode(enc_msg) == msg
