from ecc import Codec
import random


def test():
    d = 6
    c = Codec(d)
    m = c.max_msg_len

    r = random.Random(0)

    for i in range(10):
        msg = bytearray(r.getrandbits(8) for _ in range(m))
        enc_msg = c.encode(msg)

        for i in range(d // 2):
            enc_msg[r.randrange(len(msg))] = r.getrandbits(8)

        dec_msg = c.decode(enc_msg)
        assert dec_msg == msg
