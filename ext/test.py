import random
import codec

def test():
    c = codec.Codec(10, './librs.so')
    r = random.Random(0)
    for i in range(1000):
        x = bytearray([r.getrandbits(8) for _ in range(c.max_msg_len)])
        y = c.encode(x)
        for _ in range(c.distance // 2):
            i = r.randrange(len(y))
            err = r.getrandbits(8)
            y[i] = err
            z = c.decode(y)
            assert list(z) == list(x)
