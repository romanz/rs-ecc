import random
import codec as c_codec
from ecc import codec as py_codec

import logging
logging.basicConfig(level=0, format='%(message)-80s %(pathname)s:%(lineno)d')

def test():
    d = 2
    c = c_codec.Codec(d, './librs.so')
    py = py_codec.Codec(d)
    n = 10
    r = random.Random(0)
    for i in range(1000):
        x = [r.getrandbits(8) for i in range(n)]
        y_c = c.encode(x)
        y_py = py.encode(x)
        assert y_c == y_py

        z_py = py.decode(y_py)
        assert z_py == x
        z_c = c.decode(y_c)
        assert z_c == x
        return

        for _ in range(c.distance // 2):
            i = r.randrange(len(y))
            err = r.getrandbits(8)
            y[i] = err
            z = c.decode(y)
            assert list(z) == list(x)
