import galois as g
import reed_solomon as rs
import pytest

GF = g.Field(reversed([1, 0, 0, 0, 1, 1, 1, 0, 1]), p=2)
gen = rs.generator(GF, 4)


def test_gen():
    assert gen == g.Polynomial(GF, [0x40, 0x78, 0x36, 0x0f, 1])


def test_encode():
    data = [0x56, 0x34, 0x12]
    r = rs.encode(data, gen)
    assert r == [0xd9, 0x78, 0xe6, 0x37] + data

    msg = [
        236, 112, 150, 198, 198, 150, 38, 39, 6, 50, 23, 118, 71, 117, 210, 64
    ]
    f = rs.encode(msg, rs.generator(GF, 10))
    assert f == [224, 75, 253, 239, 175, 107, 19, 144, 42, 188] + msg


def test_syndrome():
    d = 10
    gen = rs.generator(GF, d)
    orig = [
        224, 75, 253, 239, 175, 107, 19, 144, 42, 188,
        236, 112, 150, 198, 198, 150, 38, 39, 6, 50, 23, 118, 71, 117, 210, 64,
    ]
    msg = orig[:]
    assert rs.syndrome(msg, gen) == [0] * gen.degree

    index = len(msg) - 1
    msg[index] = 0
    synd = rs.syndrome(msg, gen)
    assert synd == [64, 192, 93, 231, 52, 92, 228, 49, 83, 245]

    indices = [2, 8, 15, 7, 0]
    msg = orig[:]
    for i in indices:
        msg[i] = i
    assert rs.decode(msg, gen) == orig[d:]

    indices = [0, 1, 2, 23, 24, 25]
    msg = orig[:]
    for i in indices:
        msg[i] = i

    with pytest.raises(ValueError):
        rs.decode(msg, gen)

    with pytest.raises(ValueError):
        rs.decode(range(len(msg)), gen)
