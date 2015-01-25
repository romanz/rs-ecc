from ecc.galois import Field
from ecc.polynomial import Polynomial, _trim
from ecc.algorithms import divmod_slow, modulo, euclid


GF = Field(reversed([1, 0, 0, 0, 1, 1, 1, 0, 1]), p=2)


def test_logs():
    assert GF.log[7] == 198
    assert GF.log[3] == 25

    assert GF.exp[198] == 7
    assert GF.exp[25] == 3


def test_inv():
    x = 125
    assert GF.mul(x, GF.inv(x)) == 1


def test_add():
    p1 = (Polynomial(GF, [5, 3]))
    p2 = (Polynomial(GF, [4, 1, 9]))

    assert p1 + p2 == Polynomial(GF, [1, 2, 9])
    assert p1 - p2 == Polynomial(GF, [1, 2, 9])


def test_mul():
    p1 = Polynomial(GF, [1, 0, 2])
    p2 = Polynomial(GF, [1, 3])
    assert p1 * p2 == Polynomial(GF, [1, 3, 2, 6])


def test_div():
    p1 = Polynomial(GF, [1, 0, 0, 2])
    p2 = Polynomial(GF, [1])
    assert divmod_slow(p1, p2) == (p1, Polynomial(GF, [0]))


def test_mod():
    p1 = Polynomial(GF, [8, 3, 5, 2, 4, 6, 9, 7, 1])
    p2 = Polynomial(GF, [10, 11, 12])
    assert divmod_slow(p1, p2)[1] == modulo(p1, p2)


def test_trim():
    assert _trim([]) == [0]
    assert _trim([0]) == [0]
    assert _trim([1]) == [1]
    assert _trim([1, 0]) == [1]
    assert _trim([1, 2]) == [1, 2]
    assert _trim([1, 0, 0, 2, 0]) == [1, 0, 0, 2]


def test_euclid():
    h = Polynomial(GF, [1, 1])
    f1 = Polynomial(GF, [1, 0, 2]) * h
    f2 = Polynomial(GF, [1, 3]) * h
    gcd, (s, t) = list(euclid(f1, f2))[-1]
    assert s*f1 + t*f2 == gcd

    zero = Polynomial(GF, [0])
    assert divmod_slow(f1, gcd)[1] == zero
    assert divmod_slow(f2, gcd)[1] == zero
