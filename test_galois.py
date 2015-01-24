import galois as g

GF = g.Field(reversed([1, 0, 0, 0, 1, 1, 1, 0, 1]), p=2)


def test_logs():
    assert GF.log[7] == 198
    assert GF.log[3] == 25

    assert GF.exp[198] == 7
    assert GF.exp[25] == 3


def test_inv():
    x = 125
    assert GF.mul(x, GF.inv(x)) == 1


def test_add():
    p1 = (g.Polynomial(GF, [5, 3]))
    p2 = (g.Polynomial(GF, [4, 1, 9]))

    assert p1 + p2 == g.Polynomial(GF, [1, 2, 9])
    assert p1 - p2 == g.Polynomial(GF, [1, 2, 9])


def test_mul():
    p1 = g.Polynomial(GF, [1, 0, 2])
    p2 = g.Polynomial(GF, [1, 3])
    assert p1 * p2 == g.Polynomial(GF, [1, 3, 2, 6])


def test_div():
    p1 = g.Polynomial(GF, [1, 0, 0, 2])
    p2 = g.Polynomial(GF, [1])
    assert g.divide(p1, p2) == (p1, g.Polynomial(GF, [0]))


def test_trim():
    assert g._trim([]) == [0]
    assert g._trim([0]) == [0]
    assert g._trim([1]) == [1]
    assert g._trim([1, 0]) == [1]
    assert g._trim([1, 2]) == [1, 2]
    assert g._trim([1, 0, 0, 2, 0]) == [1, 0, 0, 2]


def test_euclid():
    h = g.Polynomial(GF, [1, 1])
    f1 = g.Polynomial(GF, [1, 0, 2]) * h
    f2 = g.Polynomial(GF, [1, 3]) * h
    gcd, (s, t) = list(g.euclid(f1, f2))[-1]
    assert s*f1 + t*f2 == gcd

    zero = g.Polynomial(GF, [0])
    assert g.divide(f1, gcd)[1] == zero
    assert g.divide(f2, gcd)[1] == zero
