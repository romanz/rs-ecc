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
