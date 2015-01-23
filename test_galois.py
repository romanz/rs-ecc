import galois as g

def test_add():
    GF = g.Field([1, 0, 0, 0, 1, 1, 1, 0, 1], p=2)
    p1 = (g.Polynomial(GF, [5, 3]))
    p2 = (g.Polynomial(GF, [4, 1]))
    assert g.add(p1, p2) == [1, 2]
