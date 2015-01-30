from .polynomial import Polynomial


def modulo(f, g):
    assert f.field == g.field
    assert g
    field = f.field
    if f.degree < g.degree:
        return f

    g = g * Polynomial(field, [field.inv(g.coeffs[-1])])
    m = g.degree + 1

    r = f.coeffs[:]

    for i in range(len(r), -1, -1):
        coeffs = r[i-m:i]
        if len(coeffs) < m:
            continue
        h = Polynomial(field, [coeffs[-1]])
        p = Polynomial(field, coeffs) - h * g

        assert p.degree < g.degree
        r[i-m:i] = p.coeffs

    return Polynomial(field, r)


def divmod_slow(f, g):
    assert f.field == g.field
    assert f.degree >= g.degree
    assert g

    field = f.field
    r = f.copy()
    factor = field.inv(g.coeffs[-1])
    g = g * Polynomial(field, [factor])

    q = Polynomial(field, [0])
    m = len(f.coeffs)
    n = len(g.coeffs) - 1

    for i in range(m-1, n-1, -1):
        if i < len(r.coeffs):
            h = Polynomial(field, [r.coeffs[i]], power=i-n)
            q = q + h
            r = r - h * g

    assert r + q * g == f

    q = q * Polynomial(field, [factor])
    return q, r


def euclid(a, b):
    assert a.field == b.field
    assert a.degree > b.degree
    P = lambda *coeffs: Polynomial(a.field, coeffs)
    r = [a, b]
    s = [P(1), P(0)]
    t = [P(0), P(1)]
    q = []
    while True:
        assert s[-1] * a + t[-1] * b == r[-1]
        yield r[-1], (s[-1], t[-1])
        q_, r_ = divmod_slow(r[-2], r[-1])
        if not r_:
            break
        q.append(q_)
        r.append(r_)
        s.append(s[-2] - q_ * s[-1])
        t.append(t[-2] - q_ * t[-1])


def gcd(a, b):
    return list(euclid(a, b))[-1][0]
