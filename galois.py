class Field(object):
    def __init__(self, generator, p):
        for i in range(2, int(p ** 0.5)+1):
            assert p % i != 0

        generator = list(generator)
        n = len(generator) - 1
        q = p ** n
        f = [1] + [0] * (n - 1)

        elements = []
        values = []
        to_poly = {}
        from_poly = {}
        for i in range(q):
            element = tuple(f)
            elements.append(element)
            value = sum(c * (p ** i) for i, c in enumerate(element))
            to_poly[value] = element
            from_poly[element] = value
            values.append(value)

            f = [0] + f  # f'[x] = x * f[x]
            c = f[-1]
            for i, g_i in enumerate(generator):
                f[i] = (f[i] - g_i * c) % p
            assert f.pop(-1) == 0

        assert values[0] == values.pop(-1)  # verify circularity
        assert len(set(values)) == q - 1  # verify uniqueness

        zero_element = tuple([0] * n)
        to_poly[0] = zero_element
        from_poly[zero_element] = 0

        self.exp = dict(enumerate(values + values))
        self.log = dict((v, i) for i, v in enumerate(values))

        self.n = n
        self.p = p
        self.q = q
        self.generator = generator
        self.from_poly = from_poly
        self.to_poly = to_poly

    def __eq__(self, other):
        return self.p == other.p and self.generator == other.generator

    def __repr__(self):
        return 'GF({0}^{1})'.format(self.p, self.n)

    def add(self, x, y):
        pairs = zip(self.to_poly[x], self.to_poly[y])
        return self.from_poly[tuple((x + y) % self.p for x, y in pairs)]

    def sub(self, x, y):
        pairs = zip(self.to_poly[x], self.to_poly[y])
        return self.from_poly[tuple((x - y) % self.p for x, y in pairs)]

    def times(self, x, n):
        return self.from_poly[tuple((c * n) % self.p for c in self.to_poly[x])]

    def mul(self, x, y):
        if x == 0 or y == 0:
            return 0
        return self.exp[self.log[x] + self.log[y]]

    def div(self, x, y):
        return self.mul(x, self.inv(y))

    def inv(self, x):
        if x == 0:
            raise ZeroDivisionError()
        order = self.q - 1
        return self.exp[order - self.log[x]]


class Polynomial(object):

    __slots__ = ['field', 'coeffs']

    def __init__(self, field, coeffs, power=0):
        self.field = field
        self.coeffs = _trim([0] * power + list(coeffs))

    def __nonzero__(self):
        return self.coeffs != [0]

    @property
    def degree(self):
        return len(self.coeffs) - 1

    def copy(self):
        return Polynomial(self.field, list(self.coeffs))

    def __repr__(self):
        return '<{0} {1}>'.format(self.field, self.coeffs)

    def __eq__(self, other):
        return self.field == other.field and self.coeffs == other.coeffs

    def __add__(self, other):
        assert self.field == other.field
        f, g = _pad(self.coeffs, other.coeffs)
        coeffs = [self.field.add(pair[0], pair[1]) for pair in zip(f, g)]
        return Polynomial(self.field, coeffs)

    def __sub__(self, other):
        assert self.field == other.field
        f, g = _pad(self.coeffs, other.coeffs)
        coeffs = [self.field.sub(pair[0], pair[1]) for pair in zip(f, g)]
        return Polynomial(self.field, coeffs)

    def eval(self, x):
        result = 0
        for c in reversed(self.coeffs):  # from higher to lower powers
            result = self.field.mul(result, x)
            result = self.field.add(result, c)
        return result

    def deriv(self):
        coeffs = [self.field.times(c, i) for i, c in enumerate(self.coeffs)]
        return Polynomial(self.field, coeffs[1:])

    def __mul__(self, other):
        result = [0] * (len(self.coeffs) + len(other.coeffs))
        for i, f_i in enumerate(self.coeffs):
            for j, g_j in enumerate(other.coeffs):
                k = i + j
                delta = self.field.mul(f_i, g_j)
                result[k] = self.field.add(result[k], delta)
        return Polynomial(self.field, result)

    def normalize(self):
        factor = self.field.inv(self.coeffs[-1])
        return self * Polynomial(self.field, [factor])


def _pad(f, g):
    n = max(len(f), len(g))
    f = f + [0] * (n - len(f))
    g = g + [0] * (n - len(g))
    return f, g


def _trim(f):
    i = 1
    for i in range(len(f), 0, -1):
        if f[i-1]:
            break
    f = f[:i]
    if not f:
        f = [0]
    return f


def divide(f, g):
    assert f.field == g.field
    assert f.degree >= g.degree
    assert g

    field = f.field
    r = f.copy()
    factor = field.inv(g.coeffs[-1])

    q = Polynomial(field, [0])
    m = len(f.coeffs)
    n = len(g.coeffs) - 1

    for i in range(m-1, n-1, -1):
        if i < len(r.coeffs):
            c = field.mul(r.coeffs[i], factor)
            h = Polynomial(field, [c], power=i-n)
            q = q + h
            r = r - h * g

    assert r + q * g == f
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
        q_, r_ = divide(r[-2], r[-1])
        if not r_:
            break
        q.append(q_)
        r.append(r_)
        s.append(s[-2] - q_ * s[-1])
        t.append(t[-2] - q_ * t[-1])


def gcd(a, b):
    return list(euclid(a, b))[-1][0]
