class Field(object):
    def __init__(self, generator, p):
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
        coeffs = tuple((x + y) % self.p for x, y in pairs)
        return self.from_poly[coeffs]

    def sub(self, x, y):
        pairs = zip(self.to_poly[x], self.to_poly[y])
        coeffs = tuple((x - y) % self.p for x, y in pairs)
        return self.from_poly[coeffs]

    def mul(self, x, y):
        if x == 0 or y == 0:
            return 0
        return self.exp[self.log[x] + self.log[y]]


class Polynomial(object):

    __slots__ = ['field', 'coeffs']

    def __init__(self, field, coeffs):
        self.field = field
        self.coeffs = list(coeffs)

    def __repr__(self):
        return '<{0} {1}>'.format(self.field, self.coeffs)


def trim(f):
    i = 1
    for i in range(len(f), 0, -1):
        if f[i-1]:
            break
    return f[:i]


def add(f, g):
    assert f.field == g.field
    field = f.field
    f = f.coeffs
    g = g.coeffs
    n = max(len(f), len(g))
    f = f + [0] * (n - len(f))
    g = g + [0] * (n - len(g))
    coeffs = [field.add(pair[0], pair[1]) for pair in zip(f, g)]
    return trim(coeffs)
