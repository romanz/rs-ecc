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

        to_poly[0] = tuple([0] * n)
        from_poly[to_poly[0]] = 0

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
