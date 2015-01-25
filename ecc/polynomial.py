class Polynomial(object):

    __slots__ = ['field', 'coeffs']

    def __init__(self, field, coeffs, power=0):
        self.field = field
        self.coeffs = _trim([0] * power + list(coeffs))

    def __nonzero__(self):
        return self.coeffs != [0]

    __bool__ = __nonzero__

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
