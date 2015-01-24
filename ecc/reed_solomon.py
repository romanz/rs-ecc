import galois


def generator(field, n):
    ''' g[x] = (x - g^0)(x - g^1)...(x - g^(n-1)) '''
    g = galois.Polynomial(field, [1])
    for i in range(n):
        c = field.sub(0, field.exp[i])  # -g^i
        g = g * galois.Polynomial(field, [c, 1])
    assert g.degree == n
    return g


def encode(msg, gen):
    padded = galois.Polynomial(gen.field, list(msg), power=gen.degree)
    q, r = galois.divide(padded, gen)
    assert q * gen + r == padded
    return r.coeffs + msg


def decode(msg, gen):
    assert len(msg) < gen.field.q  # must be <= (p^n-1)
    synd = syndrome(msg, gen)
    if any(synd):
        locator, evaluator = solve(gen.field, synd)
        indices = search(locator, range(len(msg)))
        if indices is None:
            raise ValueError('Cannot find errors')
        for j, delta in correct(indices, locator, evaluator):
            msg[j] = gen.field.add(msg[j], delta)

    prefix_len = gen.degree
    return msg[prefix_len:]


def syndrome(msg, gen):
    field = gen.field
    degree = gen.degree
    p = galois.Polynomial(field, msg)
    return [p.eval(field.exp[i]) for i in range(degree)]


def solve(field, synd):
    S = galois.Polynomial(field, list(synd))
    a = galois.Polynomial(field, [1], power=S.degree+1)
    b = S

    k = a.degree // 2
    for r, (s, t) in galois.euclid(a, b):
        if t.degree <= k and r.degree < k:
            break
    else:
        raise AssertionError()

    c = field.inv(t.coeffs[0])
    c = galois.Polynomial(field, [c])
    locator = c * t
    evaluator = c * r

    assert locator.eval(0) == 1
    assert galois.gcd(locator, evaluator).normalize().coeffs == [1]
    # Locator[x] Syndrome[x] = Evaluator[x] {mod x^(d-1)}
    assert galois.divide(locator * S, a)[1] == evaluator
    return locator, evaluator


def search(p, indices):
    field = p.field
    # q[x] == 0 iff p[1/x] == 0
    q = galois.Polynomial(field, reversed(p.coeffs))
    indices = [i for i in indices if q.eval(field.exp[i]) == 0]
    if len(indices) == p.degree:
        return indices


def correct(indices, locator, evaluator):
    locator_deriv = locator.deriv()
    field = locator.field

    for j in indices:
        alpha = field.exp[j]
        inv_alpha = field.inv(alpha)
        num = field.mul(evaluator.eval(inv_alpha), alpha)
        den = locator_deriv.eval(inv_alpha)
        yield j, field.div(num, den)
