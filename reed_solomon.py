import galois


def generator(field, n):
    ''' g[x] = (x - g^0)(x - g^1)...(x - g^(n-1)) '''
    g = galois.Polynomial(field, [1])
    for i in range(n):
        c = field.sub(0, field.exp[i])  # -g^i
        g = g * galois.Polynomial(field, [c, 1])
    return g


def encode(msg, gen):
    padded = galois.Polynomial(gen.field, list(msg), power=gen.order)
    q, r = galois.divide(padded, gen)
    assert q * gen + r == padded
    return msg + r.coeffs
