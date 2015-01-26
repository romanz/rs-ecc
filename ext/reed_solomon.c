#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define N (8)
#define Q (1L << (N))

typedef unsigned int uint;

typedef uint symbol_t;

struct field_t {
	symbol_t exp[Q];
	symbol_t log[Q];
	symbol_t inv[Q];
};

struct polynomial_t {
	uint length;
	symbol_t coeffs[Q];
	const struct field_t *field;
};

uint field_size() {
    return Q;
}

void poly_print(const struct polynomial_t p, const char *fmt) {
    printf("[");
    printf(fmt, p.coeffs[0]);
    for (int i = 1; i < p.length; ++i) {
        printf(", ");
        printf(fmt, p.coeffs[i]);
    }
    printf("]");
}

symbol_t symb_mult(const struct field_t *field, symbol_t x, symbol_t y) {
	if (x == 0 || y == 0) {
		return 0;
	}
	int exponent = field->log[x] + field->log[y];
	return field->exp[exponent % (Q-1)];
}


symbol_t symb_inv(const struct field_t *field, symbol_t x) {
	assert(x);
	return field->inv[x];
}

struct polynomial_t poly_mult(struct polynomial_t f, struct polynomial_t g) {
    assert(f.field);
    assert(f.field == g.field);
    const struct field_t *field = f.field;

    struct polynomial_t result = {f.length + g.length - 1, {0}, field};

    for (int i = 0; i < f.length; ++i) {
        for (int j = 0; j < g.length; ++j) {
            int k = i + j;
            result.coeffs[k] ^= symb_mult(field, f.coeffs[i], g.coeffs[j]);
        }
    }
    return result;
}

struct polynomial_t poly_add(struct polynomial_t f, struct polynomial_t g) {
    assert(f.field);
    assert(f.field == g.field);
    const struct field_t *field = f.field;

    uint length = (f.length > g.length) ? f.length : g.length;
    struct polynomial_t result = {length, {0}, field};

    for (int i = 0; i < length; ++i) {
        result.coeffs[i] = f.coeffs[i] ^ g.coeffs[i];
    }
    return result;
}

struct polynomial_t poly_scale(struct polynomial_t f, symbol_t c) {
    assert(f.field);
    const struct field_t *field = f.field;

    for (int i = 0; i < f.length; ++i) {
        f.coeffs[i] = symb_mult(field, f.coeffs[i], c);
    }
    return f;
}

void poly_divmod(struct polynomial_t f, struct polynomial_t g,
                 struct polynomial_t *q_ptr, struct polynomial_t *r_ptr) {
    assert(f.field);
    assert(f.field == g.field);

    const struct field_t *field = f.field;

    symbol_t factor = symb_inv(field, g.coeffs[g.length - 1]);
    g = poly_scale(g, factor);

    struct polynomial_t q = {0, {0}, field};
    struct polynomial_t r = f;

    for (int i = f.length - 1; i >= g.length - 1; --i) {
        int offset = i - g.length + 1;
        assert(offset >= 0);

        symbol_t c = r.coeffs[i];
        struct polynomial_t g_c = poly_scale(g, c);
        for (int j = 0; j < g.length; ++j) {
            r.coeffs[offset + j] ^= g_c.coeffs[j];
        }
        assert(r.coeffs[i] == 0);
        r.length -= 1;
        q.coeffs[offset] = c;
        q.length += 1;
    }
    if (r_ptr) {
        *r_ptr = r;
    }
    if (q_ptr) {
        *q_ptr = q;
    }
}

symbol_t poly_eval(struct polynomial_t f, symbol_t x) {
    symbol_t res = 0;
    symbol_t z = 1;
    for (int i = 0; i < f.length; ++i) {
        res ^= symb_mult(f.field, f.coeffs[i], z);
        z = symb_mult(f.field, z, x);
    }
    return res;
}

struct polynomial_t syndrome(struct polynomial_t f, uint size) {
    struct polynomial_t res = {size, {0}, f.field};
    for (int i = 0; i < size; ++i) {
        res.coeffs[i] = poly_eval(f, f.field->exp[i]);
    }
    return res;
}

struct field_t *field_alloc() {
	return calloc(1, sizeof(struct field_t));
}

void field_free(struct field_t *f) {
	free(f);
}

struct polynomial_t *polynomial_alloc(const struct field_t *f) {
	struct polynomial_t *p = calloc(1, sizeof(struct polynomial_t));
    if (p) {
        p->field = f;
    }
    return p;
}

void polynomial_free(struct polynomial_t *p) {
	free(p);
}

///////////////////////////////////////////////////////////////////////////////


int rs_field(struct field_t *field, symbol_t primitive_poly) {
	assert(field);
	symbol_t alpha = 1;

	for (int i = 0; i < Q-1; ++i) {
		field->exp[i] = alpha;

		assert(field->log[alpha] == 0);
		field->log[alpha] = i;

		alpha = alpha << 1;
		if (alpha & (1 << N)) {
			alpha = alpha ^ primitive_poly;
		}
	}
	assert(alpha == 1);
	field->exp[Q-1] = 1;

	for (int x = 1; x < Q-1; ++x) {
		uint e = Q - 1 - field->log[x];
		symbol_t y = field->exp[e];
		assert(symb_mult(field, x, y) == 1);
		field->inv[x] = y;
	}
	return 0;
}

int rs_generator(const struct field_t *field, uint distance,
                 struct polynomial_t *gen) {
	assert(field);
	assert(gen);

	struct polynomial_t result = {1, {1}, field};

	for (int i = 0; i < distance; ++i) {
		symbol_t alpha = field->exp[i];
		struct polynomial_t p = {2, {alpha, 1}, field};
		result = poly_mult(result, p);
	}

	*gen = result;
	return 0;
}

int is_zero(const struct polynomial_t p) {
    for (int i = 0; i < p.length; ++i) {
        if (p.coeffs[i]) {
            return 0;
        }
    }
    return 1;
}

int rs_encode(const struct polynomial_t *gen, struct polynomial_t *msg) {
    assert(gen);
    assert(gen->field);
    assert(msg);
    msg->field = gen->field;

	uint offset = gen->length - 1;
    assert(msg->length + offset < Q);

    for (int i = msg->length - 1; i >= 0; --i) {
        msg->coeffs[i + offset] = msg->coeffs[i];
        msg->coeffs[i] = 0;
    }
    msg->length += offset;

    struct polynomial_t res = {0};
    poly_divmod(*msg, *gen, NULL, &res);
    for (int i = 0; i < msg->length; ++i) {
        res.coeffs[i + offset] = msg->coeffs[i + offset];
    }
    res.length = msg->length;
    *msg = res;
    return 0;
}

int key_equation_solver(const struct polynomial_t synd,
                        struct polynomial_t *locator_ptr,
                        struct polynomial_t *evaluator_ptr) {

    uint n = synd.length;
    struct polynomial_t a = {n + 1, {0}, synd.field};
    a.coeffs[n] = 1;
    struct polynomial_t b = synd;
    assert(a.length > b.length);

    struct polynomial_t one = {1, {1}, synd.field};
    struct polynomial_t zero = {1, {0}, synd.field};

    struct polynomial_t r[2] = {a, b};
    struct polynomial_t s[2] = {one, zero};
    struct polynomial_t t[2] = {zero, one};
    struct polynomial_t q_ = zero;
    struct polynomial_t r_ = zero;

    uint k = n / 2;
    while (1) {
        if ((t[1].length - 1 <= k) && (r[1].length - 1 < k)) {
            struct polynomial_t locator = t[1];
            struct polynomial_t evaluator = r[1];

            symbol_t c = symb_inv(synd.field, locator.coeffs[0]);
            *locator_ptr = poly_scale(locator, c);
            *evaluator_ptr = poly_scale(evaluator, c);
        }

        poly_divmod(r[0], r[1], &q_, &r_);
        if (is_zero(r_)) {
            break;
        }
        r[0] = r[1];
        r[1] = r_;

        struct polynomial_t s_ = poly_add(s[0], poly_mult(q_, s[1]));
        struct polynomial_t t_ = poly_add(t[0], poly_mult(q_, t[1]));

        s[0] = s[1];
        s[1] = s_;

        t[0] = t[1];
        t[1] = t_;
    }

    return 1;
}

int rs_decode(const struct polynomial_t *gen, struct polynomial_t *msg) {
    assert(gen);
    assert(gen->field);
    assert(msg);
    msg->field = gen->field;

    struct polynomial_t synd = syndrome(*msg, gen->length - 1);
    if (is_zero(synd)) {
	    return 0;
    }

    struct polynomial_t locator = {0};
    struct polynomial_t evaluator = {0};
    int error = 0;

    error = key_equation_solver(synd, &locator, &evaluator);
    assert(error == 0);

    return 1;
}
