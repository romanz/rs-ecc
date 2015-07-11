#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define N (16)
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

#define dprintf(X)

void poly_print(const struct polynomial_t *p, const char *fmt) {
	assert(p);
    dprintf(("["));
    dprintf((fmt, p->coeffs[0]));
    for (int i = 1; i < p->length; ++i) {
        dprintf((", "));
        dprintf((fmt, p->coeffs[i]));
    }
    dprintf(("]"));
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

void poly_trim(struct polynomial_t *p) {
    assert(p);
    while (p->length > 1) {
        if (p->coeffs[p->length - 1] != 0) {
            break;
        }
        p->length -= 1;
    }
}

void poly_mult(const struct polynomial_t *f, const struct polynomial_t *g,
    		   struct polynomial_t *result) {
    assert(f);
    assert(f->length);
    assert(g);
    assert(g->length);
    assert(result);
    assert(f->field);
    assert(f->field == g->field);
    const struct field_t *field = f->field;

    result->length = f->length + g->length - 1;
    assert(result->length <= Q);
    for (int i = 0; i < result->length; ++i) {
    	result->coeffs[i] = 0;
    }
    result->field = field;

    for (int i = 0; i < f->length; ++i) {
        for (int j = 0; j < g->length; ++j) {
            int k = i + j;
            result->coeffs[k] ^= symb_mult(field, f->coeffs[i], g->coeffs[j]);
        }
    }
}

void poly_add(const struct polynomial_t *f, const struct polynomial_t *g,
    		  struct polynomial_t *result) {
    assert(f);
    assert(g);
    assert(result);
    assert(f->field);
    assert(f->field == g->field);
    const struct field_t *field = f->field;

    uint length = (f->length > g->length) ? f->length : g->length;

    for (int i = 0; i < length; ++i) {
    	symbol_t f_i = (i < f->length) ? f->coeffs[i] : 0;
    	symbol_t g_i = (i < g->length) ? g->coeffs[i] : 0;
        result->coeffs[i] = f_i ^ g_i;
    }
    result->field = field;
    result->length = length;
}

void poly_scale(const struct polynomial_t *f, symbol_t c,
				struct polynomial_t *res) {
	assert(f);
	assert(res);
    assert(f->field);

    for (int i = 0; i < f->length; ++i) {
        res->coeffs[i] = symb_mult(f->field, f->coeffs[i], c);
    }
    res->field = f->field;
    res->length = f->length;
}

symbol_t poly_eval(const struct polynomial_t *f, symbol_t x) {
	assert(f);
    symbol_t res = 0;
    symbol_t z = 1;
    for (int i = 0; i < f->length; ++i) {
        res ^= symb_mult(f->field, f->coeffs[i], z);
        z = symb_mult(f->field, z, x);
    }
    return res;
}

void syndrome(const struct polynomial_t *f, uint size,
			  struct polynomial_t *res) {
    assert(res);
    assert(f);
    assert(f->field);
    for (int i = 0; i < size; ++i) {
        res->coeffs[i] = poly_eval(f, f->field->exp[i]);
    }
    res->field = f->field;
    res->length = size;
    poly_trim(res);
}

struct field_t *field_alloc() {
	struct field_t *f = calloc(1, sizeof(struct field_t));
	assert(f);
	dprintf(("field: %p\n", f));
	return f;
}

void field_free(struct field_t *f) {
	free(f);
}

struct polynomial_t *poly_alloc(const struct field_t *f) {
	struct polynomial_t *p = calloc(1, sizeof(struct polynomial_t));
	dprintf(("poly: %p\n", p));
	assert(p);
    if (p) {
    	p->length = 1;
        p->field = f;
    }
    return p;
}

void poly_free(struct polynomial_t *p) {
	free(p);
}

void poly_divmod(const struct polynomial_t *f, const struct polynomial_t *g,
                 struct polynomial_t *q, struct polynomial_t *r) {
	assert(f);
	assert(g);
    assert(f->field);
    assert(f->field == g->field);
    assert(r);
    const struct field_t *field = f->field;

    symbol_t factor = symb_inv(field, g->coeffs[g->length - 1]);

    struct polynomial_t *g_norm = poly_alloc(field);
    poly_scale(f, factor, r);
    poly_scale(g, factor, g_norm);
    struct polynomial_t *g_scaled = poly_alloc(field);

    for (int i = f->length - 1; i >= g->length - 1; --i) {
    	poly_print(r, "%d");
    	dprintf(("\n"));
        int offset = i - g->length + 1;
        assert(offset >= 0);

        symbol_t c = r->coeffs[i];
        poly_scale(g_norm, c, g_scaled);
        for (int j = 0; j < g->length; ++j) {
            r->coeffs[offset + j] ^= g_scaled->coeffs[j];
        }
        assert(r->coeffs[i] == 0);
        r->length -= 1;
        if (q) {
        	q->coeffs[offset] = c;
        	q->length += 1;
    	}
    }
    poly_trim(r);  // remainder must be trimmed.
    if (q) {
        q->field = field;
    }

    poly_free(g_norm);
    poly_free(g_scaled);
}

///////////////////////////////////////////////////////////////////////////////


int rs_field(struct field_t *field, symbol_t primitive_poly) {
	assert(field);
	symbol_t alpha = 1;

	for (int i = 0; i < Q-1; ++i) {
		assert(i < Q);
		field->exp[i] = alpha;

		assert(field->log[alpha] == 0);
		assert(alpha < Q);
		field->log[alpha] = i;

		alpha = alpha << 1;
		if (alpha & (1 << N)) {
			alpha = alpha ^ primitive_poly;
		}
	}
	assert(alpha == 1);
	field->exp[Q-1] = 1;

	for (int x = 1; x < Q; ++x) {
		uint e = Q - 1 - field->log[x];
		symbol_t y = field->exp[e];
		assert(symb_mult(field, x, y) == 1);
		assert(x < Q);
		field->inv[x] = y;
	}
	return 0;
}

int is_zero(const struct polynomial_t *p) {
	assert(p);
    for (int i = 0; i < p->length; ++i) {
        if (p->coeffs[i]) {
            return 0;
        }
    }
    return 1;
}

void poly_copy(const struct polynomial_t *src, struct polynomial_t *dst) {
	assert(src);
	assert(dst);

	dst->length = src->length;
	dst->field = src->field;
	for (int i = 0; i < src->length; ++i) {
		dst->coeffs[i] = src->coeffs[i];
	}
}

struct polynomial_t *poly_clone(const struct polynomial_t *src) {
	assert(src);
	assert(src->field);

	struct polynomial_t *dst = poly_alloc(src->field);
	assert(dst);
	dst->length = src->length;
	for (int i = 0; i < src->length; ++i) {
		dst->coeffs[i] = src->coeffs[i];
	}
	return dst;
}

int rs_generator(const struct field_t *field, uint distance,
                 struct polynomial_t *result) {
	dprintf(("field: %p\n", field));
	dprintf(("result: %p\n", result));
	assert(field);
	assert(result);

	result->length = 1;
	result->coeffs[0] = 1;
	result->field = field;

	struct polynomial_t *left = poly_alloc(field);
	struct polynomial_t *right = poly_alloc(field);

	for (int i = 0; i < distance; ++i) {
		poly_copy(result, left);
		right->coeffs[0] = field->exp[i];;
		right->coeffs[1] = 1;
		right->length = 2;
		poly_mult(left, right, result);
	}

	poly_print(result, "%d");
	poly_free(left);
	poly_free(right);

	return 0;
}

int rs_encode(const struct polynomial_t *gen, struct polynomial_t *msg) {
    assert(gen);
    assert(gen->field);
    assert(msg);
    msg->field = gen->field;

	uint distance = gen->length - 1;
    assert(msg->length + distance < Q);

    // Multiply by msg x^d
    for (int i = msg->length - 1; i >= 0; --i) {
        msg->coeffs[i + distance] = msg->coeffs[i];
        msg->coeffs[i] = 0;
    }
    msg->length += distance;

    // Compute remainder by dividing msg by generator.
    struct polynomial_t *r = poly_alloc(msg->field);
    poly_divmod(msg, gen, NULL, r);
    for (int i = 0; i < distance; ++i) {
        msg->coeffs[i] = r->coeffs[i];
    }
    poly_free(r);
    return 0;
}

void poly_deriv(const struct polynomial_t *p, struct polynomial_t *result) {
	assert(p);
	assert(result);
    for (int i = 1; i < p->length; ++i) {
        result->coeffs[i-1] = p->coeffs[i] * (i % 2);
    }
    result->length = p->length - 1;
    result->field = p->field;
}

int key_equation_solver(const struct polynomial_t *synd,
                        struct polynomial_t *locator,
                        struct polynomial_t *evaluator) {
	assert(synd);
	assert(synd->field);
	assert(locator);
	assert(evaluator);
	const struct field_t *field = synd->field;

    uint n = synd->length;
    struct polynomial_t *a = poly_alloc(field);
    a->length = n + 1;
    a->coeffs[n] = 1;

    const struct polynomial_t *b = synd;
    assert(a->length > b->length);

    struct polynomial_t *r[2] = {a, poly_clone(b)};

    struct polynomial_t *s[2] = {poly_alloc(field), poly_alloc(field)};
    s[0]->coeffs[0] = 1;

    struct polynomial_t *t[2] = {poly_alloc(field), poly_alloc(field)};
    s[1]->coeffs[0] = 1;

    struct polynomial_t *q_ = poly_alloc(field);
    struct polynomial_t *r_ = poly_alloc(field);

    struct polynomial_t *tmp = poly_alloc(field);

    uint k = n / 2;
    int success = 0;
    while (1) {
        if ((t[1]->length - 1 <= k) && (r[1]->length - 1 < k)) {
            symbol_t c = symb_inv(field, locator->coeffs[0]);
            poly_scale(t[1], c, locator);
            poly_scale(r[1], c, evaluator);
            success = 1;
            break;
        }
        poly_trim(r[1]);
        poly_divmod(r[0], r[1], q_, r_);

        if (is_zero(r_)) {
            break;
        }
        poly_copy(r[1], r[0]);
        poly_copy(r_, r[1]);

        poly_mult(q_, s[1], tmp);
        poly_add(s[0], tmp, tmp);
        poly_trim(tmp);  // tmp = s[0] - q * s[1]
        poly_copy(s[1], s[0]);
        poly_copy(tmp, s[1]);

        poly_mult(q_, t[1], tmp);
        poly_add(t[0], tmp, tmp);
        poly_trim(tmp);  // tmp = t[0] - q * t[1]
        poly_copy(t[1], t[0]);
        poly_copy(tmp, t[1]);
    }

    for (int i = 0; i < 2; ++i) {
    	poly_free(r[i]);
    	poly_free(s[i]);
    	poly_free(t[i]);
    }
    poly_free(q_);
    poly_free(r_);
    poly_free(tmp);
    return success;
}

int correct(const struct polynomial_t *gen, struct polynomial_t *msg) {
    assert(gen);
    assert(gen->field);
    assert(msg);
    const struct field_t *field = gen->field;
    msg->field = field;

    uint distance = gen->length - 1;
    struct polynomial_t *synd = poly_alloc(field);
    syndrome(msg, distance, synd);
    if (is_zero(synd)) {
        return 0;
    }
    dprintf(("syndrome = "));
    poly_print(synd, "%d");

    struct polynomial_t *locator = poly_alloc(field);
    struct polynomial_t *evaluator = poly_alloc(field);

    int error = key_equation_solver(synd, locator, evaluator);
    if (error) {
        return -1;
    }

    dprintf(("\nlocator = "));
    poly_print(locator, "%d");
    dprintf(("\nevaluator = "));
    poly_print(evaluator, "%d");
    dprintf(("\n"));

    struct polynomial_t *locator_deriv = poly_alloc(field);
    poly_deriv(locator, locator_deriv);

    uint corrections = 0;
    for (int j = 0; j < msg->length; ++j) {
        symbol_t alpha = field->exp[j];
        symbol_t inv_alpha = symb_inv(field, alpha);
        dprintf(("%d %d => %d\n", j, alpha, inv_alpha));
        if (poly_eval(locator, inv_alpha) != 0) {
            continue;
        }
        symbol_t num = symb_mult(field, poly_eval(evaluator, inv_alpha), alpha);
        symbol_t den = poly_eval(locator_deriv, inv_alpha);

        msg->coeffs[j] ^= symb_mult(field, num, symb_inv(field, den));
        corrections += 1;
    }

    syndrome(msg, distance, synd);
    if (is_zero(synd)) {
        return corrections;
    } else {
        return -1;  // ECC failure
    }
}

int rs_decode(const struct polynomial_t *gen, struct polynomial_t *msg) {
    int result = correct(gen, msg);
    uint distance = gen->length - 1;
    if (result >= 0) {  // ECC success
        for (int i = distance; i < msg->length; ++i) {
            msg->coeffs[i - distance] = msg->coeffs[i];
        }
        msg->length -= distance;
    }
    return result;
}
