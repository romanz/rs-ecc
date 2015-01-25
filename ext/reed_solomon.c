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

void poly_divmod(struct polynomial_t f, struct polynomial_t g, struct polynomial_t *r_, struct polynomial_t *q_) {
	assert(r_);
	assert(q_);
}

struct field_t *field_alloc() {
	return malloc(sizeof(struct field_t));
}

void field_free(struct field_t *f) {
	free(f);
}

struct field_t *polynomial_alloc() {
	return malloc(sizeof(struct polynomial_t));
}

void polynomial_free(struct polynomial_t *p) {
	free(p);
}

void polynomial_print(const struct polynomial_t *p, const char *fmt) {
	assert(p);
	printf("[");
	printf(fmt, p->coeffs[0]);
	for (int i = 1; i < p->length; ++i) {
		printf(", ");
		printf(fmt, p->coeffs[i]);
	}
	printf("]");
}

int rs_init_field(struct field_t *field, symbol_t primitive_poly) {
	assert(field);
	symbol_t alpha = 1;
	for (int i = 0; i < Q-1; ++i) {
		field->exp[i] = alpha;
		// printf("%5d: %5d\n", i, alpha);

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

int rs_generator(const struct field_t *field, uint distance, struct polynomial_t *gen) {
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

int rs_encode(const struct polynomial_t *gen, struct polynomial_t *msg) {
	return 0;
}

int rs_decode(const struct polynomial_t *gen, struct polynomial_t *msg) {
	return 0;
}
