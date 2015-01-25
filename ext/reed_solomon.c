#include <stdlib.h>

#define N (8)
#define Q (1L << (N))

typedef unsigned int uint;

typedef uint symbol_t;

static const symbol_t PRIMITIVE_POLYNOMIAL = 0x11d;

struct field_t {
	symbol_t exp[Q];
	symbol_t log[Q];
	symbol_t inv[Q];
};

struct polynomial_t {
	uint length;
	symbol_t coeffs[Q];
};

int rs_init_field(struct field_t *field) {
	return 0;
}

int rs_generator(const struct field_t *field, size_t distance, struct polynomial_t *gen) {
	return 0;
}

int rs_encode(const struct field_t *field, const struct polynomial_t *gen, struct polynomial_t *msg) {
	return 0;
}

int rs_decode(const struct field_t *field, const struct polynomial_t *gen, struct polynomial_t *msg) {
	return 0;
}
