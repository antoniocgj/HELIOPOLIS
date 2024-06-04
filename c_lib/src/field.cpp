#include <sys/random.h>
#include <stdint.h>
#include <assert.h>
#include <immintrin.h>
#include <fri_c.h>

#define PRIME	((1 << 16) + 1)
#define QNR		3
#define TESTS	1000

static fp_t mod, rcp;
static unsigned int len;

void fp2_he_in(fp2_he_t a, const fp_he_t in[2]);
void fp4_he_in(fp4_he_t a, const fp_he_t in[4]);
void fp8_he_in(fp8_he_t a, const fp_he_t in[8]);
void fp16_he_in(fp16_he_t a, fp_he_t in[16]);

void fp_ini(fp_t prime) {
	mod = prime;
	len = 32 - __builtin_clzll(mod);
	rcp = 1 + ((uint64_t)1 << (2 * len)) / mod;
}


#define NUM_tmps 8
__thread fp_he_t * t;
__thread fp2_he_t t2[NUM_tmps];
__thread fp4_he_t t4[NUM_tmps];
__thread fp8_he_t t8[NUM_tmps];
static __thread bool tmp_setup = false;
void fp_setup_temps(uint64_t n){
	if(tmp_setup) return;
	t = tlwe_alloc_sample_array(NUM_tmps, n);
	for (size_t i = 0; i < NUM_tmps; i++){
		TLWE * t2_lwe = tlwe_alloc_sample_array(2, n);
		TLWE * t4_lwe = tlwe_alloc_sample_array(4, n);
		TLWE * t8_lwe = tlwe_alloc_sample_array(8, n);
		fp2_he_in(t2[i], t2_lwe);
		fp4_he_in(t4[i], t4_lwe);
		fp8_he_in(t8[i], t8_lwe);
	}
	tmp_setup = true;
}

// fp_t fp_mul(const fp_t x, const fp_t y) {
// 	uint64_t t = (uint64_t)x * y;
// 	uint32_t lo, hi;
// 	lo = (uint32_t) t;
// 	hi = t >> 32;
// 	hi = (hi << (32 - len)) | (lo >> len);
// 	lo = lo & ((1LL << len) - 1);

// 	uint32_t q = ((uint64_t)hi * rcp) >> len;
// 	uint32_t q0 = ((uint32_t)q*mod) & ((1L << len) - 1);
// 	uint32_t q1 = ((uint64_t)q*mod)>>len;
// 	q0 = (lo - q0);
// 	q1 = (hi - q1) - (lo < q0);
// 	q0 = (q0 - q1*mod) & ((1L << len) - 1);
// 	while (q0 >= mod) q0 -= mod;
// 	return q0;
// }

fp_t fp_mul(const fp_t x, const fp_t y) {
	return intel::hexl::MultiplyMod(x, y, mod);
}

// fp_t fp_mul_he(const fp_t x, const fp_he_t y) {
// 	return fp_mul(x, *y);
// }

fp_t fp_sqr(const fp_t a) {
	return fp_mul(a, a);
}

fp_t fp_add(const fp_t a, const fp_t b) {
	fp_t ret;
	ret = a + b;
	if (mod <= ret)
		ret = ret - mod;
	return ret;
}

// fp_t fp_add_he(const fp_he_t a, const fp_he_t b) {
// 	fp_t ret;
// 	ret = *a + *b;
// 	if (mod <= ret)
// 		ret = ret - mod;
// 	return ret;
// }

fp_t fp_sub(const fp_t a, const fp_t b) {
	fp_t ret;
	ret = mod - b;
	ret = ret + a;
	if (mod <= ret)
		ret = ret - mod;
	return ret;
}

fp_t fp_rnd(void) {
	uint32_t ret;
	getrandom(&ret, sizeof(ret), GRND_NONBLOCK);
	return (ret % mod);
}

fp_t fp_exp(const fp_t a, uint32_t exp) {
	fp_t ret = 1, tmp = a;
	while (exp) {
		if (exp & 1) {
			ret = fp_mul(ret, tmp);
		}
		tmp = fp_sqr(tmp);
		exp >>= 1;
	}
	return ret;
}

fp_t fp_inv(const fp_t a) {
	return fp_exp(a, mod - 2);
}

int fp_is_sqr(const fp_t a) {
	return fp_exp(a, (mod - 1)/2) == 1;
}

fp_t fp_srt(const fp_t a) {
	int f = 0, m;
	fp_t c, e, t0, t1;

	if (!fp_is_sqr(a)) {
		return 0;
	}

	/* Find a quadratic non-residue modulo p, that is a number t2
		* such that (t2 | p) = t2^((p - 1)/2)!= 1. */
	do {
		t1 = fp_rnd();
	} while (fp_is_sqr(t1));

	/* Write p - 1 as (e * 2^f), odd e. */
	e = mod - 1;
	while (e % 2 == 0) {
		e = e >> 1;
		f++;
	}

	/* Compute t2 = t2^e. */
	t1 = fp_exp(t1, e);

	/* Compute t1 = a^e, c = a^((e + 1)/2) = a^(e/2 + 1), odd e. */
	e = e >> 1;
	t0 = fp_exp(a, e);
	e = fp_mul(t0, a);
	t0 = fp_sqr(t0);
	t0 = fp_mul(t0, a);
	c = e;

	while (1) {
		if (t0 == 1) {
			break;
		}
		e = t0;
		for (m = 0; (m < f) && (t0 != 1); m++) {
			t0 = fp_sqr(t0);
		}
		t0 = e;
		for (int i = 0; i < f - m - 1; i++) {
			t1 = fp_sqr(t1);
		}
		c = fp_mul(c, t1);
		t1 = fp_sqr(t1);
		t0 = fp_mul(t0, t1);
		f = m;
	}
	return c;
}

void fp2_mul(fp2_t c, const fp2_t a, const fp2_t b) {
	fp_t t0, t1, t2, t3, t4;

	/* t2 = a_0 + a_1, t1 = b_0 + b_1. */
	t2 = fp_add(a[0], a[1]);
	t1 = fp_add(b[0], b[1]);

	/* t3 = (a_0 + a_1) * (b_0 + b_1). */
	t3 = fp_mul(t2, t1);

	/* t0 = a_0 * b_0, t4 = a_1 * b_1. */
	t0 = fp_mul(a[0], b[0]);
	t4 = fp_mul(a[1], b[1]);

	/* t2 = (a_0 * b_0) + (a_1 * b_1). */
	t2 = fp_add(t0, t4);

	/* t1 = (a_0 * b_0) + i^2 * (a_1 * b_1). */
	t1 = fp_add(t0, t4);
	for (int i = 1; i < QNR; i++) {
		t1 = fp_add(t1, t4);
	}

	/* c_0 = t1 mod p. */
	c[0] = t1;

	/* c_1 = t3 - t2. */
	c[1] = fp_sub(t3, t2);
}

void fp2_mul_he(fp2_he_t c, const fp2_t a, const fp2_he_t b) {
	fp_t t2;

	/* t2 = a_0 + a_1, t1 = b_0 + b_1. */
	t2 = fp_add(a[0], a[1]);

	// t1 = fp_add_he(b[0], b[1]);
	tlwe_add(t[1], b[0], b[1]);

	/* t3 = (a_0 + a_1) * (b_0 + b_1). */
	// t3 = fp_mul(t2, t1);
	tlwe_scale(t[3], t[1], t2);

	/* t0 = a_0 * b_0, t4 = a_1 * b_1. */
	// t0 = fp_mul_he(a[0], b[0]);
	tlwe_scale(c[0], b[0], a[0]);
	// t4 = fp_mul_he(a[1], b[1]);
	tlwe_scale(t[4], b[1], a[1]);


	/* t2 = (a_0 * b_0) + (a_1 * b_1). */
	// t2 = fp_add(t0, t4);
	tlwe_add(t[2], c[0], t[4]);

	/* t1 = (a_0 * b_0) + i^2 * (a_1 * b_1). */
	// t1 = fp_add(t0, t4);
	// for (int i = 1; i < QNR; i++) {
	// 	t1 = fp_add(t1, t4);
	// }
	tlwe_scale_addto(c[0], t[4], 3);

	/* c_0 = t1 mod p. */
	// c[0] = t1;

	/* c_1 = t3 - t2. */
	tlwe_sub(c[1], t[3], t[2]);
	// c[1] = fp_sub(t3, t2);
}

void fp2_mul_nor(fp2_t c, const fp2_t a) {
	/* (a_0 + a_1 * i) * i = (a_1 * i^2) + a_0 * i. */
	fp_t t = a[0];
	c[0] = a[1];
	for (int i = 1; i < QNR; i++) {
		c[0] = fp_add(c[0], a[1]);
	}
	c[1] = t;
}

void fp2_mul_nor_he(fp2_he_t c, const fp2_he_t a) {
	/* (a_0 + a_1 * i) * i = (a_1 * i^2) + a_0 * i. */
	tlwe_copy(t[0], a[0]);
	tlwe_scale(c[0], a[1], 3);
	tlwe_copy(c[1], t[0]);
}

void fp2_sqr(fp2_t c, const fp2_t a) {
	return fp2_mul(c, a, a);
}

void fp2_add(fp2_t c, const fp2_t a, const fp2_t b) {
	c[0] = fp_add(a[0], b[0]);
	c[1] = fp_add(a[1], b[1]);
}

void fp2_add_he(fp2_he_t c, const fp2_he_t a, const fp2_he_t b) {
	tlwe_add(c[0], a[0], b[0]);
	tlwe_add(c[1], a[1], b[1]);
}

void fp2_sub_he(fp2_he_t c, const fp2_he_t a, const fp2_he_t b) {
	tlwe_sub(c[0], a[0], b[0]);
	tlwe_sub(c[1], a[1], b[1]);
}

void fp2_sub(fp2_t c, const fp2_t a, const fp2_t b) {
	c[0] = fp_sub(a[0], b[0]);
	c[1] = fp_sub(a[1], b[1]);
}

void fp2_neg(fp2_t c, const fp2_t a) {
	c[0] = mod - a[0];
	c[1] = mod - a[1];
}

void fp2_rnd(fp2_t a) {
	a[0] = fp_rnd();
	a[1] = fp_rnd();
}

void fp2_inv(fp2_t c, const fp2_t a) {
	fp_t t0, t1;

	/* t0 = a_0^2, t1 = a_1^2. */
	t0 = fp_sqr(a[0]);
	t1 = fp_sqr(a[1]);

	/* t1 = 1/(a_0^2 + a_1^2). */
	t1 = fp_mul(t1, QNR);
	t0 = fp_sub(t0, t1);

	t1 = fp_inv(t0);

	/* c_0 = a_0/(a_0^2 + a_1^2). */
	c[0] = fp_mul(a[0], t1);
	/* c_1 = - a_1/(a_0^2 + a_1^2). */
	c[1] = mod - fp_mul(a[1], t1);
}

int fp2_cmp(const fp2_t a, const fp2_t b) {
	return (a[0] == b[0]) && (a[1] == b[1]);
}

void fp4_mul(fp4_t c, const fp4_t a, const fp4_t b) {
	fp2_t t0, t1, t2;

	/* Karatsuba algorithm. */

	/* t0 = a_0 * b_0. */
	fp2_mul(t0, a[0], b[0]);
	/* t1 = a_1 * b_1. */
	fp2_mul(t1, a[1], b[1]);
	/* t2 = b_0 + b_1. */
	fp2_add(t2, b[0], b[1]);

	/* c_1 = a_0 + a_1. */
	fp2_add(c[1], a[0], a[1]);

	/* c_1 = (a_0 + a_1) * (b_0 + b_1) */
	fp2_mul(c[1], c[1], t2);
	fp2_sub(c[1], c[1], t0);
	fp2_sub(c[1], c[1], t1);

	/* c_0 = a_0b_0 + v * a_1b_1. */
	fp2_mul_nor(t2, t1);
	fp2_add(c[0], t0, t2);
}

void fp4_mul_he(fp4_he_t c, const fp4_t a, const fp4_he_t b) {
	fp2_t t0;
	/* Karatsuba algorithm. */
	
	/* t0 = a_0 * b_0. */
	fp2_mul_he(t2[0], a[0], b[0]);
	/* t1 = a_1 * b_1. */
	fp2_mul_he(t2[1], a[1], b[1]);
	/* t2 = b_0 + b_1. */
	fp2_add_he(t2[2], b[0], b[1]);

	/* c_1 = a_0 + a_1. */
	fp2_add(t0, a[0], a[1]);

	/* c_1 = (a_0 + a_1) * (b_0 + b_1) */
	fp2_mul_he(c[1], t0, t2[2]);
	fp2_sub_he(c[1], c[1], t2[0]);
	fp2_sub_he(c[1], c[1], t2[1]);

	/* c_0 = a_0b_0 + v * a_1b_1. */
	fp2_mul_nor_he(t2[1], t2[1]);
	fp2_add_he(c[0], t2[0], t2[1]);
}

void fp4_mul_nor(fp4_t c, const fp4_t a) {
	fp2_t t0;

	/* (a_0 + a_1 * v) * v = a_0 * v + a_1 * v^2 */
	t0[0] = a[0][0];
	t0[1] = a[0][1];
	fp2_mul_nor(c[0], a[1]);
	c[1][0] = t0[0];
	c[1][1] = t0[1];
}

void fp4_mul_nor_he(fp4_he_t c, const fp4_he_t a) {
	/* (a_0 + a_1 * v) * v = a_0 * v + a_1 * v^2 */
	tlwe_copy(t2[0][0], a[0][0]);
	tlwe_copy(t2[0][1], a[0][1]);
	fp2_mul_nor_he(c[0], a[1]);
	tlwe_copy(c[1][0], t2[0][0]);
	tlwe_copy(c[1][1], t2[0][1]);
}

void fp4_sqr(fp4_t c, const fp4_t a) {
	return fp4_mul(c, a, a);
}

void fp4_add(fp4_t c, const fp4_t a, const fp4_t b) {
	fp2_add(c[0], a[0], b[0]);
	fp2_add(c[1], a[1], b[1]);
}

void fp4_add_he(fp4_he_t c, const fp4_he_t a, const fp4_he_t b) {
	fp2_add_he(c[0], a[0], b[0]);
	fp2_add_he(c[1], a[1], b[1]);
}

void fp4_sub(fp4_t c, const fp4_t a, const fp4_t b) {
	fp2_sub(c[0], a[0], b[0]);
	fp2_sub(c[1], a[1], b[1]);
}

void fp4_sub_he(fp4_he_t c, const fp4_he_t a, const fp4_he_t b) {
	fp2_sub_he(c[0], a[0], b[0]);
	fp2_sub_he(c[1], a[1], b[1]);
}

void fp4_neg(fp4_t c, const fp4_t a) {
	fp2_neg(c[0], a[0]);
	fp2_neg(c[1], a[1]);
}

void fp4_rnd(fp4_t a) {
	fp2_rnd(a[0]);
	fp2_rnd(a[1]);
}

void fp4_inv(fp4_t c, const fp4_t a) {
	fp2_t t0;
	fp2_t t1;

	fp2_sqr(t0, a[0]);
	fp2_sqr(t1, a[1]);
	fp2_mul_nor(t1, t1);
	fp2_sub(t0, t0, t1);
	fp2_inv(t0, t0);

	fp2_mul(c[0], a[0], t0);
	fp2_neg(c[1], a[1]);
	fp2_mul(c[1], c[1], t0);
}

int fp4_cmp(const fp4_t a, const fp4_t b) {
	return fp2_cmp(a[0], b[0]) && fp2_cmp(a[1], b[1]);
}

void fp8_mul(fp8_t c, const fp8_t a, const fp8_t b) {
	fp4_t t0, t1, t2;

	/* Karatsuba algorithm. */

	/* t0 = a_0 * b_0. */
	fp4_mul(t0, a[0], b[0]);
	/* t1 = a_1 * b_1. */
	fp4_mul(t1, a[1], b[1]);
	/* t2 = b_0 + b_1. */
	fp4_add(t2, b[0], b[1]);

	/* c_1 = a_0 + a_1. */
	fp4_add(c[1], a[0], a[1]);

	/* c_1 = (a_0 + a_1) * (b_0 + b_1) */
	fp4_mul(c[1], c[1], t2);
	fp4_sub(c[1], c[1], t0);
	fp4_sub(c[1], c[1], t1);

	/* c_0 = a_0b_0 + v * a_1b_1. */
	fp4_mul_nor(t2, t1);
	fp4_add(c[0], t0, t2);
}

void fp8_mul_he(fp8_he_t c, const fp8_t a, const fp8_he_t b) {
	fp4_t t0;

	/* Karatsuba algorithm. */

	/* t0 = a_0 * b_0. */
	fp4_mul_he(t4[0], a[0], b[0]);
	/* t1 = a_1 * b_1. */
	fp4_mul_he(t4[1], a[1], b[1]);
	/* t2 = b_0 + b_1. */
	fp4_add_he(t4[2], b[0], b[1]);

	/* c_1 = a_0 + a_1. */
	fp4_add(t0, a[0], a[1]);

	/* c_1 = (a_0 + a_1) * (b_0 + b_1) */
	fp4_mul_he(c[1], t0, t4[2]);
	fp4_sub_he(c[1], c[1], t4[0]);
	fp4_sub_he(c[1], c[1], t4[1]);

	/* c_0 = a_0b_0 + v * a_1b_1. */
	fp4_mul_nor_he(t4[1], t4[1]);
	fp4_add_he(c[0], t4[0], t4[1]);
}

void fp8_mul_nor(fp8_t c, const fp8_t a) {
	fp4_t t0;

	/* (a_0 + a_1 * v) * v = a_0 * v + a_1 * v^2 */
	t0[0][0] = a[0][0][0];
	t0[0][1] = a[0][0][1];
	t0[1][0] = a[0][1][0];
	t0[1][1] = a[0][1][1];
	fp4_mul_nor(c[0], a[1]);
	c[1][0][0] = t0[0][0];
	c[1][0][1] = t0[0][1];
	c[1][1][0] = t0[1][0];
	c[1][1][1] = t0[1][1];
}

void fp8_mul_nor_he(fp8_he_t c, const fp8_he_t a) {
	/* (a_0 + a_1 * v) * v = a_0 * v + a_1 * v^2 */
	tlwe_copy(t4[0][0][0], a[0][0][0]);
	tlwe_copy(t4[0][0][1], a[0][0][1]);
	tlwe_copy(t4[0][1][0], a[0][1][0]);
	tlwe_copy(t4[0][1][1], a[0][1][1]);
	fp4_mul_nor_he(c[0], a[1]);
	tlwe_copy(c[1][0][0], t4[0][0][0]);
	tlwe_copy(c[1][0][1], t4[0][0][1]);
	tlwe_copy(c[1][1][0], t4[0][1][0]);
	tlwe_copy(c[1][1][1], t4[0][1][1]);
}

void fp8_sqr(fp8_t c, const fp8_t a) {
	return fp8_mul(c, a, a);
}

void fp8_add(fp8_t c, const fp8_t a, const fp8_t b) {
	fp4_add(c[0], a[0], b[0]);
	fp4_add(c[1], a[1], b[1]);
}

void fp8_add_he(fp8_he_t c, const fp8_he_t a, const fp8_he_t b) {
	fp4_add_he(c[0], a[0], b[0]);
	fp4_add_he(c[1], a[1], b[1]);
}

void fp8_sub_he(fp8_he_t c, const fp8_he_t a, const fp8_he_t b) {
	fp4_sub_he(c[0], a[0], b[0]);
	fp4_sub_he(c[1], a[1], b[1]);
}

void fp8_sub(fp8_t c, const fp8_t a, const fp8_t b) {
	fp4_sub(c[0], a[0], b[0]);
	fp4_sub(c[1], a[1], b[1]);
}

void fp8_neg(fp8_t c, const fp8_t a) {
	fp4_neg(c[0], a[0]);
	fp4_neg(c[1], a[1]);
}

void fp8_rnd(fp8_t a) {
	fp4_rnd(a[0]);
	fp4_rnd(a[1]);
}

void fp8_inv(fp8_t c, const fp8_t a) {
	fp4_t t0;
	fp4_t t1;

	fp4_sqr(t0, a[0]);
	fp4_sqr(t1, a[1]);
	fp4_mul_nor(t1, t1);
	fp4_sub(t0, t0, t1);
	fp4_inv(t0, t0);

	fp4_mul(c[0], a[0], t0);
	fp4_neg(c[1], a[1]);
	fp4_mul(c[1], c[1], t0);
}

int fp8_cmp(const fp8_t a, const fp8_t b) {
	return fp4_cmp(a[0], b[0]) && fp4_cmp(a[1], b[1]);
}

void fp16_mul(fp16_t c, const fp16_t a, const fp16_t b) {
	fp8_t t0, t1, t2;

	/* Karatsuba algorithm. */

	/* t0 = a_0 * b_0. */
	fp8_mul(t0, a[0], b[0]);
	/* t1 = a_1 * b_1. */
	fp8_mul(t1, a[1], b[1]);
	/* t2 = b_0 + b_1. */
	fp8_add(t2, b[0], b[1]);

	/* c_1 = a_0 + a_1. */
	fp8_add(c[1], a[0], a[1]);

	/* c_1 = (a_0 + a_1) * (b_0 + b_1) */
	fp8_mul(c[1], c[1], t2);
	fp8_sub(c[1], c[1], t0);
	fp8_sub(c[1], c[1], t1);

	/* c_0 = a_0b_0 + v * a_1b_1. */
	fp8_mul_nor(t2, t1);
	fp8_add(c[0], t0, t2);
}

void fp16_mul_he(fp16_he_t c, const fp16_t a, const fp16_he_t b) {
	fp8_t t0;

	/* Karatsuba algorithm. */

	/* t0 = a_0 * b_0. */
	fp8_mul_he(t8[0], a[0], b[0]);
	/* t1 = a_1 * b_1. */
	fp8_mul_he(t8[1], a[1], b[1]);
	/* t2 = b_0 + b_1. */
	fp8_add_he(t8[2], b[0], b[1]);

	/* c_1 = a_0 + a_1. */
	fp8_add(t0, a[0], a[1]);

	/* c_1 = (a_0 + a_1) * (b_0 + b_1) */
	fp8_mul_he(c[1], t0, t8[2]);
	fp8_sub_he(c[1], c[1], t8[0]);
	fp8_sub_he(c[1], c[1], t8[1]);

	/* c_0 = a_0b_0 + v * a_1b_1. */
	fp8_mul_nor_he(t8[1], t8[1]);
	fp8_add_he(c[0], t8[0], t8[1]);
}

void fp16_get_powers(fp16_t * out, fp16_t in, uint32_t size){
	memset(out[0], 0, sizeof(out[0]));
	out[0][0][0][0][0] = 1;
	for (size_t i = 1; i < size; i++){
		fp16_mul(out[i], out[i-1], in); 
	}
}


void fp16_sqr(fp16_t c, const fp16_t a) {
	return fp16_mul(c, a, a);
}

void fp16_add(fp16_t c, const fp16_t a, const fp16_t b) {
	fp8_add(c[0], a[0], b[0]);
	fp8_add(c[1], a[1], b[1]);
}

void fp16_add_he(fp16_he_t c, const fp16_he_t a, const fp16_he_t b) {
	fp8_add_he(c[0], a[0], b[0]);
	fp8_add_he(c[1], a[1], b[1]);
}

void fp16_sub(fp16_t c, const fp16_t a, const fp16_t b) {
	fp8_sub(c[0], a[0], b[0]);
	fp8_sub(c[1], a[1], b[1]);
}

void fp16_rnd(fp16_t a) {
	fp8_rnd(a[0]);
	fp8_rnd(a[1]);
}

void fp16_inv(fp16_t c, const fp16_t a) {
	fp8_t t0;
	fp8_t t1;

	fp8_sqr(t0, a[0]);
	fp8_sqr(t1, a[1]);
	fp8_mul_nor(t1, t1);
	fp8_sub(t0, t0, t1);
	fp8_inv(t0, t0);

	fp8_mul(c[0], a[0], t0);
	fp8_neg(c[1], a[1]);
	fp8_mul(c[1], c[1], t0);
}

int fp16_cmp(const fp16_t a, const fp16_t b) {
	return fp8_cmp(a[0], b[0]) && fp8_cmp(a[1], b[1]);
}

void fp16_in(fp16_t a, fp_t in[16]) {
	a[0][0][0][0] = in[0];
	a[0][0][0][1] = in[1];
	a[0][0][1][0] = in[2];
	a[0][0][1][1] = in[3];
	a[0][1][0][0] = in[4];
	a[0][1][0][1] = in[5];
	a[0][1][1][0] = in[6];
	a[0][1][1][1] = in[7];
	a[1][0][0][0] = in[8];
	a[1][0][0][1] = in[9];
	a[1][0][1][0] = in[10];
	a[1][0][1][1] = in[11];
	a[1][1][0][0] = in[12];
	a[1][1][0][1] = in[13];
	a[1][1][1][0] = in[14];
	a[1][1][1][1] = in[15];
}

void fp16_in_u64(fp16_t a, uint64_t in[16]) {
	a[0][0][0][0] = (uint32_t) in[0];
	a[0][0][0][1] = (uint32_t) in[1];
	a[0][0][1][0] = (uint32_t) in[2];
	a[0][0][1][1] = (uint32_t) in[3];
	a[0][1][0][0] = (uint32_t) in[4];
	a[0][1][0][1] = (uint32_t) in[5];
	a[0][1][1][0] = (uint32_t) in[6];
	a[0][1][1][1] = (uint32_t) in[7];
	a[1][0][0][0] = (uint32_t) in[8];
	a[1][0][0][1] = (uint32_t) in[9];
	a[1][0][1][0] = (uint32_t) in[10];
	a[1][0][1][1] = (uint32_t) in[11];
	a[1][1][0][0] = (uint32_t) in[12];
	a[1][1][0][1] = (uint32_t) in[13];
	a[1][1][1][0] = (uint32_t) in[14];
	a[1][1][1][1] = (uint32_t) in[15];
}

void fp16_in_const64(fp16_t a, uint64_t in) {
	a[0][0][0][0] = (uint32_t) in;
	a[0][0][0][1] = 0;
	a[0][0][1][0] = 0;
	a[0][0][1][1] = 0;
	a[0][1][0][0] = 0;
	a[0][1][0][1] = 0;
	a[0][1][1][0] = 0;
	a[0][1][1][1] = 0;
	a[1][0][0][0] = 0;
	a[1][0][0][1] = 0;
	a[1][0][1][0] = 0;
	a[1][0][1][1] = 0;
	a[1][1][0][0] = 0;
	a[1][1][0][1] = 0;
	a[1][1][1][0] = 0;
	a[1][1][1][1] = 0;
}

void fp2_he_in(fp2_he_t a, const fp_he_t in[2]) {
	a[0] = in[0];
	a[1] = in[1];
}

void fp4_he_in(fp4_he_t a, const fp_he_t in[4]) {
	a[0][0] = in[0];
	a[0][1] = in[1];
	a[1][0] = in[2];
	a[1][1] = in[3];
}

void fp8_he_in(fp8_he_t a, const fp_he_t in[8]) {
	a[0][0][0] = in[0];
	a[0][0][1] = in[1];
	a[0][1][0] = in[2];
	a[0][1][1] = in[3];
	a[1][0][0] = in[4];
	a[1][0][1] = in[5];
	a[1][1][0] = in[6];
	a[1][1][1] = in[7];
}

void fp16_he_in(fp16_he_t a, fp_he_t in[16]) {
	a[0][0][0][0] = in[0];
	a[0][0][0][1] = in[1];
	a[0][0][1][0] = in[2];
	a[0][0][1][1] = in[3];
	a[0][1][0][0] = in[4];
	a[0][1][0][1] = in[5];
	a[0][1][1][0] = in[6];
	a[0][1][1][1] = in[7];
	a[1][0][0][0] = in[8];
	a[1][0][0][1] = in[9];
	a[1][0][1][0] = in[10];
	a[1][0][1][1] = in[11];
	a[1][1][0][0] = in[12];
	a[1][1][0][1] = in[13];
	a[1][1][1][0] = in[14];
	a[1][1][1][1] = in[15];
}

void fp16_out(fp_t out[16], fp16_t a) {
	out[0]  = a[0][0][0][0];
	out[1]  = a[0][0][0][1];
	out[2]  = a[0][0][1][0];
	out[3]  = a[0][0][1][1];
	out[4]  = a[0][1][0][0];
	out[5]  = a[0][1][0][1];
	out[6]  = a[0][1][1][0];
	out[7]  = a[0][1][1][1];
	out[8]  = a[1][0][0][0];
	out[9]  = a[1][0][0][1];
	out[10] = a[1][0][1][0];
	out[11] = a[1][0][1][1];
	out[12] = a[1][1][0][0];
	out[13] = a[1][1][0][1];
	out[14] = a[1][1][1][0];
	out[15] = a[1][1][1][1];
}

uint64_t mod_switch(uint64_t v, uint64_t p, uint64_t q);

void fp16_print(fp16_t a){
	printf("[%u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u]\n", a[0][0][0][0], a[0][0][0][1], a[0][0][1][0], a[0][0][1][1], a[0][1][0][0], a[0][1][0][1], a[0][1][1][0], a[0][1][1][1], a[1][0][0][0], a[1][0][0][1], a[1][0][1][0], a[1][0][1][1], a[1][1][0][0], a[1][1][0][1], a[1][1][1][0], a[1][1][1][1]);
}

#define DECRYPT(X) (mod_switch(torus2int(tlwe_phase(X, key), 52), 1ULL<<52, p)%p)
void fp16_he_print(fp16_he_t a, TLWE_Key key, uint64_t p){
	printf("[%lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu]\n", DECRYPT(a[0][0][0][0]), DECRYPT(a[0][0][0][1]), DECRYPT(a[0][0][1][0]), DECRYPT(a[0][0][1][1]), DECRYPT(a[0][1][0][0]), DECRYPT(a[0][1][0][1]), DECRYPT(a[0][1][1][0]), DECRYPT(a[0][1][1][1]), DECRYPT(a[1][0][0][0]), DECRYPT(a[1][0][0][1]), DECRYPT(a[1][0][1][0]), DECRYPT(a[1][0][1][1]), DECRYPT(a[1][1][0][0]), DECRYPT(a[1][1][0][1]), DECRYPT(a[1][1][1][0]), DECRYPT(a[1][1][1][1]));
}



// void test() {
// 	fp_t a, b, c, d, e;

// 	for (int i = 0; i < TESTS; i++) {
// 		a = fp_rnd();
// 		b = fp_sqr(a);
// 		c = fp_mul(a, a);
// 		assert(b == c);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		a = fp_rnd();
// 		b = fp_rnd();
// 		c = fp_rnd();
// 		d = fp_mul(b, c);
// 		e = fp_mul(a, d);
// 		d = fp_mul(a, b);
// 		c = fp_mul(c, d);
// 		assert(c == e);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		a = fp_rnd();
// 		b = fp_inv(a);
// 		assert(fp_mul(a, b) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		do {
// 			a = fp_rnd();
// 		} while (!fp_is_sqr(a) || a == 0);
// 		b = fp_srt(a);
// 		b = fp_sqr(b);
// 		assert(a == b);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp_he_t t;
// 		a = fp_rnd();
// 		t = &a;
// 		b = fp_sqr(a);
// 		c = fp_mul_he(a, t);
// 		assert(b == c);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp_he_t t;
// 		a = fp_rnd();
// 		t = &a;
// 		b = fp_rnd();
// 		c = fp_rnd();
// 		d = fp_mul(b, c);
// 		e = fp_mul_he(d, t);
// 		d = fp_mul_he(b, t);
// 		c = fp_mul(c, d);
// 		assert(c == e);
// 	}
// }

// void test2() {
// 	fp2_t a, b, c, d, e;

// 	for (int i = 0; i < TESTS; i++) {
// 		fp2_rnd(a);
// 		fp2_sqr(b, a);
// 		fp2_mul(c, a, a);
// 		assert(fp2_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp2_rnd(a);
// 		fp2_rnd(b);
// 		fp2_rnd(c);
// 		fp2_mul(d, b, c);
// 		fp2_mul(e, a, d);
// 		fp2_mul(d, a, b);
// 		fp2_mul(c, c, d);
// 		assert(fp2_cmp(c, e) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp2_rnd(a);
// 		fp2_inv(b, a);
// 		fp2_mul(c, a, b);
// 		fp2_mul(c, c, a);
// 		assert(fp2_cmp(a, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp2_he_t t;
// 		fp2_rnd(a);
// 		t[0] = &a[0];
// 		t[1] = &a[1];
// 		fp2_sqr(b, a);
// 		fp2_mul_he(c, a, t);
// 		assert(fp2_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp2_he_t t;
// 		fp2_rnd(a);
// 		t[0] = &a[0];
// 		t[1] = &a[1];
// 		fp2_rnd(b);
// 		fp2_rnd(c);
// 		fp2_mul(d, b, c);
// 		fp2_mul_he(e, d, t);
// 		fp2_mul_he(d, b, t);
// 		fp2_mul(c, c, d);
// 		assert(fp2_cmp(c, e) == 1);
// 	}
// }

// void test4() {
// 	fp4_t a, b, c, d, e;

// 	for (int i = 0; i < TESTS; i++) {
// 		fp4_rnd(a);
// 		fp4_sqr(b, a);
// 		fp4_mul(c, a, a);
// 		assert(fp4_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp4_rnd(a);
// 		fp4_rnd(b);
// 		fp4_rnd(c);
// 		fp4_mul(d, b, c);
// 		fp4_mul(e, a, d);
// 		fp4_mul(d, a, b);
// 		fp4_mul(c, c, d);
// 		assert(fp4_cmp(c, e) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp4_rnd(a);
// 		fp4_inv(b, a);
// 		fp4_mul(c, a, b);
// 		fp4_mul(c, c, a);
// 		assert(fp4_cmp(a, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp4_he_t t;
// 		fp4_rnd(a);
// 		t[0][0] = &a[0][0];
// 		t[0][1] = &a[0][1];
// 		t[1][0] = &a[1][0];
// 		t[1][1] = &a[1][1];
// 		fp4_sqr(b, a);
// 		fp4_mul_he(c, a, t);
// 		assert(fp4_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp4_he_t t;
// 		fp4_rnd(a);
// 		t[0][0] = &a[0][0];
// 		t[0][1] = &a[0][1];
// 		t[1][0] = &a[1][0];
// 		t[1][1] = &a[1][1];
// 		fp4_rnd(b);
// 		fp4_rnd(c);
// 		fp4_mul(d, b, c);
// 		fp4_mul_he(e, d, t);
// 		fp4_mul_he(d, b, t);
// 		fp4_mul(c, c, d);
// 		assert(fp4_cmp(c, e) == 1);
// 	}
// }

// void test8() {
// 	fp8_t a, b, c, d, e;

// 	for (int i = 0; i < TESTS; i++) {
// 		fp8_rnd(a);
// 		fp8_sqr(b, a);
// 		fp8_mul(c, a, a);
// 		assert(fp8_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp8_rnd(a);
// 		fp8_rnd(b);
// 		fp8_rnd(c);
// 		fp8_mul(d, b, c);
// 		fp8_mul(e, a, d);
// 		fp8_mul(d, a, b);
// 		fp8_mul(c, c, d);
// 		assert(fp8_cmp(c, e) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp8_rnd(a);
// 		fp8_inv(b, a);
// 		fp8_mul(c, a, b);
// 		fp8_mul(c, c, a);
// 		assert(fp8_cmp(a, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp8_he_t t;
// 		fp8_rnd(a);
// 		t[0][0][0] = &a[0][0][0];
// 		t[0][0][1] = &a[0][0][1];
// 		t[0][1][0] = &a[0][1][0];
// 		t[0][1][1] = &a[0][1][1];
// 		t[1][0][0] = &a[1][0][0];
// 		t[1][0][1] = &a[1][0][1];
// 		t[1][1][0] = &a[1][1][0];
// 		t[1][1][1] = &a[1][1][1];
// 		fp8_sqr(b, a);
// 		fp8_mul_he(c, a, t);
// 		assert(fp8_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp8_he_t t;
// 		fp8_rnd(a);
// 		t[0][0][0] = &a[0][0][0];
// 		t[0][0][1] = &a[0][0][1];
// 		t[0][1][0] = &a[0][1][0];
// 		t[0][1][1] = &a[0][1][1];
// 		t[1][0][0] = &a[1][0][0];
// 		t[1][0][1] = &a[1][0][1];
// 		t[1][1][0] = &a[1][1][0];
// 		t[1][1][1] = &a[1][1][1];
// 		fp8_rnd(b);
// 		fp8_rnd(c);
// 		fp8_mul(d, b, c);
// 		fp8_mul_he(e, d, t);
// 		fp8_mul_he(d, b, t);
// 		fp8_mul(c, c, d);
// 		assert(fp8_cmp(c, e) == 1);
// 	}
// }

// void test16() {
// 	fp16_t a, b, c, d, e;

// 	for (int i = 0; i < TESTS; i++) {
// 		fp16_rnd(a);
// 		fp16_sqr(b, a);
// 		fp16_mul(c, a, a);
// 		assert(fp16_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp16_rnd(a);
// 		fp16_rnd(b);
// 		fp16_rnd(c);
// 		fp16_mul(d, b, c);
// 		fp16_mul(e, a, d);
// 		fp16_mul(d, a, b);
// 		fp16_mul(c, c, d);
// 		assert(fp16_cmp(c, e) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp16_rnd(a);
// 		fp16_inv(b, a);
// 		fp16_mul(c, a, b);
// 		fp16_mul(c, c, a);
// 		assert(fp16_cmp(a, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp_t t[16];
// 		fp16_rnd(a);
// 		fp16_out(t, a);
// 		fp16_in(b, t);
// 		assert(fp16_cmp(a, b) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp16_he_t t;
// 		fp16_rnd(a);
// 		t[0][0][0][0] = &a[0][0][0][0];
// 		t[0][0][0][1] = &a[0][0][0][1];
// 		t[0][0][1][0] = &a[0][0][1][0];
// 		t[0][0][1][1] = &a[0][0][1][1];
// 		t[0][1][0][0] = &a[0][1][0][0];
// 		t[0][1][0][1] = &a[0][1][0][1];
// 		t[0][1][1][0] = &a[0][1][1][0];
// 		t[0][1][1][1] = &a[0][1][1][1];
// 		t[1][0][0][0] = &a[1][0][0][0];
// 		t[1][0][0][1] = &a[1][0][0][1];
// 		t[1][0][1][0] = &a[1][0][1][0];
// 		t[1][0][1][1] = &a[1][0][1][1];
// 		t[1][1][0][0] = &a[1][1][0][0];
// 		t[1][1][0][1] = &a[1][1][0][1];
// 		t[1][1][1][0] = &a[1][1][1][0];
// 		t[1][1][1][1] = &a[1][1][1][1];
// 		fp16_sqr(b, a);
// 		fp16_mul_he(c, a, t);
// 		assert(fp16_cmp(b, c) == 1);
// 	}

// 	for (int i = 0; i < TESTS; i++) {
// 		fp16_he_t t;
// 		fp16_rnd(a);
// 		t[0][0][0][0] = &a[0][0][0][0];
// 		t[0][0][0][1] = &a[0][0][0][1];
// 		t[0][0][1][0] = &a[0][0][1][0];
// 		t[0][0][1][1] = &a[0][0][1][1];
// 		t[0][1][0][0] = &a[0][1][0][0];
// 		t[0][1][0][1] = &a[0][1][0][1];
// 		t[0][1][1][0] = &a[0][1][1][0];
// 		t[0][1][1][1] = &a[0][1][1][1];
// 		t[1][0][0][0] = &a[1][0][0][0];
// 		t[1][0][0][1] = &a[1][0][0][1];
// 		t[1][0][1][0] = &a[1][0][1][0];
// 		t[1][0][1][1] = &a[1][0][1][1];
// 		t[1][1][0][0] = &a[1][1][0][0];
// 		t[1][1][0][1] = &a[1][1][0][1];
// 		t[1][1][1][0] = &a[1][1][1][0];
// 		t[1][1][1][1] = &a[1][1][1][1];
// 		fp16_rnd(b);
// 		fp16_rnd(c);
// 		fp16_mul(d, b, c);
// 		fp16_mul_he(e, d, t);
// 		fp16_mul_he(d, b, t);
// 		fp16_mul(c, c, d);
// 		assert(fp16_cmp(c, e) == 1);
// 	}

// }

// int main(void) {

// 	fp_ini(PRIME);

// 	test();
// 	test2();
// 	test4();
// 	test8();
// 	test16();
// }