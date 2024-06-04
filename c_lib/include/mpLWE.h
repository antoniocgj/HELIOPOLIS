#include <x86intrin.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <mosfhet.h>

typedef struct _LWE_params{
  uint64_t n, k, N, d;
  double rlwe_sigma;
} * LWE_params;

typedef struct _LWE{
  uint64_t ** a, * b;
  uint64_t d;
  LWE_params params;
} * LWE;

typedef struct _Polynomial{
  uint64_t ** coeffs;
  uint64_t N, d;
} * Polynomial;

typedef struct _RLWE{
  Polynomial * a, b;
  LWE_params params;
} * RLWE;

#ifdef __cplusplus
extern "C" { 
#endif 
LWE_params new_lwe_params(uint64_t n, uint64_t k, uint64_t N, uint64_t d, double rlwe_sigma);
LWE new_LWE(LWE_params params);
void LWE_zero(LWE c);
LWE * new_LWE_array(uint64_t size, LWE_params params);
void free_LWE(LWE c);
void free_LWE_array(uint64_t size, LWE * c);
void free_MP_polynomial(Polynomial p);
void LWE_propagate_carry(LWE c);
void rlwe_inner_prod(LWE * out, RLWE in, uint64_t ** alpha, uint64_t size, uint64_t pack_size);
void LWE_scale_addto(LWE out, LWE in, uint64_t scale);
Polynomial new_polynomial(uint64_t N, uint64_t d);
void polynomial_mul_by_xai(Polynomial out, Polynomial in, uint64_t a);
void polynomial_negate(Polynomial out, Polynomial in);
void polynomial_add(Polynomial out, Polynomial a, Polynomial b);
void polynomial_rnd(Polynomial poly);
void polynomial_scale(Polynomial out, Polynomial in, uint64_t scale);
void polynomial_scale_addto(Polynomial out, Polynomial in, uint64_t scale);
void polynomial_zero(Polynomial poly);
void polynomial_propagate_carry(Polynomial p);
void polynomial_mul_addto_sparse_polynomial(Polynomial out, Polynomial a, uint64_t * b, uint64_t size);
RLWE new_rlwe(LWE_params params);
RLWE * new_rlwe_array(uint64_t size, LWE_params params);
void rlwe_get_pos_key(uint64_t * out, uint64_t size, uint64_t * in, uint64_t N, uint64_t p);
void rlwe_encrypt(RLWE c, uint64_t * m, uint64_t p, uint64_t * key, uint64_t key_size);
void rlwe_scale(RLWE out, RLWE in, uint64_t scale);
void rlwe_scale_addto(RLWE out, RLWE in, uint64_t scale);
void rlwe_zero(RLWE c);
void rlwe_add(RLWE out, RLWE a, RLWE b);
void free_rlwe_array(uint64_t size, RLWE * c);
void polynomial_drop_digits(Polynomial p, uint64_t num_digits);
void rlwe_drop_digits(RLWE c, uint64_t num_digits);
void rlwe_propagate_carry(RLWE c);
uint64_t mod_switch(uint64_t v, uint64_t p, uint64_t q);

#ifdef __cplusplus
}
#endif 