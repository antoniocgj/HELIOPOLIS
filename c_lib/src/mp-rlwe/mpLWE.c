#include <mpLWE.h>

uint64_t double2int(double x){
  return ((uint64_t) ((int64_t) x));
}

LWE_params new_lwe_params(uint64_t n, uint64_t k, uint64_t N, uint64_t d, double rlwe_sigma){
  LWE_params res = (LWE_params) safe_malloc(sizeof(*res));
  res->n = n;
  res->N = N;
  res->k = k;
  res->d = d;
  res->rlwe_sigma = rlwe_sigma;
  return res;
}

Polynomial new_polynomial(uint64_t N, uint64_t d){
  Polynomial res = (Polynomial) safe_malloc(sizeof(*res));
  res->coeffs = (uint64_t **) safe_malloc(sizeof(*res->coeffs)*d);
  res->d = d;
  res->N = N;
  for (size_t i = 0; i < d; i++){
    res->coeffs[i] = (uint64_t *) safe_aligned_malloc(sizeof(**(res->coeffs))*N);
    memset(res->coeffs[i], 0, sizeof(*(res->coeffs))*N);
  }
  return res;
}

LWE new_LWE(LWE_params params){
  const uint64_t n = params->n, d = params->d;
  LWE res = (LWE) safe_malloc(sizeof(*res));
  res->a = (uint64_t **) safe_malloc(sizeof(*res->a)*d);
  res->d = d;
  for (size_t i = 0; i < d; i++){
    res->a[i] = (uint64_t *) safe_aligned_malloc(sizeof(**(res->a))*n);
    memset(res->a[i], 0, sizeof(*(res->a))*n);
  }
  const uint64_t d8 = ((d>>3) + 1) << 3;
  res->b = (uint64_t *) safe_aligned_malloc(sizeof(res->b)*d8);
  memset(res->b, 0, sizeof(res->b)*d8);
  res->params = params;
  return res;
}

LWE * new_LWE_array(uint64_t size, LWE_params params){
  LWE * res = (LWE *) safe_malloc(sizeof(LWE)*size);
  for (size_t i = 0; i < size; i++){
    res[i] = new_LWE(params);
  }
  return res;
}

void free_LWE_array(uint64_t size, LWE * c){
  for (size_t i = 0; i < size; i++){
    free_LWE(c[i]);
  }
  free(c);
}

void free_LWE(LWE c){
  for (size_t i = 0; i < c->d; i++){
    free(c->a[i]);
  }
  free(c->a);
  free(c->b);
  free(c);
}


// void lwe_extract(LWE out, RLWE in, int idx){
//   const int N = in->b->N, d = in->b->d;
//   for (size_t i = 0; i < d; i++){
//     for (size_t j = 0; j <= idx; j++){
//       out->a[i][j] = in->a[0]->coeffs[i][idx - j];
//     }
//     for (size_t j = idx + 1; j < N; j++){
//       out->a[i][j] = -in->a[0]->coeffs[i][N + idx - j];
//     }
//     out->b[i] = in->b->coeffs[i][idx];
//   }
// }

void rlwe_inner_prod(LWE * out, RLWE in, uint64_t ** alpha, uint64_t size, uint64_t pack_size){
  LWE tmp = new_LWE(out[0]->params);
  Polynomial a_neg = new_polynomial(in->a[0]->N, in->a[0]->d);
  polynomial_negate(a_neg, in->a[0]);
  polynomial_propagate_carry(a_neg);
  const int N = in->b->N, d = in->b->d;
  assert(pack_size <= N);
  assert(out[0]->params->n == in->b->N);
  for (size_t idx = 0; idx < pack_size; idx++){
    // mul by x^idx
    for (size_t i = 0; i < d; i++){
      for (size_t j = 0; j <= idx; j++){
        tmp->a[i][j] = in->a[0]->coeffs[i][idx - j];
      }
      for (size_t j = idx + 1; j < N; j++){
        tmp->a[i][j] = a_neg->coeffs[i][N + idx - j];
      }
      tmp->b[i] = in->b->coeffs[i][idx];
    }
    /* scale and add */
    for (size_t i = 0; i < size; i++){
      if(idx == 0) LWE_zero(out[i]);
      LWE_scale_addto(out[i], tmp, alpha[idx][i]);
      if((i&((1<<10)-1)) == 0 && i > 0){
        LWE_propagate_carry(out[i]);
      }
    }
  }
  free_MP_polynomial(a_neg);
  free_LWE(tmp);
}


void free_MP_polynomial(Polynomial p){
  for (size_t i = 0; i < p->d; i++){
    free(p->coeffs[i]);
  }
  free(p->coeffs);
  free(p);
}

// out = in*X^a mod (X^N + 1)
void polynomial_mul_by_xai(Polynomial out, Polynomial in, uint64_t a){
  const int N = out->N;
  assert(out->d == in->d);
  assert(out->N == in->N);
  a &= ((N<<1) - 1); // a % 2N
  if (!a) return;
  for (size_t j = 0; j < in->d; j++){
    if (a < N) {
      for (int i = 0; i < a; i++) out->coeffs[j][i] = -in->coeffs[j][i - a + N];
      for (int i = a; i < N; i++) out->coeffs[j][i] = in->coeffs[j][i - a];
    }else{
      for (int i = 0; i < a - N; i++) out->coeffs[j][i] = in->coeffs[j][i - a + 2*N];
      for (int i = a - N; i < N; i++) out->coeffs[j][i] = -in->coeffs[j][i - a + N];
    }
  }
}


void polynomial_negate(Polynomial out, Polynomial in){
  __m512i neg_mask = _mm512_set1_epi64(0x000FFFFFFFFFFFFFULL);
  __m512i neg_mask2 = _mm512_set1_epi64(1);
  for (size_t j = 0; j < in->d; j++){
    __m512i * inv = (__m512i *) in->coeffs[j];
    __m512i * outv = (__m512i *) out->coeffs[j];
    for (size_t i = 0; i < in->N/8; i++){
      outv[i] = _mm512_xor_epi64(inv[i], neg_mask);
      if(j == in->d-1) outv[i] = _mm512_add_epi64(outv[i], neg_mask2);
    }
  }
}

void polynomial_add(Polynomial out, Polynomial a, Polynomial b){
  for (size_t j = 0; j < a->d; j++){
    __m512i * av = (__m512i *) a->coeffs[j];
    __m512i * bv = (__m512i *) b->coeffs[j];
    __m512i * outv = (__m512i *) out->coeffs[j];
    for (size_t i = 0; i < a->N/8; i++){
      outv[i] = _mm512_add_epi64(av[i], bv[i]);
    }
  }
}

void polynomial_drop_digits(Polynomial p, uint64_t num_digits){
  for (size_t i = 1; i <= num_digits; i++){
    free(p->coeffs[p->d - i]);
  }
  p->d -= num_digits;
}

void polynomial_rnd(Polynomial poly){
  __m512i and_mask = _mm512_set1_epi64(0x000FFFFFFFFFFFFFULL);
  for (size_t j = 0; j < poly->d; j++){
    generate_random_bytes(sizeof(uint64_t)*poly->N, (uint8_t *) poly->coeffs[j]);
    __m512i * polyv = (__m512i *) poly->coeffs[j];
    for (size_t i = 0; i < poly->N/8; i++){
      polyv[i] = _mm512_and_epi64(polyv[i], and_mask);
    }
  }
}

// polynomials should not have carry bits
void polynomial_scale(Polynomial out, Polynomial in, uint64_t scale){
  __m512i zero = _mm512_setzero_si512();
  __m512i scalev = _mm512_set1_epi64(scale);
  for (size_t j = 0; j < in->d; j++){
    __m512i * inv = (__m512i *) in->coeffs[j];
    __m512i * outv_prev = (__m512i *) out->coeffs[j - 1];
    __m512i * outv = (__m512i *) out->coeffs[j];
    for (size_t i = 0; i < in->N/8; i++){
      if(j > 0) outv_prev[i] = _mm512_madd52hi_epu64 (outv_prev[i], inv[i], scalev);
      outv[i] = _mm512_madd52lo_epu64 (zero, inv[i], scalev);
    }
  }
}

void LWE_scale_addto(LWE out, LWE in, uint64_t scale){
  const __m512i scalev = _mm512_set1_epi64(scale);
  const __m512i idx_b_mask = _mm512_set_epi64(7, 7, 6, 5, 4, 3, 2, 1);

  // a d=0
  __m512i * inv = (__m512i *) in->a[0];
  __m512i * outv = (__m512i *) out->a[0];
  for (size_t i = 0; i < in->params->n/8; i+=4){
    outv[i] = _mm512_madd52lo_epu64 (outv[i], inv[i], scalev);
    outv[i+1] = _mm512_madd52lo_epu64 (outv[i+1], inv[i+1], scalev);
    outv[i+2] = _mm512_madd52lo_epu64 (outv[i+2], inv[i+2], scalev);
    outv[i+3] = _mm512_madd52lo_epu64 (outv[i+3], inv[i+3], scalev);
  }
  for (size_t j = 1; j < in->d; j++){
    // a
    __m512i * inv = (__m512i *) in->a[j];
    __m512i * outv_prev = (__m512i *) out->a[j - 1];
    __m512i * outv = (__m512i *) out->a[j];
    for (size_t i = 0; i < in->params->n/8; i+=4){
      outv_prev[i] = _mm512_madd52hi_epu64 (outv_prev[i], inv[i], scalev);
      outv_prev[i+1] = _mm512_madd52hi_epu64 (outv_prev[i+1], inv[i+1], scalev);
      outv_prev[i+2] = _mm512_madd52hi_epu64 (outv_prev[i+2], inv[i+2], scalev);
      outv_prev[i+3] = _mm512_madd52hi_epu64 (outv_prev[i+3], inv[i+3], scalev);
      outv[i] = _mm512_madd52lo_epu64 (outv[i], inv[i], scalev);
      outv[i+1] = _mm512_madd52lo_epu64 (outv[i+1], inv[i+1], scalev);
      outv[i+2] = _mm512_madd52lo_epu64 (outv[i+2], inv[i+2], scalev);
      outv[i+3] = _mm512_madd52lo_epu64 (outv[i+3], inv[i+3], scalev);
    }
  }
  // b
  // todo: loop for digits > 8
  __m512i * invb = (__m512i *) in->b;
  __m512i * outvb = (__m512i *) out->b;
  __m512i tmp = _mm512_permutexvar_epi64 (idx_b_mask, invb[0]);
  outvb[0] = _mm512_madd52lo_epu64 (outvb[0], invb[0], scalev);
  outvb[0] = _mm512_mask_madd52hi_epu64 (outvb[0], 0x7F, tmp, scalev);
}

void LWE_scale_addto2(LWE out, LWE in, uint64_t scale){
  const __m512i scalev = _mm512_set1_epi64(scale);
  const __m512i idx_b_mask = _mm512_set_epi64(7, 7, 6, 5, 4, 3, 2, 1);
  for (size_t j = 0; j < in->d; j++){
    // a
    __m512i * inv = (__m512i *) in->a[j];
    __m512i * outv_prev = (__m512i *) out->a[j - 1];
    __m512i * outv = (__m512i *) out->a[j];
    for (size_t i = 0; i < in->params->n/8; i++){
      if(j > 0) outv_prev[i] = _mm512_madd52hi_epu64 (outv_prev[i], inv[i], scalev);
      outv[i] = _mm512_madd52lo_epu64 (outv[i], inv[i], scalev);
    }
  }
  // b
  // todo: loop for digits > 8
  __m512i * invb = (__m512i *) in->b;
  __m512i * outvb = (__m512i *) out->b;
  __m512i tmp = _mm512_permutexvar_epi64 (idx_b_mask, invb[0]);
  outvb[0] = _mm512_madd52lo_epu64 (outvb[0], invb[0], scalev);
  outvb[0] = _mm512_mask_madd52hi_epu64 (outvb[0], 0x7F, tmp, scalev);
}

void polynomial_sp_scale_mp(Polynomial out, Polynomial in, __m512i * scale){
  __m512i zero = _mm512_setzero_si512();
  for (size_t j = 0; j < in->d; j++){
    __m512i * inv_sp = (__m512i *) in->coeffs[in->d - 1];
    __m512i * outv_prev = (__m512i *) out->coeffs[j - 1];
    __m512i * outv = (__m512i *) out->coeffs[j];
    for (size_t i = 0; i < in->N/8; i++){
      if(j > 0) outv_prev[i] = _mm512_madd52hi_epu64(outv_prev[i], scale[j], inv_sp[i]);
      outv[i] = _mm512_madd52lo_epu64(zero, scale[j], inv_sp[i]);
    }
  }
}

void polynomial_scale_addto(Polynomial out, Polynomial in, uint64_t scale){
  __m512i scalev = _mm512_set1_epi64(scale);
  for (size_t j = 0; j < in->d; j++){
    __m512i * inv = (__m512i *) in->coeffs[j];
    __m512i * outv_prev = (__m512i *) out->coeffs[j - 1];
    __m512i * outv = (__m512i *) out->coeffs[j];
    for (size_t i = 0; i < in->N/8; i++){
      if(j > 0) outv_prev[i] = _mm512_madd52hi_epu64 (outv_prev[i], inv[i], scalev);
      outv[i] = _mm512_madd52lo_epu64 (outv[i], inv[i], scalev);
    }
  }
}

void polynomial_zero(Polynomial poly){
  for (size_t j = 0; j < poly->d; j++){
    memset(poly->coeffs[j], 0, sizeof(uint64_t)*poly->N);
  }
}

void LWE_zero(LWE c){
  for (size_t j = 0; j < c->d; j++){
    memset(c->a[j], 0, sizeof(uint64_t)*c->params->n);
  }
  memset(c->b, 0, sizeof(uint64_t)*c->d);
}


void polynomial_propagate_carry(Polynomial p){
  __m512i mod_mask = _mm512_set1_epi64(0x000FFFFFFFFFFFFFULL);
  for (int64_t j = p->d - 1; j > 0; j--){
    __m512i * pv = (__m512i *) p->coeffs[j];
    __m512i * pv_next = (__m512i *) p->coeffs[j - 1];
    for (size_t i = 0; i < p->N/8; i++){
      // propagate carry 
      // this can be cheaper in enchange for some bits of precision using _mm512_madd52hi_epu64
      pv_next[i] = _mm512_add_epi64(pv_next[i], _mm512_srli_epi64(pv[i], 52));
      pv[i] = _mm512_and_epi64(pv[i], mod_mask); // remove MSB
    }
  }
  __m512i * pv = (__m512i *) p->coeffs[0];
  for (size_t i = 0; i < p->N/8; i++){
    pv[i] = _mm512_and_epi64(pv[i], mod_mask); // remove MSB
  }
}

void LWE_propagate_carry(LWE c){
  __m512i mod_mask = _mm512_set1_epi64(0x000FFFFFFFFFFFFFULL);
  for (int64_t j = c->d - 1; j > 0; j--){
    // a
    __m512i * pv = (__m512i *) c->a[j];
    __m512i * pv_next = (__m512i *) c->a[j - 1];
    for (size_t i = 0; i < c->params->n/8; i++){
      // propagate carry 
      // this can be cheaper in enchange for some bits of precision using _mm512_madd52hi_epu64
      pv_next[i] = _mm512_add_epi64(pv_next[i], _mm512_srli_epi64(pv[i], 52));
      pv[i] = _mm512_and_epi64(pv[i], mod_mask); // remove MSB
    }
    // b
    c->b[j-1] += (c->b[j] >> 52);
    c->b[j] &= 0x000FFFFFFFFFFFFFULL;
  }
  __m512i * pv = (__m512i *) c->a[0];
  for (size_t i = 0; i < c->params->n/8; i++){
    pv[i] = _mm512_and_epi64(pv[i], mod_mask); // remove MSB
  }
  c->b[0] &= 0x000FFFFFFFFFFFFFULL;
}

// out = a*b mod (X^N + 1)
void polynomial_mul_addto_sparse_polynomial(Polynomial out, Polynomial a, uint64_t * b, uint64_t size){
  const uint64_t N = a->N, N_mask = N - 1;
  Polynomial a_neg = new_polynomial(N, a->d);
  polynomial_negate(a_neg, a);
  polynomial_propagate_carry(a_neg);
  for (size_t i = 0; i < size; i++){
    for (size_t j = 0; j < N; j++){
      const uint64_t pos = (b[i] + j)&N_mask, sign = (((b[i]&N_mask) + j)>=N)^(b[i]>>63);
      for (size_t k = 0; k < a->d; k++){
        if(sign) out->coeffs[k][j] += a_neg->coeffs[k][pos];
        else out->coeffs[k][j] += a->coeffs[k][pos];
      }      
    }
  }
  polynomial_propagate_carry(out);
  free_MP_polynomial(a_neg);
}

RLWE new_rlwe(LWE_params params){
  RLWE res = (RLWE) safe_malloc(sizeof(*res));
  res->params = params;
  res->a = (Polynomial *) safe_malloc(params->k*sizeof(Polynomial *));
  for (size_t i = 0; i < params->k; i++){
    res->a[i] = new_polynomial(params->N, params->d);
  }
  res->b = new_polynomial(params->N, params->d);
  return res;
}

void free_rlwe(RLWE c){
  for (size_t i = 0; i < c->params->k; i++){
    free_MP_polynomial(c->a[i]);
  }  
  free_MP_polynomial(c->b);
  free(c->a);
  free(c);
}

RLWE * new_rlwe_array(uint64_t size, LWE_params params){
  RLWE * res = (RLWE *) safe_malloc(sizeof(RLWE)*size);
  for (size_t i = 0; i < size; i++){
    res[i] = new_rlwe(params);
  }
  return res;
}

void free_rlwe_array(uint64_t size, RLWE * c){
  for (size_t i = 0; i < size; i++){
    free_rlwe(c[i]);
  }
  free(c);
}

// returns the key in positional notation with sign bit
void rlwe_get_pos_key(uint64_t * out, uint64_t size, uint64_t * in, uint64_t N, uint64_t p){
  const uint64_t sign_mask = 1ULL << 63, p_half = p == 0? sign_mask : p>>1;
  uint64_t idx = 0;
  if(in[0] != 0){
    out[idx] = 0;
    if(in[0] > p_half) out[idx] |= sign_mask;
    idx++;
  }
  for (size_t i = 1; i < N; i++){
    if(in[i] != 0){
      out[idx] = N - i;
      if(in[i] < p_half) out[idx] |= sign_mask;
      assert(idx < size);
      idx++;
    }
  }
}


uint64_t array32_bit_slice52(uint64_t * array, uint64_t start){
  uint64_t res = 0;
  const uint64_t word_idx = start/32, word_shift = 52 - 32 + (start%32), mod_mask = (1ULL<<52)-1;
  res |= array[word_idx] << word_shift;
  if(word_shift < 32){
    res |= array[word_idx + 1] >> (32 - word_shift);
  }else{
    res |= array[word_idx + 1] << (word_shift - 32);
    res |= array[word_idx + 2] >> (64 - word_shift);
  }
  res &= mod_mask;
  // printf("delta: %lx\n", res);
  return res;
}

__m512i * _delta = NULL;
uint64_t _p = 0, _d = 0;
void setup_mod_switch_delta(uint64_t d, uint64_t p) {
  if(p == _p && _d >= d){ 
    return;
  }else{
    if(_delta) free(_delta);
    _delta = (__m512i *) safe_aligned_malloc(sizeof(__m512i)*d);
  }
  assert(p < (1ULL<<32));
  uint64_t delta32[2*d], rem = 1;
  memset(delta32, 0, sizeof(uint64_t)*2*d);
  for (size_t i = 0; i < 2*d; i++){
    delta32[i] = (rem<<32)/p;
    rem = (rem<<32)%p;
    // printf("delta32: %lx\n", delta32[i]);
  }
  for (size_t i = 0; i < d; i++){
    _delta[i] = _mm512_set1_epi64(array32_bit_slice52(delta32, 52*i));
  }
  _p = p;
  _d = d;
}

void rlwe_encrypt(RLWE c, uint64_t * m, uint64_t p, uint64_t * key, uint64_t key_size){
  assert(p < (1ULL<<52));
  polynomial_zero(c->b);
  setup_mod_switch_delta(c->b->d, p);
  for (size_t i = 0; i < c->b->N; i++){
    c->b->coeffs[c->b->d - 1][i] = m[i];
  }
  polynomial_sp_scale_mp(c->b, c->b, _delta);
  for (size_t i = 0; i < c->params->k; i++){
    polynomial_rnd(c->a[i]);
    polynomial_mul_addto_sparse_polynomial(c->b, c->a[i], &key[key_size*i], key_size);
  }
  for (size_t i = 0; i < c->b->N; i++){
    c->b->coeffs[c->b->d - 1][i] += double2int(generate_normal_random(c->params->rlwe_sigma));
  }
  polynomial_propagate_carry(c->b);
}

void rlwe_scale(RLWE out, RLWE in, uint64_t scale){
  for (size_t i = 0; i < in->params->k; i++){
    polynomial_scale(out->a[i], in->a[i], scale);
  }  
  polynomial_scale(out->b, in->b, scale);
}

void rlwe_scale_addto(RLWE out, RLWE in, uint64_t scale){
  for (size_t i = 0; i < in->params->k; i++){
    polynomial_scale_addto(out->a[i], in->a[i], scale);
  }  
  polynomial_scale_addto(out->b, in->b, scale);
}

void rlwe_zero(RLWE c){
  for (size_t i = 0; i < c->params->k; i++){
    polynomial_zero(c->a[i]);
  }  
  polynomial_zero(c->b);
}

void rlwe_propagate_carry(RLWE c){
  for (size_t i = 0; i < c->params->k; i++){
    polynomial_propagate_carry(c->a[i]);
  }  
  polynomial_propagate_carry(c->b);
}


// out = a + b
void rlwe_add(RLWE out, RLWE a, RLWE b){
  for (size_t i = 0; i < b->params->k; i++){
    polynomial_add(out->a[i], a->a[i], b->a[i]);
  }  
  polynomial_add(out->b, a->b, b->b);
}

void rlwe_drop_digits(RLWE c, uint64_t num_digits){
  for (size_t i = 0; i < c->params->k; i++){
    polynomial_drop_digits(c->a[i], num_digits);
  }  
  polynomial_drop_digits(c->b, num_digits);
}