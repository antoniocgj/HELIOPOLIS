#include "rns-rlwe.h"

LWE_Key lwe_alloc_key(uint64_t n){
  LWE_Key key;
  key = (LWE_Key) safe_malloc(sizeof(*key));
  key->s = (uint64_t *) safe_aligned_malloc(sizeof(*key->s)*n);
  key->n = n;
  key->q = 0;
  return key;
}

LWE lwe_alloc_sample(uint64_t n, uint64_t q){
  LWE c;
  c = (LWE) safe_malloc(sizeof(*c));
  c->a = (uint64_t *) safe_aligned_malloc(sizeof(*c->a)*n);
  c->n = n;
  c->q = q;
  return c;
}

void free_lwe_sample(LWE c){
  free(c->a);
  free(c);
}


LWE_Key lwe_new_key(uint64_t n, uint64_t q, double sec_sigma, double err_sigma){
  LWE_Key key = lwe_alloc_key(n);
  for (size_t i = 0; i < n; i++){
    key->s[i] = (q + double2int(generate_normal_random(sec_sigma)))%q;
  }
  key->n = n;
  key->q = q;
  key->sigma = err_sigma;
  return key;
}

// Generates a sparse ternary key with Hamming Weight h, balanced
LWE_Key lwe_new_sparse_ternary_key(uint64_t n, uint64_t q, uint64_t h, double err_sigma){
  LWE_Key key = lwe_alloc_key(n);
  key->n = n;
  key->q = q;
  gen_sparse_ternary_array_modq(key->s, key->n, h, q);
  key->sigma = err_sigma;
  return key;
}

void lwe_sample(LWE c, uint64_t m, LWE_Key key){
  generate_random_bytes(c->n*sizeof(*c->a), (uint8_t *) c->a);
  array_reduce_mod_N(c->a, c->a, key->n, key->q);
  array_mod_switch(c->a, c->a, next_power_of_2(key->q), key->q, key->n);
  c->b = (key->q + double2int(generate_normal_random(key->sigma)))%key->q;
  c->q = key->q;
  for (size_t i = 0; i < key->n; i++){
    const uint64_t as_i = intel::hexl::MultiplyMod(c->a[i], key->s[i], key->q); 
    c->b = intel::hexl::AddUIntMod(c->b, as_i, key->q);
  }
  c->b = intel::hexl::AddUIntMod(c->b, m, key->q);
  assert(c->b < key->q);
}

LWE lwe_new_sample(uint64_t m, LWE_Key key){
  LWE c = lwe_alloc_sample(key->n, key->q);
  lwe_sample(c, m, key);
  return c;
}

LWE lwe_new_trivial_sample(uint64_t m, uint64_t n, uint64_t q){
  LWE c = lwe_alloc_sample(n, q);
  memset(c->a, 0, sizeof(*c->a)*n);
  c->b = m;
  return c;
}

uint64_t lwe_phase(LWE c, LWE_Key key){
  uint64_t m = c->b;
  for (size_t i = 0; i < key->n; i++){
    const uint64_t as_i = intel::hexl::MultiplyMod(c->a[i], key->s[i], key->q); 
    m = intel::hexl::SubUIntMod(m, as_i, key->q);
  }
  return m;
}

void lwe_subto(LWE out, LWE in){
  assert(out->q == in->q);
  assert(out->n == in->n);
  intel::hexl::EltwiseSubMod(out->a, out->a, in->a, in->n, in->q);
  out->b = intel::hexl::SubUIntMod(out->b, in->b, in->q);
}

LWE_KS_Key lwe_new_KS_key(LWE_Key out_key, LWE_Key in_key, uint64_t t, uint64_t base_bit){
  const uint64_t base = 1ULL << base_bit, bit_size = ((uint64_t)log2(in_key->q)) + 1, q = in_key->q;
  LWE_KS_Key res;
  res = (LWE_KS_Key) safe_malloc(sizeof(*res));
  res->base_bit = base_bit;
  res->t = t;
  res->s = (LWE ***) safe_malloc(sizeof(LWE**) * in_key->n);
  for (size_t i = 0; i < in_key->n; i++){
    res->s[i] = (LWE **) safe_malloc(sizeof(LWE*) * t);
    for (size_t j = 0; j < t; j++){
      res->s[i][j] = (LWE*) safe_malloc(sizeof(LWE) * (base - 1));
      for (size_t k = 0; k < base - 1; k++){
        const uint64_t si = (q + in_key->s[i]) % q;
        const uint64_t dec_key = intel::hexl::MultiplyMod(si, (k + 1), q);
        const uint64_t dec_key2 = intel::hexl::MultiplyMod(dec_key, (1ULL << (bit_size - (j + 1) * base_bit)), q);
        if(dec_key2 < q) res->s[i][j][k] = lwe_new_sample(dec_key2, out_key);
      }
    }
  }
  return res;
}

void lwe_keyswitch(LWE out, LWE in, LWE_KS_Key ks_key){
  const uint64_t bit_size = ((uint64_t)log2(in->q)) + 1;
  const uint64_t prec_offset = ks_key->t < bit_size ? 1ULL << (bit_size - (1 + ks_key->base_bit * ks_key->t)) : 0;

  const uint64_t mask = (1ULL << ks_key->base_bit) - 1;
  assert(out->n == ks_key->s[0][0][0]->n);
  
  memset(out->a, 0, sizeof(*out->a)*out->n); // zero a
  out->b = in->b;

  for (size_t i = 0; i < in->n; i++) {
    const uint64_t ai = in->a[i] + prec_offset;
    for (size_t j = 0; j < ks_key->t; j++) {
      const uint64_t aij = (ai >> (bit_size - (j + 1) * ks_key->base_bit)) & mask;
      if (aij != 0) lwe_subto(out, ks_key->s[i][j][aij - 1]);
    }
  }
}
