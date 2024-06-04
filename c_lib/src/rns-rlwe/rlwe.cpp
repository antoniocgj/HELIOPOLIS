#include "rns-rlwe.h"
// RLWE RNS functions

RNS_RLWE_Key rlwe_alloc_RNS_key(uint64_t N, uint64_t l, uint64_t log_mod, intel::hexl::NTT ** ntt, double sigma){
  RNS_RLWE_Key res;
  res = (RNS_RLWE_Key) safe_malloc(sizeof(*res));
  res->log_mod = log_mod;
  res->sigma = sigma;
  res->N = N;
  res->l = l;
  res->s = polynomial_new_int_polynomial(N);
  res->s_RNS = polynomial_new_RNS_polynomial(N, l, ntt);
  return res;
}

void free_rlwe_RNS_key(RNS_RLWE_Key key){
  free_polynomial(key->s);
  free_RNS_polynomial(key->s_RNS);
  free(key);
}

RNS_RLWE_Key rlwe_new_RNS_key(uint64_t N, uint64_t l, uint64_t log_mod, intel::hexl::NTT ** ntt, double sigma){
  RNS_RLWE_Key res = rlwe_alloc_RNS_key(N, l, log_mod, ntt, sigma);
  const uint64_t mod_mask = (1ULL << log_mod) - 1;
  generate_random_bytes(N*sizeof(uint64_t), (uint8_t *) res->s->coeffs);
  for (size_t i = 0; i < N; i++) res->s->coeffs[i] &= mod_mask; 
  polynomial_to_RNS(res->s_RNS, res->s);
  return res;
}

// void rlwe_RNSc_extract_lwe_key(LWE_Key out, RNS_RLWE_Key in){
//   const uint64_t Q = in->Q;
//   assert(Q == out->n);
//   out->q = in->s_RNS->ntt[0]->GetModulus();
//   for (size_t i = 0; i < Q; i++){
//     out->s[i] = in->s->coeffs[i];
//   }
// }


// negacyclic extract
void rlwe_RNSc_extract_lwe(uint64_t * out, RNSc_RLWE in, uint64_t idx){
  const uint64_t N = in->a->N;
  assert(in->a->l == 1);
  for (size_t i = 0; i <= idx; i++){
    out[i] = in->a->coeffs[0][idx - i];
  }
  for (size_t i = idx + 1; i < N; i++){
    out[i] = -in->a->coeffs[0][N + idx - i];
  }
  out[N] = in->b->coeffs[0][idx];
}

RNS_RLWE_Key rlwe_new_RNS_gaussian_key(uint64_t N, uint64_t l, double key_sigma, intel::hexl::NTT ** ntt, double sigma){
  RNS_RLWE_Key res = rlwe_alloc_RNS_key(N, l, ceil(key_sigma), ntt, sigma);
  // clwe_sample_gaussian(res->s->coeffs, Q, key_sigma);
  polynomial_to_RNS(res->s_RNS, res->s);
  return res;
}

RNS_RLWE_Key rlwe_get_RNS_key_from_array(uint64_t N, uint64_t l, uint64_t * array, intel::hexl::NTT ** ntt, double sigma){
  RNS_RLWE_Key res = rlwe_alloc_RNS_key(N, l, 32, ntt, sigma);
  for (size_t i = 0; i < N; i++) res->s->coeffs[i] = array[i];
  polynomial_to_RNS(res->s_RNS, res->s);
  return res;
}

RNS_RLWE rlwe_alloc_RNS_sample(uint64_t N, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_RLWE res;
  res = (RNS_RLWE) safe_malloc(sizeof(*res));
  res->a = polynomial_new_RNS_polynomial(N, l, ntt);
  res->b = polynomial_new_RNS_polynomial(N, l, ntt);
  return res;
}

RNSc_RLWE rlwe_alloc_RNSc_sample(uint64_t N, uint64_t l, intel::hexl::NTT ** ntt){
  RNSc_RLWE res;
  res = (RNSc_RLWE) safe_malloc(sizeof(*res));
  res->a = (RNSc_Polynomial) polynomial_new_RNS_polynomial(N, l, ntt);
  res->b = (RNSc_Polynomial) polynomial_new_RNS_polynomial(N, l, ntt);
  return res;
}


RNS_RLWE * rlwe_alloc_RNS_sample_array(uint64_t size, uint64_t N, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_RLWE * res;
  res = (RNS_RLWE *) safe_malloc(size*sizeof(*res));
  for (size_t i = 0; i < size; i++){
    res[i] = rlwe_alloc_RNS_sample(N, l, ntt);
  }
  return res;
}

RNS_RLWE * rlwe_alloc_RNS_sample_array2(uint64_t size, RNS_RLWE c){
  RNS_RLWE * res;
  res = (RNS_RLWE *) safe_malloc(size*sizeof(*res));
  for (size_t i = 0; i < size; i++){
    res[i] = rlwe_alloc_RNS_sample(c->b->N, c->b->l, c->b->ntt);
  }
  return res;
}

void free_RNS_rlwe_array(uint64_t size, RNS_RLWE * v){
  for (size_t i = 0; i < size; i++){
    free_rlwe_RNS_sample(v[i]);
  }
  free(v);
}


void free_RNS_rlwe_sample(RNS_RLWE c){
  free_RNS_polynomial(c->a);
  free_RNS_polynomial(c->b);
  free(c);
}

void rlwe_copy_RNS_sample(RNS_RLWE out, RNS_RLWE in){
  polynomial_copy_RNS_polynomial(out->a, in->a);
  polynomial_copy_RNS_polynomial(out->b, in->b);
}

void rlwe_copy_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in){
  rlwe_copy_RNS_sample((RNS_RLWE) out, (RNS_RLWE) in);
}

void free_rlwe_RNS_sample(void * p){
  const RNS_RLWE pp = (RNS_RLWE) p;
  free_RNS_polynomial(pp->a);
  free_RNS_polynomial(pp->b);
  free(pp);
}

void rlwe_RNSc_sample_of_zero(RNSc_RLWE out, RNS_RLWE_Key key){
  const uint64_t N = key->N, l = key->l;
  RNS_Polynomial tmp = polynomial_new_RNS_polynomial(N, l, key->s_RNS->ntt);
  polynomial_gen_random_RNSc_polynomial(out->a);
  polynomial_RNSc_to_RNS(tmp, out->a);
  polynomial_mul_RNS_polynomial2(out->b, key->s_RNS, tmp);
  polynomial_RNSc_add_noise(out->b, out->b, key->sigma);
  free_RNS_polynomial(tmp);
}


void rlwe_scale_RNSc_rlwe(RNSc_RLWE c, uint64_t scale){
  for (size_t j = 0; j < c->a->l; j++){
    const uint64_t mod = c->a->ntt[j]->GetModulus();
    for (size_t k = 0; k < c->a->N; k++){
      c->a->coeffs[j][k] = intel::hexl::MultiplyMod(c->a->coeffs[j][k], scale, mod);
      c->b->coeffs[j][k] = intel::hexl::MultiplyMod(c->b->coeffs[j][k], scale, mod);
    } 
  }
}

// out += in*scale
void rlwe_scale_RNS_rlwe_addto(RNS_RLWE out, RNS_RLWE in, uint64_t scale){
  for (size_t j = 0; j < in->a->l; j++){
    const uint64_t mod = in->a->ntt[j]->GetModulus();
    intel::hexl::EltwiseFMAMod(out->a->coeffs[j], in->a->coeffs[j], scale, out->a->coeffs[j], in->a->N, mod, 1);
    intel::hexl::EltwiseFMAMod(out->b->coeffs[j], in->b->coeffs[j], scale, out->b->coeffs[j], in->b->N, mod, 1);
  }
}


void rlwe_RNSc_mod_switch(RNSc_RLWE c, uint64_t q){
  const uint64_t p = c->a->ntt[0]->GetModulus();
  array_mod_switch(c->a->coeffs[0], c->a->coeffs[0], p, q, c->a->N);
  array_mod_switch(c->b->coeffs[0], c->b->coeffs[0], p, q, c->a->N);
}

// Return Q_0 * X^m
RNSc_RLWE rlwe_new_RNSc_sample(RNS_RLWE_Key key, uint64_t m){
  RNSc_RLWE res = rlwe_new_RNSc_sample_of_zero(key);
  const uint64_t p = res->b->ntt[0]->GetModulus();
  res->b->coeffs[0][m] = intel::hexl::AddUIntMod(res->b->coeffs[0][m], 1, p);
  return res;
}

RNS_RLWE rlwe_new_RNS_sample2(RNS_RLWE_Key key, uint64_t * m, uint64_t p){
  RNSc_RLWE res = rlwe_new_RNSc_sample_of_zero(key);
  const uint64_t q = res->b->ntt[0]->GetModulus();
  uint64_t qi_prod = 1;
  for (size_t i = 1; i < res->a->l; i++) 
    qi_prod = intel::hexl::MultiplyMod(qi_prod, res->a->ntt[i]->GetModulus(), q);
  for (size_t i = 0; i < key->N; i++){
    const uint64_t val = intel::hexl::MultiplyMod(mod_switch(m[i], p, q), qi_prod, q);
    res->b->coeffs[0][i] = intel::hexl::AddUIntMod(res->b->coeffs[0][i], val, q);
  }
  rlwe_RNSc_to_RNS((RNS_RLWE) res, res);
  return (RNS_RLWE) res;
}


RNS_RLWE rlwe_new_RNS_sample(RNS_RLWE_Key key, uint64_t * m, uint64_t p){
  RNSc_RLWE res = rlwe_new_RNSc_sample_of_zero(key);
  uint64_t qi_prod = 1;
  for (size_t i = 0; i < res->a->l; i++) 
    qi_prod = intel::hexl::MultiplyMod(qi_prod, (res->a->ntt[i]->GetModulus())%p, p);
  for (size_t j = 0; j < key->l; j++){
    const uint64_t q = res->b->ntt[j]->GetModulus();
    const uint64_t minus_p_inv_q = q - intel::hexl::InverseMod(p, q);
    for (size_t i = 0; i < key->N; i++){
      const uint64_t val = intel::hexl::MultiplyMod(m[i], qi_prod, p);
      const uint64_t val2 = intel::hexl::MultiplyMod(val, minus_p_inv_q, q);
      res->b->coeffs[j][i] = intel::hexl::AddUIntMod(res->b->coeffs[j][i], val2, q);
    }
  }
  rlwe_RNSc_to_RNS((RNS_RLWE) res, res);
  return (RNS_RLWE) res;
}

void rlwe_RNS_phase(RNS_Polynomial out, RNS_RLWE in, RNS_RLWE_Key key){
  polynomial_mul_RNS_polynomial(out, in->a, key->s_RNS);
  polynomial_sub_RNS_polynomial(out, in->b, out);
}

void rlwe_RNS_mul_by_poly(RNS_RLWE out, RNS_RLWE in, RNS_Polynomial poly){
  polynomial_mul_RNS_polynomial(out->a, in->a, poly);
  polynomial_mul_RNS_polynomial(out->b, in->b, poly);
}

RNSc_RLWE rlwe_new_RNSc_sample_of_zero(RNS_RLWE_Key key){
  RNSc_RLWE res = (RNSc_RLWE) rlwe_alloc_RNS_sample(key->N, key->l, key->s_RNS->ntt);
  rlwe_RNSc_sample_of_zero(res, key);
  return res;
}

void rlwe_RNS_sample_of_zero(RNS_RLWE out, RNS_RLWE_Key key){
  rlwe_RNSc_sample_of_zero((RNSc_RLWE) out, key);
  rlwe_RNSc_to_RNS(out, (RNSc_RLWE) out);
}

RNS_RLWE rlwe_new_RNS_sample_of_zero(RNS_RLWE_Key key){
  RNS_RLWE res = rlwe_alloc_RNS_sample(key->N, key->l, key->s_RNS->ntt);
  rlwe_RNS_sample_of_zero(res, key);
  return res;
}

RNS_RLWE rlwe_new_RNS_trivial_sample_of_zero(uint64_t N, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_RLWE res = rlwe_alloc_RNS_sample(N, l, ntt);
  for (size_t i = 0; i < res->a->l; i++){
    memset(res->a->coeffs[i], 0, sizeof(uint64_t)*res->a->N);
    memset(res->b->coeffs[i], 0, sizeof(uint64_t)*res->b->N);
  }
  return res;
}

void rlwe_RNS_trivial_sample_of_zero(RNS_RLWE out){
  for (size_t i = 0; i < out->a->l; i++){
    memset(out->a->coeffs[i], 0, sizeof(uint64_t)*out->a->N);
    memset(out->b->coeffs[i], 0, sizeof(uint64_t)*out->b->N);
  }
}

void rlwe_shrink_RNSc_sample(RNSc_RLWE c){
  polynomial_base_reduce_RNSc(c->a);
  polynomial_base_reduce_RNSc(c->b);
}

void rlwe_automorphism_RNSc(RNSc_RLWE out, RNSc_RLWE in, uint64_t gen, RNS_RLWE_KS_Key ksk){
  RNSc_RLWE tmp = (RNSc_RLWE) rlwe_alloc_RNS_sample(out->a->N, out->a->l, out->a->ntt);
  polynomial_RNSc_permute(tmp->a, in->a, gen);  
  polynomial_RNSc_permute(tmp->b, in->b, gen);
  rlwe_RNSc_keyswitch(out, tmp, ksk);
  free_rlwe_RNS_sample(tmp);
}

void rlwe_addto_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in){
  polynomial_add_RNSc_polynomial(out->a, out->a, in->a);
  polynomial_add_RNSc_polynomial(out->b, out->b, in->b);
}


void rlwe_RNSc_to_RNS(RNS_RLWE out, RNSc_RLWE in){
  polynomial_RNSc_to_RNS(out->a, in->a);
  polynomial_RNSc_to_RNS(out->b, in->b);
}

void rlwe_RNS_to_RNSc(RNSc_RLWE out, RNS_RLWE in){
  polynomial_RNS_to_RNSc(out->a, in->a);
  polynomial_RNS_to_RNSc(out->b, in->b);
}

RNS_RLWE_KS_Key rlwe_new_RNS_ks_key(RNS_RLWE_Key out_key, RNS_RLWE_Key in_key){
  const uint64_t l = in_key->l;
  RNS_RLWE_KS_Key res;
  res = (RNS_RLWE_KS_Key) safe_malloc(sizeof(*res));
  res->s = (RNS_RLWE*) safe_malloc(sizeof(RNS_RLWE)*l);
  res->l = l;
  for (size_t i = 0; i < l; i++){
    res->s[i] = (RNS_RLWE) rlwe_new_RNSc_sample_of_zero(out_key);
    const uint64_t p = res->s[i]->b->ntt[i]->GetModulus();
    for (size_t j = 0; j < in_key->N; j++){
      res->s[i]->b->coeffs[i][j] = intel::hexl::AddUIntMod(res->s[i]->b->coeffs[i][j], (p+in_key->s->coeffs[j])%p, p);
    }
    rlwe_RNSc_to_RNS(res->s[i], (RNSc_RLWE) res->s[i]);
  }
  return res;
}

void naive_full_mul(uint64_t * out, uint64_t * in1, uint64_t * in2, size_t N){
  memset(out, 0, sizeof(uint64_t)*N*2);
  for (size_t i = 0; i < N; i++){
    for (size_t j = 0; j < N; j++){
      out[i + j] += in1[i]*in2[j];
    }
  }
}

void free_rlwe_RNS_ks_key(RNS_RLWE_KS_Key key){
  for (size_t i = 0; i < key->l; i++){
    free_rlwe_RNS_sample(key->s[i]);
  }
  free(key->s);
  free(key);
}

RNS_RLWE_KS_Key rlwe_new_RNS_automorphism_key(RNS_RLWE_Key key, uint64_t gen){
  RNS_RLWE_Key perm_key = rlwe_alloc_RNS_key(key->N, key->l, key->log_mod, key->s_RNS->ntt, key->sigma);
  polynomial_int_permute_mod_Q(perm_key->s, key->s, gen);
  // It's not necessary to permute s_RNS
  RNS_RLWE_KS_Key res = rlwe_new_RNS_ks_key(key, perm_key);
  free_rlwe_RNS_key(perm_key);
  return res;
}


void rlwe_RNSc_mul_by_xai(RNSc_RLWE out, RNSc_RLWE in, uint64_t a){
  polynomial_RNSc_mul_by_xai(out->a, in->a, a);
  polynomial_RNSc_mul_by_xai(out->b, in->b, a);
}

void rlwe_RNSc_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key ksk){
  assert(in != out);
  rlwe_RNS_trivial_sample_of_zero((RNS_RLWE)out);
  RNSc_Polynomial tmp = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->N, in->a->l, in->a->ntt);
  RNSc_Polynomial tmp2 = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->N, in->a->l, in->a->ntt);
  const uint64_t l = in->a->l;
  for (size_t i = 0; i < l; i++){
    const uint64_t p = in->a->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in->a->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], ksk->s[i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], ksk->s[i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }
  polynomial_RNS_to_RNSc(out->a, (RNS_Polynomial) out->a);
  polynomial_RNS_to_RNSc(out->b, (RNS_Polynomial) out->b);
  polynomial_add_RNSc_polynomial(out->b, out->b, in->b);
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
}

void rlwe_RNSc_priv_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key kska, RNS_RLWE_KS_Key kskb){
  rlwe_RNS_trivial_sample_of_zero((RNS_RLWE)out);
  RNSc_Polynomial tmp = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->N, in->a->l, in->a->ntt);
  RNSc_Polynomial tmp2 = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->N, in->a->l, in->a->ntt);
  const uint64_t l = in->a->l;
  // a
  for (size_t i = 0; i < l; i++){
    const uint64_t p = in->a->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in->a->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kska->s[i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kska->s[i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }
  // b
  for (size_t i = 0; i < l; i++){
    const uint64_t p = in->a->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in->b->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kskb->s[i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kskb->s[i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }
  // reduce
  polynomial_RNS_to_RNSc(out->a, (RNS_Polynomial) out->a);
  polynomial_RNS_to_RNSc(out->b, (RNS_Polynomial) out->b);
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
}