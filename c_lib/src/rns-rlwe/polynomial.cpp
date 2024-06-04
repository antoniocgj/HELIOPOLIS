#include "rns-rlwe.h"

RNS_Polynomial polynomial_new_RNS_polynomial(uint64_t N, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_Polynomial res;
  res = (RNS_Polynomial) safe_malloc(sizeof(*res));
  res->coeffs = (uint64_t **) safe_malloc(sizeof(uint64_t*) * l);
  for (size_t i = 0; i < l; i++){
    res->coeffs[i] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t) * N);
  }
  res->N = N;
  res->l = l;
  res->ntt = ntt;
  return res;
}

RNS_Polynomial * polynomial_new_RNS_polynomial_array(uint64_t size, uint64_t N, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_Polynomial * res;
  res = (RNS_Polynomial *) safe_malloc(sizeof(RNS_Polynomial)*size);
  for (size_t i = 0; i < size; i++){
    res[i] = polynomial_new_RNS_polynomial(N, l, ntt);
  }
  return res;
}

void polynomial_copy_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in){
  for (size_t i = 0; i < in->l; i++){
    memcpy(out->coeffs[i], in->coeffs[i], sizeof(uint64_t)*out->N);
  }
}

void polynomial_RNS_zero(RNS_Polynomial p){
  for (size_t i = 0; i < p->l; i++){
    memset(p->coeffs[i], 0, sizeof(uint64_t)*p->N);
  }
}

void free_RNS_polynomial(void * p){
  RNS_Polynomial pp = (RNS_Polynomial) p;
  for (size_t i = 0; i < pp->l; i++){
    free(pp->coeffs[i]);
  }
  free(pp->coeffs);
  free(pp);
}

void free_RNS_polynomial_array(uint64_t size, RNS_Polynomial * p){
  for (size_t i = 0; i < size; i++){
    free_RNS_polynomial(p[i]);
  }
  free(p);
}

RNS_Polynomial * polynomial_new_array_of_RNS_polynomials(uint64_t N, uint64_t l, uint64_t size, intel::hexl::NTT ** ntt){
  RNS_Polynomial * res = (RNS_Polynomial *) safe_malloc(sizeof(RNS_Polynomial)*size);
  for (size_t i = 0; i < size; i++) res[i] = polynomial_new_RNS_polynomial(N, l, ntt);
  return res;
}

// out = RNS(in)
// Assumes ||in||_inf < min(p)^2
void polynomial_to_RNS(RNS_Polynomial out, IntPolynomial in){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->N; j++){
      out->coeffs[i][j] = (p + in->coeffs[j])%p;
    }
    out->ntt[i]->ComputeForward(out->coeffs[i], out->coeffs[i], 1, 1);
  }
}

void array_to_RNS(RNS_Polynomial out, uint64_t * in){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->N; j++){
      out->coeffs[i][j] = (p + in[j])%p;
    }
    out->ntt[i]->ComputeForward(out->coeffs[i], out->coeffs[i], 1, 1);
  }
}

void array_conj_to_RNS(RNS_Polynomial out, uint64_t * in){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt[i]->GetModulus();
    out->coeffs[i][0] = (p + in[0])%p;
    for (size_t j = 1; j < out->N; j++){
      out->coeffs[i][j] = (p - in[out->N - j])%p;
    }
    out->ntt[i]->ComputeForward(out->coeffs[i], out->coeffs[i], 1, 1);
  }
}


// todo: fix uniform distribution
void polynomial_gen_random_RNSc_polynomial(RNSc_Polynomial out){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt[i]->GetModulus();
    generate_random_bytes(sizeof(uint64_t)*out->N, (uint8_t *) out->coeffs[i]);
    array_mod_switch_from_2k(out->coeffs[i], out->coeffs[i], p, p, out->N);
  }
}

/* out = in1*in2 mod (X^N + 1) */
void polynomial_mul_RNS_polynomial2(RNSc_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseMultMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q, 1);
    out->ntt[i]->ComputeInverse(out->coeffs[i], out->coeffs[i], 1, 1);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

/* out = in1*in2 mod (X^N + 1) */
void polynomial_mul_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseMultMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q, 1);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

/* out = in1 - in2 mod (X^Q - 1) */
void polynomial_sub_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseSubMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

/* out = in1 - in2 mod (X^Q - 1) */
void polynomial_sub_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseSubMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}


void polynomial_add_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseAddMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

void polynomial_scale_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, uint64_t scale){
  const uint64_t l = in1->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseFMAMod(out->coeffs[i], in1->coeffs[i], scale, nullptr, in1->N, q, 1);
  }
}

void polynomial_RNSc_negate(RNSc_Polynomial out, RNSc_Polynomial in){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->N; j++){
      out->coeffs[i][j] = intel::hexl::SubUIntMod(0, in->coeffs[i][j], q);
    }
  }
}


void polynomial_RNSc_to_RNS(RNS_Polynomial out, RNSc_Polynomial in){
  for (size_t i = 0; i < out->l; i++){
    out->ntt[i]->ComputeForward(out->coeffs[i], in->coeffs[i], 1, 1);
  }
}

void polynomial_RNS_to_RNSc(RNSc_Polynomial out, RNS_Polynomial in){
  for (size_t i = 0; i < out->l; i++){
    out->ntt[i]->ComputeInverse(out->coeffs[i], in->coeffs[i], 1, 1);
  }
}

// out = in + e(sigma)
void polynomial_RNSc_add_noise(RNSc_Polynomial out, RNSc_Polynomial in, double sigma){
  for (size_t j = 0; j < out->N; j++){
    const double real_noise = generate_normal_random(sigma);
    for (size_t i = 0; i < out->l; i++){
      const uint64_t p = out->ntt[i]->GetModulus();
      const uint64_t noise = ((uint64_t) (p + real_noise))%p;
      out->coeffs[i][j] = intel::hexl::AddUIntMod(in->coeffs[i][j], noise, p);
    }
  }
}

/* Assumes q_i/q_j < 2 for all i, j*/
void polynomial_base_extend_RNSc(RNSc_Polynomial out, uint64_t * in, uint64_t p){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    if(q >= p){ // copy
      memcpy(out->coeffs[i], in, sizeof(uint64_t)*out->N);
    }else{ // reduce
      intel::hexl::EltwiseReduceMod(out->coeffs[i], in, out->N, q, 2, 1);
    }
  }
}

// Removes the last prime of the RNS polynomial
void polynomial_base_reduce_RNSc_wo_free(RNSc_Polynomial out){
  const uint64_t l = out->l, p = out->ntt[l-1]->GetModulus();
  const uint64_t N = out->N;
  for (size_t i = 0; i < l - 1; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    const uint64_t inv_p = intel::hexl::InverseMod(p, q);
    if(q >= p){ // copy
      for (size_t j = 0; j < N; j++){
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }else{ // reduce
      for (size_t j = 0; j < N; j++){
        const uint64_t in_j = intel::hexl::ReduceMod<2>(out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], in_j, q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }
  } 
  out->l -= 1;
}


// Removes the last prime of the RNS polynomial
void polynomial_base_reduce_RNSc(RNSc_Polynomial out){
  polynomial_base_reduce_RNSc_wo_free(out);
  free(out->coeffs[out->l]);
}

// Removes the last prime of the RNS polynomial
void polynomial_base_reduce_RNSc_and_scale(RNSc_Polynomial out, uint64_t p){
  polynomial_base_reduce_RNSc_wo_free(out);
  polynomial_scale_RNSc_polynomial(out, out, (out->ntt[out->l]->GetModulus())%p);
  free(out->coeffs[out->l]);
}



void polynomial_RNSc_permute(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t gen){
  const uint64_t N = out->N;
  assert(gen < N);
  assert(gen > 0);
  uint64_t idx = 0;
  polynomial_RNS_zero((RNS_Polynomial) out);
  for (size_t i = 0; i < N; i++){
    for (size_t j = 0; j < out->l; j++){
      out->coeffs[j][idx] = in->coeffs[j][i];
    }
    idx = intel::hexl::AddUIntMod(idx, gen, N);
  }
}

void polynomial_int_permute_mod_Q(IntPolynomial out, IntPolynomial in, uint64_t gen){
  const uint64_t N = in->N;
  uint64_t idx = 0;
  for (size_t i = 0; i < N; i++){
    out->coeffs[idx] = in->coeffs[i];
    idx = intel::hexl::AddUIntMod(idx, gen, N);
  }
}

// todo: check
void polynomial_RNSc_mul_by_xai(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t a){
  assert(in != out);
  const uint64_t N = out->N;
  assert(a < N);
  polynomial_RNS_zero((RNS_Polynomial) out);
  for (size_t j = 0; j < out->l; j++){
    for (size_t i = 0; i < N-a; i++){
      out->coeffs[j][i+a] = in->coeffs[j][i];
    }
    for (size_t i = N-a; i < N; i++){
      out->coeffs[j][i-N+a] = -in->coeffs[j][i];
    }
  }
}

void polynomial_int_decompose_i(IntPolynomial out, IntPolynomial in, uint64_t Bg_bit, uint64_t l, uint64_t q, uint64_t bit_size, uint64_t i){
  const uint64_t N = in->N;
  const uint64_t h_mask = (1UL << Bg_bit) - 1;
  const uint64_t h_bit = bit_size - (i + 1) * Bg_bit;

  uint64_t offset = 1ULL << (bit_size - l * Bg_bit - 1);

  for (size_t c = 0; c < N; c++){
    const uint64_t coeff_off = in->coeffs[c] + offset;
    out->coeffs[c] = (coeff_off>>h_bit) & h_mask;
  }
}

IntPolynomial polynomial_new_int_polynomial(uint64_t N){
  IntPolynomial res;
  res = (IntPolynomial) safe_malloc(sizeof(*res));
  res->coeffs = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t) * N);
  res->N = N;
  return res;
}

// void free_polynomial(void * p){
//   free(((IntPolynomial) p)->coeffs);
//   free(p);
// }