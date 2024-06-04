#include <fri_c.h>

void _ntt(RNS_RLWE * out, RNS_RLWE * in, uint64_t size, uint64_t * ws, uint64_t w_power, uint64_t depth, uint64_t maskMod, uint64_t m, uint64_t ntt_max_depth){
  // printf("_NTT: %lu, %lu, %lu, %lu\n", depth, ntt_max_depth, m, size);
  if(depth == ntt_max_depth || (size%m)){
    assert(out != in);
    for (size_t i = 0; i < size; i++){
      rlwe_RNS_trivial_sample_of_zero(out[i]);
      for (size_t j = 0; j < size; j++){
        rlwe_scale_RNS_rlwe_addto(out[i], in[j], ws[(w_power*i*j)&maskMod]);
      }
    }
    return;
  }
  RNS_RLWE g[size/m];
  RNS_RLWE * h[m];
  for (size_t i = 0; i < m; i++){
    for (size_t j = 0; j < size/m; j++){
      g[j] = in[m*j + i];
    }
    h[i] = rlwe_alloc_RNS_sample_array2(size/m, out[0]);
    _ntt(h[i], g, size/m, ws, m*w_power, depth+1, maskMod, m, ntt_max_depth);
  }
  for (size_t k1 = 0; k1 < size/m; k1++){
    for (size_t k2 = 0; k2 < m; k2++){
      rlwe_RNS_trivial_sample_of_zero(out[k1 + (size/m)*k2]);
      for (size_t i = 0; i < m; i++){
        rlwe_scale_RNS_rlwe_addto(out[k1 + (size/m)*k2], h[i][k1], ws[(w_power*i*k1 + w_power*(size/m)*i*k2)&maskMod]);
      }
    }
  }
  for (size_t i = 0; i < m; i++) free_RNS_rlwe_array(size/m, h[i]);
}

typedef struct _ntt_args{
  uint64_t size, w_power, maskMod, * ws;
  RNS_RLWE * out, * in; 
} ntt_args;


void * ntt_thread(void * arg_p){
  ntt_args * arg = (ntt_args *) arg_p;
  assert(arg->out != arg->in);
  for (size_t i = 0; i < arg->size; i++){
    rlwe_RNS_trivial_sample_of_zero(arg->out[i]);
    for (size_t j = 0; j < arg->size; j++){
      rlwe_scale_RNS_rlwe_addto(arg->out[i], arg->in[j], arg->ws[(arg->w_power*i*j)&arg->maskMod]);
    }
  }
  return NULL;
}


void _ntt_mt(RNS_RLWE * out, RNS_RLWE * in, uint64_t size, uint64_t * ws, uint64_t w_power, uint64_t depth, uint64_t maskMod, uint64_t m, uint64_t ntt_max_depth){
  // printf("_NTT: %lu, %lu, %lu, %lu\n", depth, ntt_max_depth, m, size);
  pthread_t threads[MAX_THREADS];
  ntt_args args[MAX_THREADS];
  if(depth == ntt_max_depth || (size%m)){
    assert(out != in);
    for (size_t i = 0; i < size; i++){
      rlwe_RNS_trivial_sample_of_zero(out[i]);
      for (size_t j = 0; j < size; j++){
        rlwe_scale_RNS_rlwe_addto(out[i], in[j], ws[(w_power*i*j)&maskMod]);
      }
    }
    return;
  }
  RNS_RLWE g[m][size/m];
  RNS_RLWE * h[m];
  for (size_t i = 0; i < m; i++){
    for (size_t j = 0; j < size/m; j++){
      g[i][j] = in[m*j + i];
    }
    h[i] = rlwe_alloc_RNS_sample_array2(size/m, out[0]);
    args[i].in = g[i];
    args[i].out = h[i];
    args[i].size = size/m;
    args[i].ws = ws;
    args[i].w_power = m*w_power;
    args[i].maskMod = maskMod;
    pthread_create(&threads[i], NULL, *ntt_thread, (void *) &(args[i]));
  }
  for (size_t i = 0; i < m; i++) pthread_join(threads[i], NULL);

  for (size_t k1 = 0; k1 < size/m; k1++){
    for (size_t k2 = 0; k2 < m; k2++){
      rlwe_RNS_trivial_sample_of_zero(out[k1 + (size/m)*k2]);
      for (size_t i = 0; i < m; i++){
        rlwe_scale_RNS_rlwe_addto(out[k1 + (size/m)*k2], h[i][k1], ws[(w_power*i*k1 + w_power*(size/m)*i*k2)&maskMod]);
      }
    }
  }
  for (size_t i = 0; i < m; i++) free_RNS_rlwe_array(size/m, h[i]);
}

void ntt(RNS_RLWE * out, RNS_RLWE * in, uint64_t root_of_unity, uint64_t p, uint64_t size, uint64_t ntt_max_depth, uint64_t m){
  uint64_t rou[size];
  rou[0] = 1;
  for (size_t i = 1; i < size; i++){
    rou[i] = (rou[i - 1]*root_of_unity)%p;
  }
  _ntt(out, in, size, rou, 1, 1, size - 1, m, ntt_max_depth);
}

void ntt_mt(RNS_RLWE * out, RNS_RLWE * in, uint64_t root_of_unity, uint64_t p, uint64_t size, uint64_t ntt_max_depth, uint64_t m){
  uint64_t rou[size];
  rou[0] = 1;
  for (size_t i = 1; i < size; i++){
    rou[i] = (rou[i - 1]*root_of_unity)%p;
  }
  _ntt_mt(out, in, size, rou, 1, 1, size - 1, m, ntt_max_depth);
}

void rlwe_extract_const_to_tlwe(TLWE out, RNSc_RLWE in){
  const uint64_t N = in->a->N;
  const uint64_t q = in->a->ntt[0]->GetModulus();
  assert(in->a->l == 1);
  out->a[0] = mod_switch(in->a->coeffs[0][0], q, 0);
  for (size_t i = 1; i < N; i++){
    out->a[i] = mod_switch(intel::hexl::SubUIntMod(0, in->a->coeffs[0][N - i], q), q, 0);
  }
  out->b = mod_switch(in->b->coeffs[0][0], q, 0);
}

void alphav_to_polynomials(RNS_Polynomial * out, uint64_t * alphav, uint64_t D){
  for (size_t k = 0; k < D; k++){
    for (size_t i = 0; i < out[k]->l; i++){
      const uint64_t p = out[k]->ntt[i]->GetModulus();
      out[k]->coeffs[i][0] = (p + alphav[k])%p;
      for (size_t j = 1; j < out[k]->N; j++){
        out[k]->coeffs[i][j] = (p - alphav[(out[k]->N - j)*D + k])%p;
      }
      out[k]->ntt[i]->ComputeForward(out[k]->coeffs[i], out[k]->coeffs[i], 1, 1);
    }
  }
}

void rlwe_inner_prod(TLWE * out, RNS_RLWE in, uint64_t ** alpha, uint64_t size){
  const uint64_t l = in->a->l;
  RNS_Polynomial enc = polynomial_new_RNS_polynomial(in->a->N, l, in->a->ntt);
  RNS_RLWE tmp = rlwe_new_RNS_trivial_sample_of_zero(in->a->N, l, in->a->ntt);
  for (size_t i = 0; i < size; i++){
    array_conj_to_RNS(enc, alpha[i]);
    rlwe_RNS_mul_by_poly(tmp, in, enc);
    rlwe_RNS_to_RNSc((RNSc_RLWE) tmp, tmp);
    for (size_t j = 0; j < l - 1; j++){
      polynomial_base_reduce_RNSc_wo_free((RNSc_Polynomial) tmp->a);
      polynomial_base_reduce_RNSc_wo_free((RNSc_Polynomial) tmp->b);
    }
    rlwe_extract_const_to_tlwe(out[i], (RNSc_RLWE) tmp);
    tmp->a->l = l; tmp->b->l = l;
  }
  free_RNS_polynomial(enc);
  free_RNS_rlwe_sample(tmp);
}

void rlwe_inner_prod2(TLWE * out, RNS_RLWE in, RNS_Polynomial * alpha_poly, uint64_t size){
  const uint64_t l = in->a->l;
  RNS_RLWE tmp = rlwe_new_RNS_trivial_sample_of_zero(in->a->N, l, in->a->ntt);
  for (size_t i = 0; i < size; i++){
    rlwe_RNS_mul_by_poly(tmp, in, alpha_poly[i]);
    rlwe_RNS_to_RNSc((RNSc_RLWE) tmp, tmp);
    for (size_t j = 0; j < l - 1; j++){
      polynomial_base_reduce_RNSc_wo_free((RNSc_Polynomial) tmp->a);
      polynomial_base_reduce_RNSc_wo_free((RNSc_Polynomial) tmp->b);
    }
    rlwe_extract_const_to_tlwe(out[i], (RNSc_RLWE) tmp);
    tmp->a->l = l; tmp->b->l = l;
  }
  free_RNS_rlwe_sample(tmp);
}

extern TLWE_Key __glb_tlwe_key;
#define DECRYPT(X) (mod_switch(tlwe_phase(X, __glb_tlwe_key), 0, 65537))
TLWE * batch_codewords(RNS_RLWE * in, uint64_t size, uint64_t * alphav, uint64_t D, uint64_t num_codewords){
  TLWE * res = tlwe_alloc_sample_array(size*D, in[0]->a->N);
  RNS_Polynomial * alpha_poly = polynomial_new_RNS_polynomial_array(D, in[0]->a->N, in[0]->a->l, in[0]->a->ntt);
  alphav_to_polynomials(alpha_poly, alphav,  D);
  // printf("Batched: \n");
  for (size_t i = 0; i < size; i++){
    rlwe_inner_prod2(&res[i*D], in[i], alpha_poly, D);
    // printf("[%lu]: ", i);
    // for (size_t j = 0; j < D; j++) printf("%lu, ", DECRYPT(res[i*D + j]));
    // printf("\n");
  }
  free_RNS_polynomial_array(D, alpha_poly);
  return res;
}

typedef struct _batch_mt_args{
  TLWE * res;
  RNS_RLWE * in;
  RNS_Polynomial * alpha_poly;
  uint64_t num_codewords, D, start, end;
} batch_mt_args;

void * batch_codewords_thread(void * arg_p){
  const batch_mt_args * args = (batch_mt_args *) arg_p;
  for (size_t i = args->start; i < args->end; i++){
    rlwe_inner_prod2(&args->res[i*args->D], args->in[i], args->alpha_poly, args->D);
  }
  return NULL;
}

TLWE * batch_codewords_mt(RNS_RLWE * in, uint64_t size, uint64_t * alphav, uint64_t D, uint64_t num_codewords, uint64_t num_threads){
  TLWE * res = tlwe_alloc_sample_array(size*D, in[0]->a->N);
  RNS_Polynomial * alpha_poly = polynomial_new_RNS_polynomial_array(D, in[0]->a->N, in[0]->a->l, in[0]->a->ntt);
  alphav_to_polynomials(alpha_poly, alphav,  D);
  pthread_t threads[MAX_THREADS];
  batch_mt_args args[MAX_THREADS];
  assert(num_threads <= MAX_THREADS);
  assert(size%num_threads == 0);
  const uint64_t n_points = size/num_threads;
  for (size_t i = 0; i < num_threads; i++){
    args[i].res = res; args[i].in = in; 
    args[i].D = D; args[i].alpha_poly = alpha_poly, args[i].num_codewords = num_codewords;
    args[i].start = i*n_points;
    args[i].end = (i+1)*n_points;
    pthread_create(&threads[i], NULL, *batch_codewords_thread, (void *) &(args[i]));
  }
  for (size_t i = 0; i < num_threads; i++) pthread_join(threads[i], NULL);
  free_RNS_polynomial_array(D, alpha_poly);
  return res;
}