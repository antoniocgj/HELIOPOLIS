#include <fri_c.h>

void repack(TRLWE_DFT out, TLWE * c, uint64_t size, PublicKey params){
  TRLWE tmp = trlwe_alloc_new_sample(out->k, out->b->N);
  trlwe_full_packing_keyswitch(tmp, c, size, params->packing_ksk);
  trlwe_to_DFT(out, tmp);
  free_trlwe(tmp);
}

void decomp_repack(TRLWE_DFT out, TLWE * c, uint64_t size, PublicKey params){
  TRLWE tmp = trlwe_alloc_new_sample(out->k, out->b->N);
  const uint64_t k = (uint64_t) ceil(((float)params->p_log)/(params->decomp_bits));
  const uint64_t base = params->decomp_bits;
  // decomposition
  for (size_t i = 0; i < size; i++){
    for (size_t j = 0; j < k; j++){
      tlwe_scale(params->tmp[i*k + j], c[i], 1ULL<<(base*j));
    }
  }
  // packing
  trlwe_full_packing_keyswitch(params->tmp_mrlwe, params->tmp, size*k, params->packing_ksk);
  // dim. reduction
  params->tmp_mrlwe->k = 1;
  trlwe_keyswitch(tmp, params->tmp_mrlwe, params->rlwe_ksk);
  params->tmp_mrlwe->k = 2;
  polynomial_addto_torus_polynomial(tmp->a[0], params->tmp_mrlwe->a[1]);
  trlwe_to_DFT(out, tmp);
  free_trlwe(tmp);
}

uint64_t getMSB(uint64_t x, uint64_t bit){
  return x>>(64-bit);
}

uint64_t signExtend(uint64_t x, uint64_t bit){
  return (uint64_t)(((int64_t)(x << (64-bit)))>>(64-bit));
}

uint64_t mod_switch_noise(uint64_t v, uint64_t p, uint64_t q) {
  const double double_q = q == 0? pow(2,64) : ((double)q); 
  double val = round((((double) v) * double_q) / ((double) p));
  const double res = val < p? val : val - p;
  const double v2 = round(res * ((double) p)) / double_q;
  const double res_v2 = v2 < double_q? v2 : v2 - double_q;
  const double diff = fabs(((double) v) - res_v2);
  const double diff_mod = diff > double_q/2? double_q - diff : diff;
  return diff_mod;
}

void recomp(uint64_t * out, uint64_t size, TRLWE_DFT in, PrivateKey sk){
  const uint64_t k = (uint64_t) ceil(((float)sk->pk->p_log)/(sk->pk->decomp_bits));
  const uint64_t log_base = sk->pk->decomp_bits, log_mask = (1ULL<<(log_base+2)) -1;
  // decrypt
  trlwe_DFT_phase(sk->tmp, in, sk->output_key);
  // recomposition
  for (size_t i = 0; i < size; i++){
    uint64_t rb = 0;
    out[i] = torus2int(sk->tmp->coeffs[i*k + k - 1], log_base);
    rb = out[i] >> (log_base - 2);
    for (int64_t j = k - 2; j >= 0; j--){
      const uint64_t phase_slice = getMSB(sk->tmp->coeffs[i*k + j], log_base + 2);
      const uint64_t slice = (phase_slice + signExtend(rb - (phase_slice&0x3), 2))&log_mask;
      out[i] |=  (slice >> 2) << (log_base*(k - j - 1));
      rb = slice >> log_base;
    }
    out[i] = mod_switch(out[i], 1ULL << (k*log_base), sk->pk->p);
  }
}

void recomp_and_print(uint64_t size, Codeword in, uint64_t index, PrivateKey sk){
  uint64_t out[size];
  recomp(out, size, in->parts[index], sk);
  printf("Recomp result: [");
  for (size_t i = 0; i < size - 1; i++){
    printf("%lu, ", out[i]);
  }
  printf("%lu]\n", out[size - 1]);
}

TRLWE_DFT * codeword_get_packed_points(Codeword cw, uint64_t * index, uint64_t size){
  TRLWE_DFT * res = (TRLWE_DFT *) safe_malloc(sizeof(TRLWE_DFT)*size);
  for (size_t i = 0; i < size; i++){
    assert(index[i] < cw->size);
    res[i] = cw->parts[index[i]];
  }
  return res;
}

uint64_t * decrypt_recomp_packed_points(TRLWE_DFT * points, uint64_t size, PrivateKey sk){
  uint64_t * out = (uint64_t *) safe_malloc(sizeof(uint64_t)*size*sk->pk->D*2);
  for (size_t i = 0; i < size; i++){
    recomp(&out[i*sk->pk->D*2], sk->pk->D*2, points[i], sk);
  }
  return out;
}

void decrypt_packed_points(uint64_t * out, TRLWE_DFT * points, uint64_t size, PrivateKey sk){
  TorusPolynomial tmp = polynomial_new_torus_polynomial(points[0]->b->N);
  const uint64_t p = sk->pk->p;
  for (size_t i = 0; i < size; i++){
    trlwe_DFT_phase(tmp, points[i], sk->packing_key);
    for (size_t j = 0; j < 2*sk->pk->D; j++){
      const uint64_t val = mod_switch(tmp->coeffs[j], 0, p); // mod_switch(torus2int(tmp->coeffs[j], 52), 1ULL<<52, p)%p; 
      out[i*sk->pk->D*2 + j] = val; 
    }
  }
  free_polynomial(tmp);
}

void decrypt_and_batch_points2(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t ** alphav, PrivateKey sk, uint64_t pack_size){
  const uint64_t N = points[0]->b->N, D = sk->pk->D;
  TorusPolynomial tmp = polynomial_new_torus_polynomial(N);
  const uint64_t p = sk->pk->p;
  assert(pack_size <= N);
  for (size_t i = 0; i < size; i++){
    trlwe_DFT_phase(tmp, points[i], sk->input_key);
    memset(&out[i*D], 0, sizeof(uint64_t)*D);
    for (size_t j = 0; j < pack_size; j++){
      const uint64_t val = mod_switch(tmp->coeffs[j], 0, p);
      for (size_t k = 0; k < D; k++){
        out[i*D + k] = (out[i*D + k]+val*alphav[j][k])%p; 
      }      
    }
  }
  free_polynomial(tmp);
}

void decrypt_and_batch_points(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t * alphav, PrivateKey sk){
  const uint64_t N = points[0]->b->N, D = sk->pk->D;
  TorusPolynomial tmp = polynomial_new_torus_polynomial(N);
  const uint64_t p = sk->pk->p;
  for (size_t i = 0; i < size; i++){
    trlwe_DFT_phase(tmp, points[i], sk->input_key);
    memset(&out[i*D], 0, sizeof(uint64_t)*D);
    #ifdef AVX512_OPT
      __m512i * outv = (__m512i *) out;
      __m512i * alphavv = (__m512i *) alphav;
      for (size_t j = 0; j < N; j++){
        const uint64_t val = tmp->coeffs[j];
        const __m512i valv = _mm512_set1_epi64(val>>12);
        outv[i*2] = _mm512_madd52lo_epu64(outv[i*2], valv, alphavv[j*2]);
        outv[i*2 + 1] = _mm512_madd52lo_epu64(outv[i*2 + 1], valv, alphavv[j*2 + 1]);
      }
    #else
    for (size_t j = 0; j < N; j++){
      const uint64_t val = mod_switch(tmp->coeffs[j], 0, p);
      intel::hexl::EltwiseFMAMod(&out[i*D], &alphav[j*D], val, &out[i*D], D, p, 1);
    }
    #endif
  }
  #ifdef AVX512_OPT
  for (size_t j = 0; j < D*size; j++){
    out[j] = mod_switch(out[j]<<12, 0, p);
  }
  #endif
  free_polynomial(tmp);
}

typedef struct _decrypt_batch_arg{
  uint64_t start, end, * alphav, * out;
  PrivateKey sk;
  TRLWE_DFT * points;
} decrypt_batch_arg;


void * decrypt_and_batch_thread(void * arg_p){
  decrypt_batch_arg * arg = (decrypt_batch_arg*) arg_p;
  const uint64_t N = arg->points[0]->b->N, D = arg->sk->pk->D;
  TorusPolynomial tmp = polynomial_new_torus_polynomial(N);
  for (size_t i = arg->start; i < arg->end; i++){
    trlwe_DFT_phase(tmp, arg->points[i], arg->sk->input_key);
    memset(&arg->out[i*D], 0, sizeof(uint64_t)*D);
    #ifdef AVX512_OPT
      __m512i * outv = (__m512i *) arg->out;
      __m512i * alphavv = (__m512i *) arg->alphav;
      for (size_t j = 0; j < N; j++){
        const uint64_t val = tmp->coeffs[j];
        const __m512i valv = _mm512_set1_epi64(val>>12);
        outv[i*2] = _mm512_madd52lo_epu64(outv[i*2], valv, alphavv[j*2]);
        outv[i*2 + 1] = _mm512_madd52lo_epu64(outv[i*2 + 1], valv, alphavv[j*2 + 1]);
      }
    #else
    const uint64_t p = arg->sk->pk->p;
    for (size_t j = 0; j < N; j++){
      const uint64_t val = mod_switch(tmp->coeffs[j], 0, p);
      intel::hexl::EltwiseFMAMod(&arg->out[i*D], &arg->alphav[j*D], val, &arg->out[i*D], D, p, 1);
    }
    #endif
  }
  free_polynomial(tmp);
  return NULL;
}

void decrypt_and_batch_points_mt(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t * alphav, PrivateKey sk, uint64_t num_threads){
  pthread_t threads[MAX_THREADS];
  decrypt_batch_arg args[MAX_THREADS];

  assert(num_threads <= MAX_THREADS);
  const uint64_t n_points = size/num_threads;
  for (size_t i = 0; i < num_threads; i++){
    args[i].alphav = alphav; 
    args[i].out = out; 
    args[i].sk = sk; 
    args[i].points = points; 
    args[i].start = i*n_points; 
    args[i].end = (i+1)*n_points;
    if(i == num_threads-1) args[i].end = size;
    pthread_create(&threads[i], NULL, *decrypt_and_batch_thread, (void *) &(args[i]));
  }
  for (size_t i = 0; i < num_threads; i++) pthread_join(threads[i], NULL);
  #ifdef AVX512_OPT
  const uint64_t p = sk->pk->p, D = sk->pk->D;
  for (size_t j = 0; j < D*size; j++){
    out[j] = mod_switch(out[j]<<12, 0, p);
  }
  #endif
}

void decrypt_codeword(uint64_t * out, Codeword c, PrivateKey sk){
  TorusPolynomial tmp = polynomial_new_torus_polynomial(c->parts[0]->b->N);
  const uint64_t p = sk->pk->p;
  for (size_t i = 0; i < c->size; i++){
    trlwe_DFT_phase(tmp, c->parts[i], sk->packing_key);
    for (size_t j = 0; j < 2*sk->pk->D; j++){
      const uint64_t val = mod_switch(torus2int(tmp->coeffs[j], 52), 1ULL<<52, p)%p; 
      out[i*sk->pk->D*2 + j] = val; 
    }
  }
  free_polynomial(tmp);
}

void decrypt_codeword2(uint64_t ** out, Codeword c, PrivateKey sk){
  TorusPolynomial tmp = polynomial_new_torus_polynomial(c->parts[0]->b->N);
  const uint64_t p = sk->pk->p, D = sk->pk->D;
  for (size_t i = 0; i < c->size; i++){
    trlwe_DFT_phase(tmp, c->parts[i], sk->packing_key);
    out[i] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*D);
    for (size_t j = 0; j < D; j++){
      const uint64_t val = mod_switch(torus2int(tmp->coeffs[j], 52), 1ULL<<52, p)%p; 
      out[i][j] = val; 
    }
    out[i+c->size] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*D);
    for (size_t j = 0; j < D; j++){
      const uint64_t val = mod_switch(torus2int(tmp->coeffs[D + j], 52), 1ULL<<52, p)%p; 
      out[i+c->size][j] = val; 
    }
  }
  free_polynomial(tmp);
}

void codeword_get_point(TLWE out, Codeword cw, uint64_t index){
  assert(false); // not implemented
  // const uint64_t N = out->n, mask = N - 1;
  // trlwe_extract_tlwe(out, cw->parts[index/N], index&mask);
}

void * repack_codeword(Codeword first_codeword, uint64_t n, PublicKey params){
  const uint64_t N = params->N, D = params->D;
  TLWE * acc = tlwe_alloc_sample_array(2*D, N);
  Codeword codeword = new_empty_repacked_codeword(n/2, 1, 512);
  for (size_t di = 0; di < 2*D; di++) tlwe_noiseless_trivial_sample(acc[di], 0); // acc <- 0
  for (size_t i = 0; i < n/2; i++){
    // +w^i
    codeword_get_point(acc[0], first_codeword, i);
    // -w^i
    codeword_get_point(acc[D], first_codeword, n/2 + i);
    // repacking
    decomp_repack(codeword->parts[i], acc, 2*D, params);
  }
  free_tlwe_array(acc,2*D);
  return codeword;
}

typedef struct _repack_args{
  Codeword out, first_codeword;
  uint64_t n, start, end;
  PublicKey params;
} repack_args;

void * repack_codeword_thread(void * arg_p){
  const repack_args * args = (repack_args *) arg_p;
  const uint64_t N = args->params->N, D = args->params->D, n = args->n;
  TLWE * acc = tlwe_alloc_sample_array(2*D, N);
  for (size_t di = 0; di < 2*D; di++) tlwe_noiseless_trivial_sample(acc[di], 0); // acc <- 0
  for (size_t i = args->start; i < args->end; i++){
    // +w^i
    codeword_get_point(acc[0], args->first_codeword, i);
    // -w^i
    codeword_get_point(acc[D], args->first_codeword, n/2 + i);
    // repacking
    decomp_repack(args->out->parts[i], acc, 2*D, args->params);
  }
  free_tlwe_array(acc,2*D);
  return NULL;
}

void * repack_codeword_mt(Codeword first_codeword, uint64_t n, PublicKey params, uint64_t num_threads){
  pthread_t threads[MAX_THREADS];
  Codeword codeword = new_empty_repacked_codeword(n/2, 1, 512);
  repack_args args[MAX_THREADS];
  const uint64_t n_points = (n/2)/num_threads;
  assert(num_threads <= MAX_THREADS);
  assert((n/2)%num_threads == 0);
  for (size_t i = 0; i < num_threads; i++){
    PublicKey new_params = copy_public_key(params);
    args[i].out = codeword; args[i].first_codeword = first_codeword; 
    args[i].n = n; args[i].params = new_params;
    args[i].start = i*n_points;
    args[i].end = (i+1)*n_points;
    pthread_create(&threads[i], NULL, *repack_codeword_thread, (void *) &(args[i]));
  }
  for (size_t i = 0; i < num_threads; i++) pthread_join(threads[i], NULL);
  return codeword;
}

void * folding_level(Codeword first_codeword, uint64_t ** consts_array, uint64_t n, uint64_t r, PublicKey params) {
  const uint64_t N = params->N, D = params->D;
  Codeword codeword = new_empty_repacked_codeword(n/4, 1, 512);
  TLWE * acc = tlwe_alloc_sample_array(2*D, N);
  TLWE tmp = tlwe_alloc_sample(N);

  for (size_t i = 0; i < n/4; i++){
    for (size_t di = 0; di < 2*D; di++) tlwe_noiseless_trivial_sample(acc[di], 0); // acc <- 0
    // +w^i
    for (size_t j = 0; j < 1ULL<<(r+1); j++){
      codeword_get_point(tmp, first_codeword, i + j*n/2);
      for (size_t di = 0; di < D; di++){
        tlwe_scale_addto(acc[di], tmp, consts_array[i + j*n/2][di]);       
      }
    }
    // -w^i
    for (size_t j = 0; j < 1ULL<<(r+1); j++){
      codeword_get_point(tmp, first_codeword, (n/4) + i + j*n/2);
      for (size_t di = 0; di < D; di++){
        tlwe_scale_addto(acc[di + D], tmp, consts_array[(n/4) + i + j*n/2][di]);
      }
    }
    // repacking
    decomp_repack(codeword->parts[i], acc, 2*D, params);
  }
  free_tlwe(tmp);
  free_tlwe_array(acc,2*D);
  return codeword;
}

// #include <debug_util.h>
extern TLWE_Key __glb_tlwe_key;
extern TRLWE_Key __glb_repack_key;
void * folding_level_LWE(TLWE * first_codeword, uint64_t ** consts_array, uint64_t n, uint64_t r, PublicKey params) {
  const uint64_t N = params->N, D = params->D;
  fp_ini(params->p);
  fp_setup_temps(N);
  Codeword codeword = new_empty_repacked_codeword(n/4, params->packing_ksk->s[0][0]->k, params->packing_ksk->s[0][0]->b->N);
  TLWE * acc_lwe = tlwe_alloc_sample_array(3*D, N);
  fp16_he_t acc[2], tmp, point;
  fp16_t alpha;
  fp16_he_in(acc[0], acc_lwe);
  fp16_he_in(acc[1], &(acc_lwe[D]));
  fp16_he_in(tmp, &(acc_lwe[2*D]));

  for (size_t i = 0; i < n/4; i++){
    for (size_t di = 0; di < 2*D; di++) tlwe_noiseless_trivial_sample(acc_lwe[di], 0); // acc <- 0
    // +w^i
    for (size_t j = 0; j < 1ULL<<(r+1); j++){
      fp16_he_in(point, &(first_codeword[D*(i + j*n/2)]));
      // fp16_he_print(point, __glb_tlwe_key, params->p);
      fp16_in_u64(alpha, consts_array[i + j*n/2]);
      // fp16_print(alpha);
      fp16_mul_he(tmp, alpha, point);
      // fp16_he_print(tmp, __glb_tlwe_key, params->p);
      fp16_add_he(acc[0], acc[0], tmp);
      // fp16_he_print(acc[0], __glb_tlwe_key, params->p);
    }
    // exit(0);
    // -w^i
    for (size_t j = 0; j < 1ULL<<(r+1); j++){
      fp16_he_in(point, &(first_codeword[D*(n/4 + i + j*n/2)]));
      fp16_in_u64(alpha, consts_array[(n/4) + i + j*n/2]);
      fp16_mul_he(tmp, alpha, point);
      fp16_add_he(acc[1], acc[1], tmp);
    }
    // repacking
    repack(codeword->parts[i], acc_lwe, 2*D, params);
  }
  free_tlwe_array(acc_lwe, 3*D);
  return codeword;
}

// multi-threaded version

typedef struct _folding_args{
  TLWE * first_codeword;
  uint64_t n, r, start, end;
  uint64_t ** consts_array;
  TRLWE_DFT * out;
  PublicKey params;
} folding_args;

void * folding_level_LWE_thread(void * arg_p){
  const folding_args * args = (folding_args *) arg_p;
  const uint64_t N = args->params->N, D = args->params->D;
  TLWE * acc_lwe = tlwe_alloc_sample_array(3*D, N);
  fp16_he_t acc[2], tmp, point;
  fp16_t alpha;
  fp16_he_in(acc[0], acc_lwe);
  fp16_he_in(acc[1], &(acc_lwe[D]));
  fp16_he_in(tmp, &(acc_lwe[2*D]));

  fp_ini(args->params->p);
  fp_setup_temps(N);

  for (size_t i = args->start; i < args->end; i++){
    for (size_t di = 0; di < 2*D; di++) tlwe_noiseless_trivial_sample(acc_lwe[di], 0); // acc <- 0
    // +w^i
    for (size_t j = 0; j < 1ULL<<(args->r+1); j++){
      fp16_he_in(point, &(args->first_codeword[D*(i + j*args->n/2)]));
      fp16_in_u64(alpha, args->consts_array[i + j*args->n/2]);
      fp16_mul_he(tmp, alpha, point);
      fp16_add_he(acc[0], acc[0], tmp);
    }
    // -w^i
    for (size_t j = 0; j < 1ULL<<(args->r+1); j++){
      fp16_he_in(point, &(args->first_codeword[D*(args->n/4 + i + j*args->n/2)]));
      fp16_in_u64(alpha, args->consts_array[(args->n/4) + i + j*args->n/2]);
      fp16_mul_he(tmp, alpha, point);
      fp16_add_he(acc[1], acc[1], tmp);
    }
    // repacking
    repack(args->out[i], acc_lwe, 2*D, args->params);
  }
  free_tlwe_array(acc_lwe, 3*D);
  return NULL;
}

void * folding_level_LWE_mt(TLWE * first_codeword, uint64_t ** consts_array, uint64_t n, uint64_t r, PublicKey params, uint64_t num_threads) {
  Codeword codeword = new_empty_repacked_codeword(n/4, params->packing_ksk->s[0][0]->k, params->packing_ksk->s[0][0]->b->N);

  pthread_t threads[MAX_THREADS];
  folding_args args[MAX_THREADS];
  assert(num_threads <= MAX_THREADS);
  assert((n/4)%num_threads == 0);
  const uint64_t n_points = (n/4)/num_threads;
  for (size_t i = 0; i < num_threads; i++){
    args[i].first_codeword = first_codeword; 
    args[i].consts_array = consts_array; 
    args[i].params = params; 
    args[i].n = n; args[i].r = r; 
    args[i].out = codeword->parts;
    args[i].start = i*n_points;
    args[i].end = (i+1)*n_points;
    pthread_create(&threads[i], NULL, *folding_level_LWE_thread, (void *) &(args[i]));
  }
  for (size_t i = 0; i < num_threads; i++) pthread_join(threads[i], NULL);
  return codeword;
}