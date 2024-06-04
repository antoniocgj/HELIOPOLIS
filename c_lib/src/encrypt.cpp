#include <fri_c.h>

Codeword new_empty_codeword(uint64_t size, uint64_t N){
  Codeword codeword = (Codeword) safe_malloc(sizeof(*codeword));
  codeword->parts = trlwe_alloc_new_DFT_sample_array(ceil(((float)size)/N), 1, N);
  codeword->size = size;
  return codeword;
}

Codeword new_empty_repacked_codeword(uint64_t size, uint64_t k, uint64_t N){
  Codeword codeword = (Codeword) safe_malloc(sizeof(*codeword));
  codeword->parts = trlwe_alloc_new_DFT_sample_array(size, k, N);
  codeword->size = size;
  return codeword;
}

Codeword encrypt_codeword(uint64_t * vector, uint64_t size, PrivateKey sk){
  const uint64_t N = sk->pk->N;
  Codeword res = new_empty_codeword(size, N);
  TRLWE tmp = trlwe_alloc_new_sample(1, N);
  for (size_t i = 0; i < ceil(((float)size)/N); i++){
    trlwe_sample(tmp, NULL, sk->input_key);
    for (size_t j = 0; j < N; j++){
      tmp->b->coeffs[j] += mod_switch(vector[i*N + j], sk->pk->p, 0);
    }
    trlwe_to_DFT(res->parts[i], tmp);
  }
  free_trlwe(tmp);
  return res;
}

RNS_RLWE * encrypt_input(uint64_t * vector, uint64_t input_size, uint64_t expansion, PrivateKey sk){
  const uint64_t N = sk->pk->N;
  assert(N == sk->rns_rlwe_key->N);
  assert(input_size%N == 0);
  RNS_RLWE * res = (RNS_RLWE *) safe_malloc(((input_size*expansion)/N)*sizeof(RNS_RLWE));
  for (size_t i = 0; i < input_size/N; i++){
    res[i] = rlwe_new_RNS_sample(sk->rns_rlwe_key, &vector[i*N], sk->pk->p);
  }
  for (size_t i = input_size/N; i < ((input_size*expansion)/N); i++){
    res[i] = rlwe_new_RNS_trivial_sample_of_zero(sk->rns_rlwe_key->N,sk->rns_rlwe_key->l,sk->rns_rlwe_key->s_RNS->ntt);
  }
  return res;
}

void mprlwe_to_trlwe(TRLWE_DFT out, RNS_RLWE in){
  const uint64_t l = in->a->l;
  TRLWE tmp2 = trlwe_alloc_new_sample(out->k, out->b->N);
  RNSc_RLWE tmp = rlwe_alloc_RNSc_sample(in->a->N, l, in->a->ntt);
  rlwe_RNS_to_RNSc(tmp, in);
  const uint64_t N = out->b->N;
  assert(N == in->b->N);
  assert(out->k == 1);
  for (size_t j = 0; j < l - 1; j++){
    polynomial_base_reduce_RNSc(tmp->a);
    polynomial_base_reduce_RNSc(tmp->b);
  }
  const uint64_t q = in->a->ntt[0]->GetModulus();
  for (size_t j = 0; j < N; j++){
    tmp2->a[0]->coeffs[j] = mod_switch(tmp->a->coeffs[0][j], q, 0);
    tmp2->b->coeffs[j] = mod_switch(tmp->b->coeffs[0][j], q, 0);
  }
  trlwe_to_DFT(out, tmp2);
  free_rlwe_RNS_sample(tmp);
  free_trlwe(tmp2);
}

void decrypt_and_print(RNS_RLWE * c, uint64_t input_size, uint64_t coeff_size, PrivateKey sk){
  const uint64_t N = sk->pk->N, p = sk->pk->p;
  TRLWE_DFT trlwe = trlwe_alloc_new_DFT_sample(1, N);
  TorusPolynomial poly = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < input_size; i++){
    mprlwe_to_trlwe(trlwe, c[i]);
    trlwe_DFT_phase(poly, trlwe, sk->input_key);
    printf("poly[%lu]: ", i);
    for (size_t j = 0; j < coeff_size; j++){
      const uint64_t val = mod_switch(torus2int(poly->coeffs[j], 52), 1ULL<<52, p)%p; 
      printf("%lu, ", val);
    }
    printf("\n");
  }
  free_trlwe(trlwe);
  free_polynomial(poly);
}

Codeword repack_init_codeword(RNS_RLWE * first_codeword, uint64_t n){
  Codeword codeword = new_empty_repacked_codeword(n, 1, first_codeword[0]->b->N);
  for (size_t i = 0; i < n; i++){
    mprlwe_to_trlwe(codeword->parts[i], first_codeword[i]);
  }
  return codeword;
}


void decrypt_poly(uint64_t * out, RNS_RLWE * c, uint64_t input_size, uint64_t idx, PrivateKey sk){
  const uint64_t N = sk->pk->N, p = sk->pk->p;
  TRLWE_DFT trlwe = trlwe_alloc_new_DFT_sample(1, N);
  TorusPolynomial poly = polynomial_new_torus_polynomial(N);
  for (size_t i = 0; i < input_size; i++){
    mprlwe_to_trlwe(trlwe, c[i]);
    trlwe_DFT_phase(poly, trlwe, sk->input_key);
    const uint64_t val = mod_switch(poly->coeffs[idx], 0, p); 
    out[i] = val;
  }
  free_trlwe(trlwe);
  free_polynomial(poly);
}


PublicKey new_public_key(PrivateKey sk, uint64_t D, uint64_t p, uint64_t decomp_bits, uint64_t t_decomp, uint64_t b_decomp, uint64_t t_dimred){
  PublicKey res = (PublicKey) safe_malloc(sizeof(*res));
  res->N = sk->input_key->s[0]->N;
  res->D = D;
  res->p = p;
  res->p_log = (uint64_t) (log2(p) + 1);
  res->packing_ksk = trlwe_new_full_packing_KS_key(sk->packing_key, sk->extracted_key, t_decomp, b_decomp);
  sk->pk = res;
  return res;
}

PublicKey copy_public_key(PublicKey pk){
  return pk; // no copy needed
  PublicKey res = (PublicKey) safe_malloc(sizeof(*res));
  res->N = pk->N;
  res->D = pk->D;
  res->p = pk->p;
  res->p_log = pk->p_log;
  res->decomp_bits = pk->decomp_bits;
  res->tmp = tlwe_alloc_sample_array(pk->N, pk->N);
  res->tmp_mrlwe = trlwe_alloc_new_sample(2, pk->tmp_mrlwe->b->N);
  res->rlwe_ksk = pk->rlwe_ksk;
  res->packing_ksk = pk->packing_ksk;
  return res;
}

TLWE_Key __glb_tlwe_key;
PrivateKey new_private_key(uint64_t in_N, uint64_t in_digits, uint64_t in_digit_size, double in_sigma_err, double in_sigma_key, uint64_t packing_N, uint64_t packing_k, double packing_sigma){
  PrivateKey res = (PrivateKey) safe_malloc(sizeof(*res));
  // input parameters
  // const uint64_t in_N = 8192;
  // const uint64_t digit_size = 49;
  // const uint64_t in_digits = 3; // ciphertext modulus q = 2^{in_digits*digit_size}
  // const double in_sigma =  1;
  const double in_sigma_tfhe = pow(2,-53);
  /* input key */
  // gen single precision key of size N and Hamming weight h
  // const uint64_t in_h = 128; // hamming weight ternary key
  // res->input_key = trlwe_new_ternary_key(in_N, 1, in_h, in_sigma_tfhe); 
  res->input_key = trlwe_new_gaussian_key(in_N, 1, in_sigma_key, in_sigma_tfhe);
  // gen RNS parameters
  auto Q = intel::hexl::GeneratePrimes(in_digits, in_digit_size, true, in_N);
  assert(Q.size() == in_digits);
  auto ntt = new_ntt_list(Q.data(), in_N, in_digits);
  res->rns_rlwe_key = rlwe_get_RNS_key_from_array(in_N, in_digits, res->input_key->s[0]->coeffs, ntt, in_sigma_err);
  // extract lwe interpretation
  res->extracted_key = tlwe_alloc_key(in_N, in_sigma_tfhe); 
  trlwe_extract_tlwe_key(res->extracted_key, res->input_key);
  __glb_tlwe_key = res->extracted_key; // key for debugging
  
  /* (re)packing key */
  res->packing_key = trlwe_new_ternary_key(packing_N, packing_k, packing_N/4, packing_sigma);

  /* output key (if different from packing) */
  res->output_key = res->packing_key;
  return res;
}
