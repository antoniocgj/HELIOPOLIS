#include <fri_c.h>
#include <iop.hpp>

// functions for the verifier

void sample_field_array2(uint64_t * res, uint64_t * hash_in, uint64_t size, uint64_t p){
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  blake3_hasher_update(&hasher, hash_in, 4*8);
  blake3_hasher_finalize(&hasher, (uint8_t *) res, size*sizeof(uint64_t));
#ifdef AVX512_OPT
  const __m512i mul_mask = _mm512_set1_epi64(p);
  const __m512i zerov = _mm512_setzero_si512();
  __m512i * resv = (__m512i *) res;
  for (size_t i = 0; i < size/8; i++){
    resv[i] = _mm512_madd52hi_epu64 (zerov, resv[i], mul_mask);
  }
#else
  for (size_t i = 0; i < size; i++){
    res[i] = mod_switch(res[i], 0, p);
  }
#endif
}

uint64_t * sample_field_array(uint64_t * hash_in, uint64_t size, uint64_t p){
  uint64_t * res = (uint64_t *) safe_aligned_malloc(size*sizeof(uint64_t));
  sample_field_array2(res, hash_in, size, p);
  return res;
}

unsigned char char_rev(unsigned char b) {
  b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
  b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
  b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
  return b;
}

uint32_t int_rev(uint32_t b){
  uint32_t a = _bswap(b);
  unsigned char * a_vec = (unsigned char *) &a;
  a_vec[0] = char_rev(a_vec[0]);
  a_vec[1] = char_rev(a_vec[1]);
  a_vec[2] = char_rev(a_vec[2]);
  a_vec[3] = char_rev(a_vec[3]);
  return a;
}

void bit_rev(uint64_t ** out, uint64_t ** in, uint64_t n){
  const uint64_t log_n = log2(n);
  for (int i = 0; i < n; i++){
    out[i] = in[int_rev(i)>>(32-log_n)];
  }
}

void calcRootsOfUnit(uint64_t * out, uint64_t min_RoU, uint64_t q, uint64_t n){
  uint64_t r = 1;
  for (size_t i = 0; i < n; i++){
    out[i] = r;
    r = intel::hexl::MultiplyMod(r, min_RoU, q);
  }
}

void _clear_ntt(uint64_t ** out, uint64_t ** in, uint64_t * ws, uint64_t p, uint64_t n, uint64_t D){
  assert(in != out);
  bit_rev(out, in, n);
  uint64_t m = 2;
  uint64_t t[D], u[D];
  while(m <= n){
    const uint64_t w_m = ws[n/m];
    uint64_t w = 1;
    for (size_t j = 0; j < m/2; j++){
      for (size_t k = 0; k < n-1; k+=m){
        intel::hexl::EltwiseFMAMod(t, out[k + j + m/2], w, nullptr, D, p, 1);
        memcpy(u, out[k + j], sizeof(uint64_t)*D);
        intel::hexl::EltwiseAddMod(out[k + j], u, t, D, p);
        intel::hexl::EltwiseSubMod(out[k + j + m/2], u, t, D, p);
      }
      w = intel::hexl::MultiplyMod(w, w_m, p);
    }
    m *= 2;
  }
}

void clear_ntt(uint64_t ** out, uint64_t ** in, uint64_t w, uint64_t p, uint64_t n, uint64_t D){
  uint64_t rous[n];
  calcRootsOfUnit(rous, w, p, n);
  _clear_ntt(out, in, rous, p, n, D);
  const uint64_t inv_n = intel::hexl::InverseMod(n, p);
  for (size_t i = 0; i < n; i++)
    intel::hexl::EltwiseFMAMod(out[i], out[i], inv_n, nullptr, D, p, 1);
}


void sample_indices(uint64_t * out, uint64_t * seed, uint64_t size, uint64_t reduced_size, uint64_t number){
  assert(number <= reduced_size);
  const uint64_t mod_mask = size - 1;
  const uint64_t mod_mask_red = reduced_size - 1;
  uint64_t counter = 0, sampled = 0, hash_seek = 1024*sizeof(uint32_t);
  uint8_t reduced_idx[reduced_size];
  memset(reduced_idx, 0, reduced_size);
  uint32_t hash_buffer[1024];
  // sample hash
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  blake3_hasher_update(&hasher, seed, 32);
  blake3_hasher_finalize(&hasher, (uint8_t *) hash_buffer, hash_seek);
  // sample indices
  while (sampled < number){
    const uint64_t idx = hash_buffer[counter++]&mod_mask;
    const uint64_t red_idx = idx & mod_mask_red;
    if(reduced_idx[red_idx] == 0){
      out[sampled++] = idx;
      reduced_idx[red_idx] = 1;
    }
    if(counter == 1024){ // gen more numbers
      blake3_hasher_finalize_seek(&hasher, hash_seek, (uint8_t *) hash_buffer, 1024*sizeof(uint32_t));
      counter = 0;
      hash_seek += 1024*sizeof(uint32_t);
    }
  }
}

bool run_lin_check_layer(uint32_t * points,  uint32_t * next_points, uint64_t * indices, uint64_t D, uint64_t num_tests, uint64_t num_rounds, uint64_t * inv_omega_powers, uint64_t ** alphas){
  fp16_t alpha, ay, by, cay, cby, inv_two, inv_omega_power, t0, t1;
  fp16_in_const64(t0, 2);
  fp16_inv(inv_two, t0);

  uint64_t omega_index = 1;
  for (size_t r = 0; r < num_rounds; r++)
  {
    fp16_in_u64(alpha, alphas[r]);
    for (size_t i = 0; i < num_tests; i++){
      fp16_in(ay, &next_points[i*2*D]);
      fp16_in(by, &next_points[i*2*D + D]);
      fp16_in(cay, &points[i*2*D]);
      fp16_in(cby, &points[i*2*D + D]);

      fp16_in_const64(inv_omega_power, inv_omega_powers[omega_index*indices[i]]);

      fp16_add(t0, ay, by); // (ay + by)
      fp16_sub(t1, ay, by); // (ay - by)
      fp16_mul(t1, t1, alpha); // (ay - by) * alpha
      fp16_mul(t1, t1, inv_omega_power); // (ay - by) * alpha * inv_omega_power
      fp16_add(t0, t0, t1); // (ay + by) + (ay - by) * alpha * inv_omega_power
      fp16_mul(t0, t0, inv_two); // inv_two * ((ay + by) + (ay - by) * alpha * inv_omega_power)

      if(fp16_cmp(t0, cay) + fp16_cmp(t0, cby) == 0){
        return false;
      }
    }
  }
  return true;
}

bool check_degree(uint64_t ** in, uint64_t D, uint64_t max_degree, uint64_t size){
  for (size_t i = max_degree; i < size; i++){
    for (size_t j = 0; j < D; j++){
      if(in[i][j] != 0) return false;
    }
  }
  return true;
}

uint64_t * gen_inv_powers(uint64_t omega, uint64_t size, uint64_t p){
  uint64_t * res = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*size);
  res[0] = 1;
  const uint64_t inv_omega = intel::hexl::InverseMod(omega, p);
  for (size_t i = 1; i < size; i++){
    res[i] = intel::hexl::MultiplyMod(res[i-1], inv_omega, p);
  }
  return res;
}

bool compare_hash(uint64_t a[4], uint64_t b[4]){
  return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]) && (a[3] == b[3]);
}

void print_vec_Fpd(char * msg, uint64_t ** vec, uint64_t size){
  printf("%s: [", msg);
  for (size_t i = 0; i < size; i++){
    printf("[");
    for (size_t j = 0; j < 15; j++){
      printf("%lu, ", vec[i][j]);
    }    
    printf("%lu], ", vec[i][15]);
  }
  printf("]\n");
}

#include <benchmark_util.h>

bool verify_benchmark(void * proof_stream_p, uint64_t size, uint64_t num_tests, uint64_t num_rounds, uint64_t omega, uint64_t expansion, PrivateKey sk){
  bool res;
  MEASURE_TIME("", 100, "\n Verifier time:", 
    res = verify(proof_stream_p, size, num_tests, num_rounds, omega, expansion, sk);
  );
  return res;
}

bool verify(void * proof_stream_p, uint64_t size, uint64_t num_tests, uint64_t num_rounds, uint64_t omega, uint64_t expansion, PrivateKey sk){
  ProofStream * proof_stream = (ProofStream *) proof_stream_p;
  proof_stream->reset_read();
  uint64_t ** merkle_roots = (uint64_t **) safe_malloc(sizeof(uint64_t *)*num_rounds);
  uint64_t ** alphas = (uint64_t **) safe_malloc(sizeof(uint64_t *)*num_rounds);


  const uint64_t D = sk->pk->D, p = sk->pk->p, N = sk->pk->N;
  uint64_t fs_hash[4], idx_mod_mask = (size - 1), omega_power = 1;
  fp16_t alpha, ay, by, cay, cby, inv_two, inv_omega_power, t0, t1;
  uint64_t * inv_omega_powers = gen_inv_powers(omega, size, p);
  fp16_in_const64(t0, 2);
  fp16_inv(inv_two, t0);

  
  merkle_roots[0] = (uint64_t *) proof_stream->pull();
  proof_stream->read_fiat_shamir(fs_hash);
  alphas[0] = sample_field_array(fs_hash, D, p);

  uint64_t * beta = sample_field_array(fs_hash, N*D, p);

  for (size_t i = 1; i < num_rounds; i++){
    merkle_roots[i] = (uint64_t *) proof_stream->pull();
    proof_stream->read_fiat_shamir(fs_hash);
    alphas[i] = sample_field_array(fs_hash, D, p);
  }

  // test last code word
  Codeword last_codeword = (Codeword) proof_stream->pull();
  const uint64_t reduced_size = last_codeword->size*2;
  assert(reduced_size == size >> (num_rounds - 1));

  // test hash
  uint64_t * last_codeword_hash = hash_full_codeword(last_codeword);
  
  // check last codeword
  if(!compare_hash(merkle_roots[num_rounds - 1], last_codeword_hash)){
    printf("last codeword is not well formed\n");
    return false;
  }
  
  uint64_t points[2*D*num_tests] __attribute__ ((aligned(64)));
  uint64_t points2[2*D*num_tests] __attribute__ ((aligned(64)));
  uint64_t points_hash[2*4*num_tests] __attribute__ ((aligned(64)));

  uint64_t * curr_points = points2, * previous_points = points, * tmp_pointer;

  uint64_t indices[num_tests];
  proof_stream->read_fiat_shamir(fs_hash);
  sample_indices(indices, fs_hash, size, reduced_size, num_tests);

  for (size_t r = 0; r < num_rounds - 1; r++){
    TRLWE_DFT * enc_points = (TRLWE_DFT *) proof_stream->pull();
    if(r == 0){ // first round: just batch
      decrypt_and_batch_points(previous_points, enc_points, 2*num_tests, beta, sk);
      hash_packed_points(points_hash, enc_points, 2*num_tests);
    }else{
      decrypt_packed_points(curr_points, enc_points, num_tests, sk);
      hash_packed_points(points_hash, enc_points, num_tests);
      fp16_in_u64(alpha, alphas[r - 1]);
      for (size_t i = 0; i < num_tests; i++){
        const uint64_t rnd_idx = indices[i]&(idx_mod_mask>>1);
        if(r == 1){ // points from round 0 are not packed together
          fp16_in_u64(ay, &previous_points[i*D]);
          fp16_in_u64(by, &previous_points[i*D + num_tests*D]);
        }else{
          fp16_in_u64(ay, &previous_points[i*2*D]);
          fp16_in_u64(by, &previous_points[i*2*D + D]);
        }
        fp16_in_u64(cay, &curr_points[i*2*D]);
        fp16_in_u64(cby, &curr_points[i*2*D + D]);

        fp16_in_const64(inv_omega_power, inv_omega_powers[omega_power*rnd_idx]);

        fp16_add(t0, ay, by); // (ay + by)
        fp16_sub(t1, ay, by); // (ay - by)
        fp16_mul(t1, t1, alpha); // (ay - by) * alpha
        fp16_mul(t1, t1, inv_omega_power); // (ay - by) * alpha * inv_omega_power
        fp16_add(t0, t0, t1); // (ay + by) + (ay - by) * alpha * inv_omega_power
        fp16_mul(t0, t0, inv_two); // inv_two * ((ay + by) + (ay - by) * alpha * inv_omega_power)

        if(fp16_cmp(t0, cay) + fp16_cmp(t0, cby) == 0){
          printf("round %lu colinearity check %lu failure\n", r, i);
          return false;
        }
      }
      tmp_pointer = previous_points;
      previous_points = curr_points; 
      curr_points = tmp_pointer;
    }

    for (size_t i = 0; i < num_tests; i++){
      MerklePath path = (MerklePath) proof_stream->pull();
      if(!Merkle::verify(merkle_roots[r], indices[i]&(idx_mod_mask>>1), path, &points_hash[4*i])){
        printf("Merkle authentication failed points a\n");
        printf("Index %lu: %lu\n", i, indices[i]&(idx_mod_mask>>1));
        return false;
      }
    }

    if(r == 0){
      for (size_t i = 0; i < num_tests; i++){
        MerklePath path = (MerklePath) proof_stream->pull();
        if(!Merkle::verify(merkle_roots[r], (indices[i]&(idx_mod_mask>>1)) + (size>>(r+1)), path, &points_hash[4*(num_tests + i)])){
          printf("Merkle authentication failed points b\n");
          return false;
        }
      }
    }else{
      idx_mod_mask >>= 1;
      omega_power <<= 1;
    }
    
  }

  // test last codeword
  uint64_t * last_codeword_dec[reduced_size];
  uint64_t * last_poly[reduced_size];
  decrypt_codeword2(last_codeword_dec, last_codeword, sk);

  assert((idx_mod_mask>>2) == ((reduced_size/2) - 1));
  // colin check for last codeword
  fp16_in_u64(alpha, alphas[num_rounds - 2]);
  for (size_t i = 0; i < num_tests; i++){
    const uint64_t prev_idx = indices[i]&(idx_mod_mask>>1);
    const uint64_t curr_idx = indices[i]&(idx_mod_mask>>2);
    if(num_rounds == 2){ // points from round 0 are not packed together
      fp16_in_u64(ay, &previous_points[i*D]);
      fp16_in_u64(by, &previous_points[i*D + num_tests*D]);
    }else{
      fp16_in_u64(ay, &previous_points[i*2*D]);
      fp16_in_u64(by, &previous_points[i*2*D + D]);
    }
    fp16_in_u64(cay, last_codeword_dec[curr_idx]);
    fp16_in_u64(cby, last_codeword_dec[curr_idx + reduced_size/2]);

    fp16_in_const64(inv_omega_power, inv_omega_powers[omega_power*prev_idx]);

    fp16_add(t0, ay, by); // (ay + by)
    fp16_sub(t1, ay, by); // (ay - by)
    fp16_mul(t1, t1, alpha); // (ay - by) * alpha
    fp16_mul(t1, t1, inv_omega_power); // (ay - by) * alpha * inv_omega_power
    fp16_add(t0, t0, t1); // (ay + by) + (ay - by) * alpha * inv_omega_power
    fp16_mul(t0, t0, inv_two); // inv_two * ((ay + by) + (ay - by) * alpha * inv_omega_power)

    if(fp16_cmp(t0, cay) + fp16_cmp(t0, cby) == 0){
      printf("last colinearity check failure\n");
      printf("Test %lu\n", i);
      return false;
    }
  }

  // degree test last codeword
  
  uint64_t last_inv_omega = inv_omega_powers[1<<(num_rounds -1)];
  clear_ntt(last_poly, last_codeword_dec, last_inv_omega, p, reduced_size, D);
  
  if(!check_degree(last_poly, D, reduced_size/expansion, reduced_size)){
    printf("last codeword does not correspond to polynomial of low enough degree\n");
    return false;
  }
  
  return true;

}