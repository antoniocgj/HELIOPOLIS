#include "rns-rlwe.h"

uint64_t next_power_of_2(uint64_t x){
  return 2ULL << intel::hexl::Log2(x);
}

uint64_t double2int(double x){
  return ((uint64_t) ((int64_t) x));
}

uint64_t mod_switch(uint64_t v, uint64_t p, uint64_t q) {
  const double double_q = q == 0? pow(2,64) : ((double)q); 
  const double double_p = p == 0? pow(2,64) : ((double)p); 
  uint64_t val = (uint64_t) round((((double) v) * double_q) / double_p);
  return val < q? val : val - q;
}

// Generates a sparse ternary array with Hamming Weight h, balanced (h/2 ones and h/2 negative ones)
void gen_sparse_ternary_array_modq(uint64_t * out, uint64_t size, uint64_t h, uint64_t q){
  memset(out, 0, sizeof(uint64_t)*size);
  uint64_t hw = 0, val = 1, * rnd_buffer;
  const uint64_t buffer_size = h*10;
  rnd_buffer = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*buffer_size);
  while(hw < h){
    generate_random_bytes(sizeof(uint64_t)*buffer_size, (uint8_t *) rnd_buffer); 
    array_mod_switch_from_2k(rnd_buffer, rnd_buffer, size, size, buffer_size);
    uint64_t i = 0; 
    while (i < buffer_size && hw < h){
      const uint64_t idx = rnd_buffer[i++];
      if(out[idx]) continue;
      out[idx] = (q + val)%q;
      val *= -1;
      hw++;
    }
  }
  free(rnd_buffer);
  #ifndef NDEBUG
  uint64_t hw_check = 0, sum_check = 0;
  for (size_t i = 0; i < size; i++){
    sum_check += out[i];
    hw_check += (out[i] != 0);
  }
  assert(hw_check == h);
  assert((sum_check%q) == 0);   
  #endif
}

void array_reduce_mod_N(uint64_t * out,  uint64_t * in, uint64_t size, uint64_t p){
  const uint64_t mask = next_power_of_2(p) - 1;
  for (size_t i = 0; i < size; i++){
    out[i] = in[i]&mask;
  } 
}

/* Mod switch the additive inverse of negative values */
/* Used to adjust negative values in Gaussian keys when represented by the inverse */
void array_additive_inverse_mod_switch(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n){
  for (size_t i = 0; i < n; i++){
    if(in[i] > p/2) out[i] = q - (p - in[i]);
    else out[i] = in[i];
  }
}

/* Switch each element mod p to mod q*/
void array_mod_switch(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n){
  for (size_t i = 0; i < n; i++){
    // out[i] = (round((((double) in[i])/((double) p))*q));
    out[i] = mod_switch(in[i], q, p);
  }
}

/* Switch each element mod next_power_of_two(p) to mod q */
void array_mod_switch_from_2k(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n){
  array_reduce_mod_N(out, in, n, p);
  p = next_power_of_2(p);
  for (size_t i = 0; i < n; i++){
    // out[i] = (round((((double) in[i])/((double) p))*q));
    out[i] = mod_switch(in[i], q, p);
  }
}

uint64_t int_mod_switch(uint64_t in, uint64_t p, uint64_t q){
  // return (round((((double) in)/((double) p))*q));
  return mod_switch(in, q, p);
}

intel::hexl::NTT ** new_ntt_list(uint64_t * primes, uint64_t N, uint64_t l){
  intel::hexl::NTT ** ntt = (intel::hexl::NTT **) safe_malloc(sizeof(intel::hexl::NTT *)*l);
  for (size_t i = 0; i < l; i++){
    ntt[i] = new intel::hexl::NTT(N, primes[i]);
  }
  return ntt;
}

// Computes (Z_q[i](Q/q[i]))**-1, for i in [0,l)
void compute_RNS_Qhat_array(uint64_t * out, uint64_t * p, uint64_t l){
  for (size_t i = 0; i < l; i++){
    out[i] = 1;
    for (size_t j = 0; j < l; j++){
      if(i!=j){
        const uint64_t inv = intel::hexl::InverseMod(p[j], p[i]);
        out[i] = intel::hexl::MultiplyMod(out[i], inv, p[i]);
      }
    }
  }
}

// out = a*b
// a and out are (size*64)-bit multi-precision integers 
// b is a 64-bit integer
void mpi_mul_int(uint64_t * out, uint64_t * a, uint64_t b, uint64_t size){
  uint64_t carry, tmp = a[0];
  out[0] = 0;
  for (size_t i = 0; i < size - 1; i++){
    out[i] += _mulx_u64 (tmp, b, (unsigned long long *)&carry);
    tmp = out[i + 1];
    out[i + 1] = carry;
  }
}

void mpi_addto(uint64_t * out, uint64_t * in, uint64_t size){

}

void RNS_compose(uint64_t * out, uint64_t * in, uint64_t * p, uint64_t * qi_hat, uint64_t l){
  // WIP
  // Z_Q(Q/q[i]) * Z_Q((Z_q[i](Q/q[i]))**-1)
  assert(l <= 32);
  uint64_t prod[32];
  for (size_t i = 0; i < l; i++){
    memset(prod, 0, sizeof(uint64_t)*l);
    prod[0] = in[0];
    // prod Q/q_i
    for (size_t j = 0; j < l; j++){
      if(i!=j){
        mpi_mul_int(prod, prod, p[j], l);
      }
    }
    // prod Qi_hat
    mpi_mul_int(prod, prod, qi_hat[i], l);
    // add to out

  }
}

double __debug_global_avg_var = 0.;
uint64_t __debug_decryption_count = 0;
// Returns m from a noisy RNS representation of Q_0*X^m, where Q_0 = Q/q_0
uint64_t __debug_get_exp_message_from_noisy_RNS(RNS_Polynomial in, uint64_t * p){
  RNSc_Polynomial tmp = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->N, in->l, in->ntt);
  polynomial_RNS_to_RNSc(tmp, in);
  for (size_t i = 0; i < in->l - 1 ; i++){
    polynomial_base_reduce_RNSc(tmp);
  }
  uint64_t res = -1, max = 0;
  double var = 0;
  for (size_t i = 0; i < in->N; i++){
    const uint64_t val = tmp->coeffs[0][i] > p[0]/2 ? p[0] - tmp->coeffs[0][i] : tmp->coeffs[0][i];
    // printf("%ld ", val);
    if(val > max){
      res = i;
      var += max*max;
      max = val;
    }else{
      var += val*val;
    }
  }
  var /= in->N;
  __debug_global_avg_var = __debug_global_avg_var*__debug_decryption_count + var;
  __debug_decryption_count += 1;
  __debug_global_avg_var /= __debug_decryption_count;
  // printf("\n");
  free_RNS_polynomial(tmp);
  return res;
}

uint64_t __debug_get_exp_message_from_noisy_RNSc(RNSc_Polynomial in, uint64_t * p){
  RNSc_Polynomial tmp = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->N, in->l, in->ntt);
  polynomial_copy_RNS_polynomial((RNS_Polynomial) tmp, (RNS_Polynomial) in);
  for (size_t i = 0; i < in->l - 1 ; i++){
    polynomial_base_reduce_RNSc(tmp);
  }
  uint64_t res = -1, max = 0;
  for (size_t i = 0; i < in->N; i++){
    const uint64_t val = tmp->coeffs[0][i] > p[0]/2 ? p[0] - tmp->coeffs[0][i]  : tmp->coeffs[0][i];
    if(val > max){
      res = i;
      max = val;
    }
  }
  free_RNS_polynomial(tmp);
  return res;
}

uint64_t mod_dist(uint64_t a, uint64_t b, uint64_t q){
  const uint64_t dist = (q + a-b)%q;
  if(dist > q/2) return q - dist;
  return dist;
}

void print_array(const char * msg, uint64_t * v, size_t size){
  printf("%s: ", msg);
  for (size_t i = 0; i < size; i++){
    printf("%ld, ", v[i]);
  }
  printf("\n");
}