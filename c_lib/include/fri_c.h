#include <mosfhet.h>
// #include <mpLWE.h>
#include <blake3.h>
#include <pthread.h>
#include <rns-rlwe.h>
// #include <debug_util.h>
#pragma once

#define MAX_THREADS 128

#ifdef __cplusplus
extern "C" { 
#endif 

typedef struct _PublicKey{
  uint64_t N, D, p, p_log, decomp_bits;
  TRLWE_KS_Key packing_ksk, rlwe_ksk;
  TLWE * tmp;
  TRLWE tmp_mrlwe;
} * PublicKey;

typedef struct _PrivateKey{
  TRLWE_Key packing_key, input_key, output_key;
  TLWE_Key extracted_key;
  TorusPolynomial tmp;
  PublicKey pk;
  RNS_RLWE_Key rns_rlwe_key;
} * PrivateKey;

typedef struct _Codeword{
  TRLWE_DFT * parts;
  uint64_t size;
} * Codeword;


// encrypt
uint64_t mod_switch_noise(uint64_t v, uint64_t p, uint64_t q);
Codeword new_empty_codeword(uint64_t size, uint64_t N);
Codeword new_empty_repacked_codeword(uint64_t size, uint64_t k, uint64_t N);
Codeword encrypt_codeword(uint64_t * vector, uint64_t size, PrivateKey sk);
RNS_RLWE * encrypt_input(uint64_t * vector, uint64_t input_size, uint64_t expansion, PrivateKey sk);
void decrypt_and_print(RNS_RLWE * c, uint64_t input_size, uint64_t coeff_size, PrivateKey sk);
Codeword repack_init_codeword(RNS_RLWE * first_codeword, uint64_t n);
void decrypt_poly(uint64_t * out, RNS_RLWE * c, uint64_t input_size, uint64_t idx, PrivateKey sk);
PublicKey new_public_key(PrivateKey sk, uint64_t D, uint64_t p, uint64_t decomp_bits, uint64_t t_decomp, uint64_t b_decomp, uint64_t t_dimred);
PublicKey copy_public_key(PublicKey pk);
PrivateKey new_private_key(uint64_t in_N, uint64_t in_digits, uint64_t in_digit_size, double in_sigma_err, double in_sigma_key, uint64_t packing_N, uint64_t packing_k, double packing_sigma);
void mprlwe_to_trlwe(TRLWE_DFT out, RNS_RLWE in);

// homNTT
void ntt(RNS_RLWE * out, RNS_RLWE * in, uint64_t root_of_unity, uint64_t p, uint64_t size, uint64_t ntt_max_depth, uint64_t m);
void ntt_mt(RNS_RLWE * out, RNS_RLWE * in, uint64_t root_of_unity, uint64_t p, uint64_t size, uint64_t ntt_max_depth, uint64_t m);
TLWE * batch_codewords(RNS_RLWE * in, uint64_t size, uint64_t * alphav, uint64_t D, uint64_t num_codewords);
TLWE * batch_codewords_mt(RNS_RLWE * in, uint64_t size, uint64_t * alphav, uint64_t D, uint64_t num_codewords, uint64_t num_threads);
void rlwe_inner_prod(TLWE * out, RNS_RLWE in, uint64_t ** alpha, uint64_t size);


// folding
void repack(TRLWE_DFT out, TLWE * c, uint64_t size, PublicKey params);
void decomp_repack(TRLWE_DFT out, TLWE * c, uint64_t size, PublicKey params);
void recomp(uint64_t * out, uint64_t size, TRLWE_DFT in, PrivateKey sk);
void recomp_and_print(uint64_t size, Codeword in, uint64_t index, PrivateKey sk);
TRLWE_DFT * codeword_get_packed_points(Codeword cw, uint64_t * index, uint64_t size);
uint64_t * decrypt_recomp_packed_points(TRLWE_DFT * points, uint64_t size, PrivateKey sk);
void decrypt_packed_points(uint64_t * out, TRLWE_DFT * points, uint64_t size, PrivateKey sk);
void decrypt_and_batch_points(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t * alphav, PrivateKey sk);
void decrypt_codeword(uint64_t * out, Codeword c, PrivateKey sk);
void hash_packed_points(uint64_t * out, TRLWE_DFT * points, uint64_t size);
void hash_codeword(uint64_t * out, Codeword c);
void hash_init_codewords(uint64_t * out, RNS_RLWE * points, uint64_t size);
void codeword_get_point(TLWE out, Codeword cw, uint64_t index);
void * repack_codeword(Codeword first_codeword, uint64_t n, PublicKey params);
void * repack_codeword_mt(Codeword first_codeword, uint64_t n, PublicKey params, uint64_t num_threads);
void * folding_level(Codeword first_codeword, uint64_t ** consts_array, uint64_t n, uint64_t r, PublicKey params);
void * folding_level_LWE(TLWE * first_codeword, uint64_t ** consts_array, uint64_t n, uint64_t r, PublicKey params);
void * folding_level_LWE_mt(TLWE * first_codeword, uint64_t ** consts_array, uint64_t n, uint64_t r, PublicKey params, uint64_t num_threads);
void decrypt_codeword2(uint64_t ** out, Codeword c, PrivateKey sk);
// void decrypt_and_batch_points2(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t * alphav, PrivateKey sk);
uint64_t * hash_full_codeword(Codeword c);


// debug functions
void recomp_and_print(uint64_t size, Codeword in, uint64_t index, PrivateKey sk);

// field functions
typedef uint32_t fp_t;
typedef fp_t fp2_t[2];
typedef fp2_t fp4_t[2];
typedef fp4_t fp8_t[2];
typedef fp8_t fp16_t[2];

typedef TLWE fp_he_t;
typedef fp_he_t fp2_he_t[2];
typedef fp2_he_t fp4_he_t[2];
typedef fp4_he_t fp8_he_t[2];
typedef fp8_he_t fp16_he_t[2];

void fp_ini(fp_t prime);
void fp_setup_temps(uint64_t n);
void fp16_he_in(fp16_he_t a, fp_he_t in[16]);
void fp16_in(fp16_t a, fp_t in[16]);
void fp16_in_u64(fp16_t a, uint64_t in[16]);
void fp16_mul_he(fp16_he_t c, const fp16_t a, const fp16_he_t b);
void fp16_add_he(fp16_he_t c, const fp16_he_t a, const fp16_he_t b);
void fp16_print(fp16_t a);
void fp16_he_print(fp16_he_t a, TLWE_Key key, uint64_t p);
void fp16_mul(fp16_t c, const fp16_t a, const fp16_t b);
void fp16_get_powers(fp16_t * out, fp16_t in, uint32_t size);
void fp16_sqr(fp16_t c, const fp16_t a);
void fp16_add(fp16_t c, const fp16_t a, const fp16_t b);
void fp16_sub(fp16_t c, const fp16_t a, const fp16_t b);
void fp16_rnd(fp16_t a);
void fp16_inv(fp16_t c, const fp16_t a);
int fp16_cmp(const fp16_t a, const fp16_t b);
void fp16_out(fp_t out[16], fp16_t a);
void fp16_in_const64(fp16_t a, uint64_t in);

// c verifier functions
uint64_t * sample_field_array(uint64_t * hash_in, uint64_t size, uint64_t p);
void sample_field_array2(uint64_t * res, uint64_t * hash_in, uint64_t size, uint64_t p);
void clear_ntt(uint64_t ** out, uint64_t ** in, uint64_t w, uint64_t p, uint64_t n, uint64_t D);
bool verify(void * proof_stream, uint64_t size, uint64_t num_tests, uint64_t num_rounds, uint64_t omega, uint64_t expansion, PrivateKey sk);
void * load_proof_stream(void ** data, uint64_t * list_hashes, uint64_t size);
bool verify_benchmark(void * proof_stream_p, uint64_t size, uint64_t num_tests, uint64_t num_rounds, uint64_t omega, uint64_t expansion, PrivateKey sk);

void sample_indices(uint64_t * out, uint64_t * seed, uint64_t size, uint64_t reduced_size, uint64_t number);
void decrypt_and_batch_points_mt(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t * alphav, PrivateKey sk, uint64_t num_threads);
void hash_packed_points_mt(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t num_threads);
void hash_codeword_mt(uint64_t * out, Codeword c, uint64_t num_threads);

#ifdef __cplusplus
}
#endif