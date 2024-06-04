#include <fri_c.h>

void hash_packed_points(uint64_t * out, TRLWE_DFT * points, uint64_t size){
  // Initialize the hasher.
  blake3_hasher hasher;
  for (size_t i = 0; i < size; i++){
    blake3_hasher_init(&hasher);
    assert(points[i]->b->N > 0);
    for (size_t j = 0; j < points[i]->k; j++){
      blake3_hasher_update(&hasher, points[i]->a[j]->coeffs, points[i]->a[j]->N*8);
    }
    blake3_hasher_update(&hasher, points[i]->b->coeffs, points[i]->b->N*8);
    blake3_hasher_finalize(&hasher, (uint8_t *) &out[4*i], BLAKE3_OUT_LEN);
  }
}

typedef struct _hash_arg{
  uint64_t * out, start, end;
  TRLWE_DFT * points;
} hash_arg;

void * hash_packed_points_thread(void * arg_p){
  hash_arg * arg = (hash_arg *) arg_p;
  blake3_hasher hasher;
  for (size_t i = arg->start; i < arg->end; i++){
    blake3_hasher_init(&hasher);
    for (size_t j = 0; j < arg->points[i]->k; j++){
      blake3_hasher_update(&hasher, arg->points[i]->a[j]->coeffs, arg->points[i]->a[j]->N*8);
    }
    blake3_hasher_update(&hasher, arg->points[i]->b->coeffs, arg->points[i]->b->N*8);
    blake3_hasher_finalize(&hasher, (uint8_t *) &arg->out[4*i], BLAKE3_OUT_LEN);
  }
  return NULL;
}

void hash_packed_points_mt(uint64_t * out, TRLWE_DFT * points, uint64_t size, uint64_t num_threads){
  pthread_t threads[MAX_THREADS];
  hash_arg args[MAX_THREADS];
  // Initialize the hasher.
  const uint64_t n_points = size/num_threads;
  for (size_t i = 0; i < num_threads; i++){
    args[i].out = out; 
    args[i].points = points; 
    args[i].start = i*n_points; 
    args[i].end = (i+1)*n_points;
    if(i == num_threads-1) args[i].end = size;
    pthread_create(&threads[i], NULL, *hash_packed_points_thread, (void *) &(args[i]));
  }
  for (size_t i = 0; i < num_threads; i++) pthread_join(threads[i], NULL);
}

void hash_codeword_mt(uint64_t * out, Codeword c, uint64_t num_threads){
  hash_packed_points_mt(out, c->parts, c->size, num_threads);
}

void hash_codeword(uint64_t * out, Codeword c){
  // Initialize the hasher.
  blake3_hasher hasher;
  for (size_t i = 0; i < c->size; i++){
    blake3_hasher_init(&hasher);
    for (size_t j = 0; j < c->parts[i]->k; j++){
      blake3_hasher_update(&hasher, c->parts[i]->a[j]->coeffs, c->parts[i]->a[j]->N*sizeof(uint64_t));
    }
    blake3_hasher_update(&hasher, c->parts[i]->b->coeffs, c->parts[i]->b->N*sizeof(uint64_t));
    blake3_hasher_finalize(&hasher, (uint8_t *) &out[4*i], BLAKE3_OUT_LEN);
  }
}

// a single hash for the entire codeword
uint64_t * hash_full_codeword(Codeword c){
  uint64_t * out = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*4);
  // Initialize the hasher.
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  for (size_t i = 0; i < c->size; i++){
    for (size_t j = 0; j < c->parts[i]->k; j++){
      blake3_hasher_update(&hasher, c->parts[i]->a[j]->coeffs, c->parts[i]->a[j]->N*sizeof(uint64_t));
    }
    blake3_hasher_update(&hasher, c->parts[i]->b->coeffs, c->parts[i]->b->N*sizeof(uint64_t));
  }
  blake3_hasher_finalize(&hasher, (uint8_t *) out, BLAKE3_OUT_LEN);
  return out;
}

void hash_init_codewords(uint64_t * out, RNS_RLWE * points, uint64_t size){
  // Initialize the hasher.
  blake3_hasher hasher;
  for (size_t i = 0; i < size; i++){
    blake3_hasher_init(&hasher);
    for (size_t j = 0; j < points[i]->a->l; j++){
      blake3_hasher_update(&hasher, points[i]->a->coeffs[j], points[i]->a->N*sizeof(uint64_t));
      blake3_hasher_update(&hasher, points[i]->b->coeffs[j], points[i]->b->N*sizeof(uint64_t));
    }
    blake3_hasher_finalize(&hasher, (uint8_t *) &out[4*i], BLAKE3_OUT_LEN);
  }
}
