#include <iop.hpp>

void Merkle::hash256(uint64_t out[4], uint64_t in[4]){
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  blake3_hasher_update(&hasher, in, 32);
  blake3_hasher_finalize(&hasher, (uint8_t *) out, BLAKE3_OUT_LEN);
}

void Merkle::hash256_2(uint64_t out[4], uint64_t in1[4], uint64_t in2[4]){
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  blake3_hasher_update(&hasher, in1, 32);
  blake3_hasher_update(&hasher, in2, 32);
  blake3_hasher_finalize(&hasher, (uint8_t *) out, BLAKE3_OUT_LEN);
}

void Merkle::commit(uint64_t * hash_array){
  memcpy(tree[log_size], hash_array, 32*size);
  
  for (int32_t i = ((int32_t) log_size) - 1; i >= 0; i--){
    for (size_t j = 0; j < (1ULL<<i); j++){
      hash256_2(&tree[i][hash_wsize*j], &tree[i+1][hash_wsize*2*j], &tree[i+1][hash_wsize*(2*j + 1)]); 
    }
  }
}

void print_vec(char * msg, uint64_t * vec, uint64_t size){
  printf("%s: [", msg);
  for (size_t i = 0; i < size; i++){
    printf("%lu, ", vec[i]);
  }
  printf("]\n");
}

void Merkle::open(MerklePath path, uint64_t index){
  // printf("Index: %lu\n", index);
  for (size_t i = log_size; i > 0; i--){
    // print_vec("\n\n\n", tree[i], 4*(1ULL<<i));
    memcpy(&path->hash_list[hash_wsize*(log_size - i)], &tree[i][hash_wsize*(index^1)], sizeof(uint64_t)*hash_wsize); 
    index>>=1;
  }
  // print_vec("\n\n\n", tree[0], 4);
  // print_vec("path: ", path->hash_list, 4*log_size);
}

MerklePath Merkle::open(uint64_t index){
  MerklePath path = (MerklePath) safe_malloc(sizeof(MerklePath));
  path->hash_list = (uint64_t *) safe_aligned_malloc(log_size*hash_wsize*sizeof(uint64_t));
  path->size = log_size;
  open(path, index);
  return path;
}

uint64_t * Merkle::get_root(){
  uint64_t * root = (uint64_t *) safe_aligned_malloc(hash_wsize*sizeof(uint64_t));
  memcpy(root, tree[0], sizeof(uint64_t)*hash_wsize);
  return root;
}

bool Merkle::verify(uint64_t * root, uint64_t index, MerklePath path, uint64_t * element_hash){
  // printf("Index: %lu\n", index);
  // print_vec("path: ", path->hash_list, 4*path->size);
  // print_vec("root: ", root, 4);

  const uint64_t hash_wsize = Merkle::hash_wsize;
  uint64_t hash[hash_wsize];
  if(index&1) Merkle::hash256_2(hash, path->hash_list, element_hash);
  else Merkle::hash256_2(hash, element_hash, path->hash_list);
  index >>= 1;

  for (size_t i = 1; i < path->size; i++){
    if(index&1) Merkle::hash256_2(hash, &path->hash_list[hash_wsize*i], hash);
    else Merkle::hash256_2(hash, hash, &path->hash_list[hash_wsize*i]);
    index >>= 1;
  }
  for (size_t i = 0; i < hash_wsize; i++){
    if(hash[i] != root[i]) return false;
  }
  return true;
}

Merkle::Merkle(uint64_t size){
  this->size = size;
  this->log_size = ceil(log2(size));
  // alloc tree
  this->tree = (uint64_t **) safe_malloc((log_size+1)*sizeof(uint64_t *));
  for (size_t i = 0; i <= log_size; i++){
    this->tree[i] = (uint64_t *) safe_aligned_malloc((1ULL<<i)*sizeof(uint64_t)*hash_wsize);
  }
}

Merkle::~Merkle()
{
  for (size_t i = 0; i <= log_size; i++){
    free(this->tree[i]);
  }
  free(this->tree);
}


// c bindings
uint64_t * merkle_get_root(uint64_t * hash_array, uint64_t size){
  Merkle mk(size);
  mk.commit(hash_array);
  return mk.get_root();
}

MerklePath merkle_open(uint64_t * hash_array, uint64_t index, uint64_t size){
  // printf("open: %p %lu %lu\n", hash_array, index, size);
  Merkle mk(size);
  mk.commit(hash_array);
  return mk.open(index);
}

bool merkle_verify(uint64_t * root, uint64_t index, MerklePath path, uint64_t * element_hash){
  return Merkle::verify(root, index, path, element_hash);
}