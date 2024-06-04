#include <iop.hpp>

void ProofStream::hash(uint64_t * out, uint64_t * in, uint64_t in_size){
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  blake3_hasher_update(&hasher, in, in_size);
  blake3_hasher_finalize(&hasher, (uint8_t *) out, BLAKE3_OUT_LEN);
}

ProofStream::ProofStream(uint64_t max_size){
  this->max_size = max_size;
  stream = (void **) safe_malloc(sizeof(void *)*max_size);
  stream_hashes = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*max_size*hash_wsize);
  read_idx = 0;
  curr_size = 0;
}

void ProofStream::reset_read(){
  read_idx = 0;
}

ProofStream::~ProofStream()
{
}

void ProofStream::push(void * obj, uint64_t * hash){
  assert(curr_size < max_size);
  stream[curr_size] = obj;
  if(hash){
    memcpy(&stream_hashes[hash_wsize*curr_size], hash, hash_wsize*sizeof(uint64_t));
  }else{
    memset(&stream_hashes[hash_wsize*curr_size], 0, hash_wsize*sizeof(uint64_t));
  }
  curr_size++;
}

void * ProofStream::pull(){
  assert(read_idx < curr_size);
  return stream[read_idx++];
}

void ProofStream::fiat_shamir(uint64_t * out){
  hash(out, stream_hashes, hash_wsize*curr_size*sizeof(uint64_t));
}

void ProofStream::read_fiat_shamir(uint64_t * out){
  assert(read_idx > 0);
  hash(out, stream_hashes, hash_wsize*read_idx*sizeof(uint64_t));
}

void ProofStream::load(void ** list, uint64_t * list_hashes, uint64_t size){
  memcpy(stream, list, sizeof(void *)*size);
  memcpy(stream_hashes, list_hashes, sizeof(uint64_t)*hash_wsize*size);
  curr_size = size;
}


// C bindings 
void * load_proof_stream(void ** data, uint64_t * list_hashes, uint64_t size){
  ProofStream * res = new ProofStream(size);
  res->load(data, list_hashes, size);
  return res;
}


void * proof_stream_setup(uint64_t max_size){
  ProofStream * ps = new ProofStream(max_size);
  return ps;
}

void proof_stream_push(void * proof_stream, void * obj, uint64_t * hash){
  ProofStream * ps = (ProofStream *) proof_stream;
  ps->push(obj, hash);
}

void * proof_stream_pull(void * proof_stream){
  ProofStream * ps = (ProofStream *) proof_stream;
  return ps->pull();
}

uint64_t * proof_stream_fs(void * proof_stream){
  ProofStream * ps = (ProofStream *) proof_stream;
  uint64_t * res =  (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*4);
  ps->fiat_shamir(res);
  return res;
}

uint64_t * proof_stream_read_fs(void * proof_stream){
  ProofStream * ps = (ProofStream *) proof_stream;
  uint64_t * res =  (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*4);
  ps->read_fiat_shamir(res);
  return res;
}


