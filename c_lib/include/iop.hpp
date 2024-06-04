#include <fri_c.h>
#pragma once

class ProofStream
{
private:
  void ** stream;
  uint64_t * stream_hashes;
  uint64_t max_size, curr_size, read_idx;
  static const uint64_t hash_wsize = 4;
  void hash(uint64_t * out, uint64_t * in, uint64_t in_size);
public:
  ProofStream(uint64_t max_size);
  ~ProofStream();
  void push(void * obj, uint64_t * hash);
  void * pull();
  void fiat_shamir(uint64_t * out);
  void read_fiat_shamir(uint64_t * out);
  void reset_read();
  void load(void ** list, uint64_t * list_hashes, uint64_t size);
};

typedef struct _MerklePath{
  uint64_t * hash_list;
  uint64_t size;
} * MerklePath;


class Merkle
{
private:
  uint64_t ** tree, size, log_size;
  static const uint64_t hash_wsize = 4;
  static void hash256(uint64_t out[4], uint64_t in[4]);
  static void hash256_2(uint64_t out[4], uint64_t in1[4], uint64_t in2[4]);
public:
  Merkle(uint64_t size);
  ~Merkle();
  void commit(uint64_t * hash_array);
  void open(MerklePath path, uint64_t index);
  MerklePath open(uint64_t index);
  uint64_t * get_root();
  static bool verify(uint64_t * root, uint64_t index, MerklePath path, uint64_t * element_hash);
};

#ifdef __cplusplus
extern "C" { 
#endif
// merkle 
uint64_t * merkle_get_root(uint64_t * hash_array, uint64_t size);
MerklePath merkle_open(uint64_t * hash_array, uint64_t index, uint64_t size);
bool merkle_verify(uint64_t * root, uint64_t index, MerklePath path, uint64_t * element_hash);

// ip
void * proof_stream_setup(uint64_t max_size);
void proof_stream_push(void * proof_stream, void * obj, uint64_t * hash);
void * proof_stream_pull(void * proof_stream);
uint64_t * proof_stream_fs(void * proof_stream);
uint64_t * proof_stream_read_fs(void * proof_stream);
#ifdef __cplusplus
}
#endif