import ctypes, os
from ctypes import c_uint64, c_void_p, c_double, POINTER, c_bool
from math import floor, log2
import subprocess

BENCH = False

class FriC:
  def __init__(self) -> None:
    self.multithreaded = False
    self.num_threads = 1
    if not os.path.exists("c_lib/lib/libfric.so"):
      print("Compiling HEXL\n\n")
      os.system("make -C c_lib hexl")
      print("\n\n\n\nCompiling FRI C library")
      output = subprocess.check_output("cat /proc/cpuinfo", shell=True)
      if("avx512ifma" in str(output)):
        os.system("make -C c_lib")
      else:
        print("AVX512ifma not detected - Compiling unoptimized version")
        os.system("make -C c_lib A_PRGN=none ENABLE_VAES=false FFT_LIB=spqlios")
    self.lib = ctypes.CDLL('c_lib/lib/libfric.so')
    # encrypt codeword
    self.lib.encrypt_codeword.argtypes = (ctypes.POINTER(c_uint64), c_uint64, c_void_p)
    self.lib.encrypt_codeword.restype = c_void_p
    # encrypt input
    self.lib.encrypt_input.argtypes = (ctypes.POINTER(c_uint64), c_uint64, c_uint64, c_void_p)
    self.lib.encrypt_input.restype = c_void_p
    # get points
    self.lib.codeword_get_packed_points.argtypes = (c_void_p, ctypes.POINTER(c_uint64), c_uint64)
    self.lib.codeword_get_packed_points.restype = c_void_p
    # decrypt points
    self.lib.decrypt_packed_points.argtypes = (POINTER(c_uint64), c_void_p, c_uint64, c_void_p)
    # hash points
    self.lib.hash_packed_points.argtypes = (ctypes.POINTER(c_uint64), c_void_p, c_uint64)
    # hash init codewords
    self.lib.hash_init_codewords.argtypes = (ctypes.POINTER(c_uint64), c_void_p, c_uint64)
    # hash codeword
    self.lib.hash_codeword.argtypes = (ctypes.POINTER(c_uint64), c_void_p)
    self.lib.hash_codeword_mt.argtypes = (ctypes.POINTER(c_uint64), c_void_p, c_uint64)
    # hash full codeword
    self.lib.hash_full_codeword.argtypes = (c_void_p, )
    self.lib.hash_full_codeword.restype = ctypes.POINTER(c_uint64)
    # folding 
    self.lib.folding_level.argtypes = (c_void_p, ctypes.POINTER(ctypes.POINTER(c_uint64)), c_uint64, c_uint64, c_void_p)
    self.lib.folding_level.restype = c_void_p
    # folding LWE input (batched version)
    self.lib.folding_level_LWE.argtypes = (c_void_p, ctypes.POINTER(ctypes.POINTER(c_uint64)), c_uint64, c_uint64, c_void_p)
    self.lib.folding_level_LWE.restype = c_void_p
    # folding LWE input (batched version) MT
    self.lib.folding_level_LWE_mt.argtypes = (c_void_p, ctypes.POINTER(ctypes.POINTER(c_uint64)), c_uint64, c_uint64, c_void_p, c_uint64)
    self.lib.folding_level_LWE_mt.restype = c_void_p
    # repack codeword 
    self.lib.repack_codeword.argtypes = (c_void_p, c_uint64, c_void_p)
    self.lib.repack_codeword.restype = c_void_p
    # repack init codeword (batched)
    self.lib.repack_init_codeword.argtypes = (c_void_p, c_uint64)
    self.lib.repack_init_codeword.restype = c_void_p
    # repack codeword MT
    self.lib.repack_codeword_mt.argtypes = (c_void_p, c_uint64, c_void_p, c_uint64)
    self.lib.repack_codeword_mt.restype = c_void_p
    # new sk
    self.lib.new_private_key.argtypes = (c_uint64, c_uint64, c_uint64, c_double, c_double, c_uint64, c_uint64, c_double)
    self.lib.new_private_key.restype = c_void_p
    # new pk
    self.lib.new_public_key.argtypes = (c_void_p, c_uint64, c_uint64, c_uint64, c_uint64, c_uint64, c_uint64)
    self.lib.new_public_key.restype = c_void_p

    # recomp_and_print
    self.lib.recomp_and_print.argtypes = (c_uint64, c_void_p, c_uint64, c_void_p)

    # NTT
    self.lib.ntt.argtypes = (c_void_p, c_void_p, c_uint64, c_uint64, c_uint64, c_uint64, c_uint64)
    self.lib.ntt_mt.argtypes = (c_void_p, c_void_p, c_uint64, c_uint64, c_uint64, c_uint64, c_uint64)

    # decrypt and print
    self.lib.decrypt_and_print.argtypes = (c_void_p, c_uint64, c_uint64, c_void_p)

    # decrypt poly
    self.lib.decrypt_poly.argtypes = (POINTER(c_uint64), c_void_p, c_uint64, c_uint64, c_void_p)

    # decrypt poly
    self.lib.decrypt_codeword.argtypes = (POINTER(c_uint64), c_void_p, c_void_p)

    # batch codewords
    self.lib.batch_codewords.argtypes = (c_void_p, c_uint64, c_void_p, c_uint64, c_uint64)
    self.lib.batch_codewords.restype = c_void_p

    # batch codewords multi-threaded
    self.lib.batch_codewords_mt.argtypes = (c_void_p, c_uint64, c_void_p, c_uint64, c_uint64, c_uint64)
    self.lib.batch_codewords_mt.restype = c_void_p

    # decrypt and batch codewords
    self.lib.decrypt_and_batch_points.argtypes = (POINTER(c_uint64), c_void_p, c_uint64, c_void_p, c_void_p)

    # sampling from hash
    self.lib.sample_field_array.argtypes = (POINTER(c_uint64), c_uint64, c_uint64)
    self.lib.sample_field_array.restype = c_void_p

    # load proof stream
    self.lib.load_proof_stream.argtypes = (POINTER(c_void_p), POINTER(c_uint64), c_uint64)
    self.lib.load_proof_stream.restype = c_void_p

    # sample indices
    self.lib.sample_indices.argtypes = (POINTER(c_uint64), POINTER(c_uint64), c_uint64, c_uint64, c_uint64)

    ### merkle functions

    # get root
    self.lib.merkle_get_root.argtypes = (POINTER(c_uint64), c_uint64)
    self.lib.merkle_get_root.restype = POINTER(c_uint64*4)

    # open
    self.lib.merkle_open.argtypes = (POINTER(c_uint64), c_uint64, c_uint64)
    self.lib.merkle_open.restype = POINTER(c_void_p)

    # verify
    self.lib.merkle_verify.argtypes = (POINTER(c_uint64), c_uint64, c_void_p, POINTER(c_uint64))
    self.lib.merkle_verify.restype = c_bool

    ### IP functions
    self.lib.proof_stream_setup.argtypes = (c_uint64,)
    self.lib.proof_stream_setup.restype = c_void_p

    self.lib.proof_stream_push.argtypes = (c_void_p, c_void_p,  POINTER(c_uint64))

    self.lib.proof_stream_pull.argtypes = (c_void_p,)
    self.lib.proof_stream_pull.restype = c_void_p

    self.lib.proof_stream_fs.argtypes = (c_void_p,)
    self.lib.proof_stream_fs.restype = POINTER(c_uint64)

    self.lib.proof_stream_read_fs.argtypes = (c_void_p,)
    self.lib.proof_stream_read_fs.restype = POINTER(c_uint64)

    self.lib.verify.argtypes = (c_void_p, c_uint64, c_uint64, c_uint64, c_uint64, c_uint64, c_void_p)
    self.lib.verify.restype = c_bool

    self.lib.verify_benchmark.argtypes = (c_void_p, c_uint64, c_uint64, c_uint64, c_uint64, c_uint64, c_void_p)
    self.lib.verify_benchmark.restype = c_bool

  def verify_native(self, proof_stream, size, num_tests, num_rounds, omega, expansion, sk):
    if(BENCH): 
      return self.lib.verify_benchmark(proof_stream, size, num_tests, num_rounds, omega, expansion, sk)
    else:
      return self.lib.verify(proof_stream, size, num_tests, num_rounds, omega, expansion, sk)


  def sample_indices(self, seed, size, reduced_size, count):
    res = (c_uint64*count)()
    self.lib.sample_indices(res, seed, size, reduced_size, count)
    return list(res)

  def merkle_commit(self, hash_array):
    return self.lib.merkle_get_root(hash_array, len(hash_array)//4)

  def merkle_open(self, hash_array, index):
    return self.lib.merkle_open(hash_array, index, len(hash_array)//4)
  
  def merkle_verify(self, root, index, path, element_hash):
    root_p = ctypes.cast(root, ctypes.POINTER(c_uint64))
    element_hash = (c_uint64*4)(*element_hash)
    return self.lib.merkle_verify(root_p, index, path, element_hash)


  def gen_private_key(self, N=4096, size=0):
    if(N == 4096):
      security = {53: 127, 54: 125, 55: 122, 56: 120, 57: 118, 58: 116, 59: 114}
      size_2_q = {2048: 55, 4096: 55, 8192: 56, 16384: 57, 32768: 58, 65536: 59}
      return self.lib.new_private_key(4096, 2, 59, 1, 1, 512, 4, c_double(2**-52))
    elif(N == 8192):
      return self.lib.new_private_key(8192, 3, 49, 1, 1, 512, 4, c_double(2**-52))
    else:
      assert(False) # We only have parameters for N in {4096, 8192}
  
  def gen_public_key(self, field, sk):
    decomp_bits = 2
    p = field.p
    D = field.D
    return self.lib.new_public_key(sk, D, p, decomp_bits, 2, 14, 28)
  
  def encrypt_input(self, input, expansion_factor, sk):
    num_numbers = len(input)
    assert(num_numbers)
    array_type = ctypes.c_uint64 * num_numbers
    return self.lib.encrypt_input(array_type(*input), num_numbers, expansion_factor, sk)
  
  def shallow_ntt(self, input, field, initial_codeword_length):
    p = field.p
    omega = int(field.primitive_nth_root(initial_codeword_length))
    
    next_power_of_2 = lambda x: 1<<int(log2(x))
    
    ntt_max_depth = 2
    m = next_power_of_2(floor((initial_codeword_length/(ntt_max_depth - 1))**(1/((ntt_max_depth - 1) + 1))))
    if(self.multithreaded):
      for i in range(3): 
        if(m%2 == 0 and m > 32): m = m//2
      self.lib.ntt_mt(input, input, omega, p, initial_codeword_length, 2, m)
    else:
      # ntt_max_depth = 3
      m = next_power_of_2(floor((initial_codeword_length/(ntt_max_depth - 1))**(1/((ntt_max_depth - 1) + 1))))
      self.lib.ntt(input, input, omega, p, initial_codeword_length, ntt_max_depth, m)

  def decrypt_poly(self, enc_codeword_array, size, idx, sk, field):
    p = (c_uint64*size)()
    self.lib.decrypt_poly(p, enc_codeword_array, size, idx, sk)
    return [field(i) for i in p]
  
  def hash_packed_points(self, points, size, return_pointer=False):
    p = (c_uint64*(size*4))()
    self.lib.hash_packed_points(p, points, size)
    if(return_pointer): return p
    return [p[4*i:4*(i+1)] for i in range(len(p)//4)]

  def hash_init_codewords(self, points, size, return_pointer=False):
    p = (c_uint64*(size*4))()
    self.lib.hash_init_codewords(p, points, size)
    if(return_pointer): return p
    return [p[4*i:4*(i+1)] for i in range(len(p)//4)]
  
  def hash_codeword(self, points, size, return_pointer=False):
    p = (c_uint64*(size*4))()
    if(self.multithreaded and size > 4):
      num_threads = min(self.num_threads, size//2)
      self.lib.hash_codeword_mt(p, points, num_threads)
    else:
      self.lib.hash_codeword(p, points)
    if(return_pointer): return p
    return [p[4*i:4*(i+1)] for i in range(len(p)//4)]
  
  def hash_full_codeword(self, points):
    return self.lib.hash_full_codeword(points)
 
  
  def decrypt_codeword(self, codeword, size, field, sk):
    D = field.D
    p = (c_uint64*(size*D))()
    self.lib.decrypt_codeword(p, codeword, sk)
    p = list(p)
    poly = [field(p[i:i + D]) for i in range(0,len(p),2*D)]
    poly += [field(p[i + D : i+2*D]) for i in range(0,len(p),2*D)]
    return poly
  
  def compare_hash(self, hash_p0, hash_p1):
    hash_0 = ctypes.cast(hash_p0, POINTER(c_uint64*4))
    hash_1 = ctypes.cast(hash_p1, POINTER(c_uint64*4))
    return list(hash_0.contents) == list(hash_1.contents)

  def sample_field_array(self, hash, size, p):
    # hash_l = (c_uint64*(4))()
    # hash_l[0] = int.from_bytes(hash[:8], "little")
    # hash_l[1] = int.from_bytes(hash[8:16], "little")
    # hash_l[2] = int.from_bytes(hash[16:24], "little")
    # hash_l[3] = int.from_bytes(hash[24:32], "little")
    return self.lib.sample_field_array(hash, size, p)
  
  # def load_proof_stream(self, proof_stream):
  #   size = len(proof_stream.objects)
  #   items = (c_uint64*(size))()
  #   hashes = (c_uint64*(size*4))()
  #   for i in range(size):
  #     items[i] = proof_stream.objects[i]
  #     hashes[4*i] = proof_stream.objects[i]
  #     hashes[4*i + 1] = proof_stream.objects[i]
  #     hashes[4*i + 2] = proof_stream.objects[i]
  #     hashes[4*i + 3] = proof_stream.objects[i]
  #   return self.lib.load_proof_stream(items, hashes, size)
  
  def decrypt_packed_points(self, points, size, field, sk):
    D = field.D
    p = (c_uint64*(size*D*2))()
    self.lib.decrypt_packed_points(p, points, size, sk)
    p = list(p)
    poly = [[field(p[i:i+D]), field(p[i+D: i+2*D])] for i in range(0,len(p),2*D)]
    return poly
  
  def decrypt_and_batch_points(self, points, size, beta_array, field, sk, num_poly):
    D = field.D
    p = (c_uint64*(size*D))()
    self.lib.decrypt_and_batch_points(p, points, size, beta_array, sk)
    p = list(p)
    n = len(p)//2
    poly = [[field(p[i:i+D]), field(p[n + i: n + i + D])] for i in range(0,n,D)]
    return poly
  
  def folding_level(self, codeword, const_array, N, r, pk):
    if(self.multithreaded and N/4 >= 4):
      num_threads = min(self.num_threads, (N//4)//2)
      return self.lib.folding_level_LWE_mt(codeword, const_array, N, r, pk, num_threads)
    else:
      return self.lib.folding_level_LWE(codeword, const_array, N, r, pk)

  def batch_codewords(self, codewords, size, beta_array, D, num_poly):
    if(self.multithreaded and size >= 4):
      num_threads = min(self.num_threads, size//2)
      return self.lib.batch_codewords_mt(codewords, size, beta_array, D, num_poly, num_threads)
    else:
      return self.lib.batch_codewords(codewords, size, beta_array, D, num_poly)

  def set_multi_threaded(self, num_thread):
    self.num_threads = num_thread
    if(num_thread > 1):
      self.multithreaded = True
    else:
      self.multithreaded = False



friC = FriC()
        