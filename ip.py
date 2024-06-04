from hashlib import blake2b
import pickle as pickle # serialization
from friC import friC
import ctypes

class ProofStream:
    def __init__( self, max_size ):
        self.max_size = int(max_size)
        self.ps_obj = friC.lib.proof_stream_setup(self.max_size)

    def push(self, obj, hash=None):
        hash_p = ctypes.cast(hash, ctypes.POINTER(ctypes.c_uint64))
        friC.lib.proof_stream_push(self.ps_obj, obj, hash_p)

    def pull(self):
        return friC.lib.proof_stream_pull(self.ps_obj)
    
    def prover_fiat_shamir(self):
        return friC.lib.proof_stream_fs(self.ps_obj)

    def verifier_fiat_shamir(self):
        return friC.lib.proof_stream_read_fs(self.ps_obj)

# class ProofStream:
#     def __init__( self ):
#         self.objects = []
#         self.hashs = []
#         self.read_index = 0
#         self.hash_wsize = 4
#         self.hash_bytes = b''

#     def push( self, obj , hash = b'0'*32):
#         self.objects += [obj]
#         self.hashs += [hash]
#         if(type(hash) == int):
#             self.hash_bytes += hash.to_bytes(32, "little")
#         else:
#             assert(len(hash) == 32)
#             self.hash_bytes += hash


#     def pull( self ):
#         assert(self.read_index < len(self.objects)), "ProofStream: cannot pull object; queue empty."
#         obj = self.objects[self.read_index]
#         self.read_index += 1
#         return obj

#     def serialize( self ):
#         assert(False)
#         return pickle.dumps(self.objects)

#     def prover_fiat_shamir( self, num_bytes=32 ):
#         return blake2b(self.hash_bytes, digest_size=num_bytes).digest()

#     def verifier_fiat_shamir( self, num_bytes=32 ):
#         return blake2b(self.hash_bytes[:self.read_index*self.hash_wsize*8], digest_size=num_bytes).digest()

#     def deserialize( self, bb ):
#         assert(False)
#         ps = ProofStream()
#         ps.objects = pickle.loads(bb)
#         return ps

