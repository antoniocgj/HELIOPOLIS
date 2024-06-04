from hashlib import blake2b
import pickle as pickle # serialization
from friC import friC

class Merkle:
    def commit(hash_array):
        return friC.merkle_commit(hash_array)

    def open(index, hash_array):
        return friC.merkle_open(hash_array, index)

    def verify(root, index, path, element_hash):
        return friC.merkle_verify(root, index, path, element_hash)


# class Merkle:
#     H = blake2b

#     def commit_( leafs ):
#         assert(len(leafs) & (len(leafs)-1) == 0), "length must be power of two"
#         if len(leafs) == 1:
#             return leafs[0]
#         else:
#             return Merkle.H(Merkle.commit_(leafs[:len(leafs)//2]) + Merkle.commit_(leafs[len(leafs)//2:]), digest_size=32).digest()

#     def commit( data_array ):
#         return Merkle.commit_([Merkle.H(bytes(da), digest_size=32).digest() for da in data_array])
    
#     def commit_hasharray(hash_array):
#         return Merkle.commit_([pickle.dumps(ha) for ha in hash_array])
    
#     def open_( index, leafs ):
#         assert(len(leafs) & (len(leafs)-1) == 0), "length must be power of two"
#         assert(0 <= index and index < len(leafs)), "cannot open invalid index"
#         if len(leafs) == 2:
#             return [leafs[1 - index]]
#         elif index < (len(leafs)/2):
#             return Merkle.open_(index, leafs[:len(leafs)//2]) + [Merkle.commit_(leafs[len(leafs)//2:])]
#         else:
#             return Merkle.open_(index - len(leafs)//2, leafs[len(leafs)//2:]) + [Merkle.commit_(leafs[:len(leafs)//2])]

#     def open( index, data_array ):
#         return Merkle.open_(index, [Merkle.H(bytes(da), digest_size=32).digest() for da in data_array])
    
#     def open_hasharray( index, hash_array ):
#         return Merkle.open_(index, [pickle.dumps(ha) for ha in hash_array])
    
#     def verify_( root, index, path, leaf ):
#         assert(0 <= index and index < (1 << len(path))), "cannot verify invalid index"
#         if len(path) == 1:
#             if index == 0:
#                 return root == Merkle.H(leaf + path[0], digest_size=32).digest()
#             else:
#                 return root == Merkle.H(path[0] + leaf, digest_size=32).digest()
#         else:
#             if index % 2 == 0:
#                 return Merkle.verify_(root, index >> 1, path[1:], Merkle.H(leaf + path[0], digest_size=32).digest())
#             else:
#                 return Merkle.verify_(root, index >> 1, path[1:], Merkle.H(path[0] + leaf, digest_size=32).digest())

#     def verify( root, index, path, data_element ):
#         return Merkle.verify_(root, index, path, Merkle.H(bytes(data_element), digest_size=32).digest())
    
#     def verify_from_hash( root, index, path, hash ):
#         return Merkle.verify_(root, index, path, pickle.dumps(hash))

