import code
from algebra import *
from merkle import *
from ip import *
from ntt import *
from binascii import hexlify, unhexlify
import math
from hashlib import blake2b
from univariate import *
from friC import *
import time, sys

class Fri:
    def __init__( self, field, omega, initial_domain_length, expansion_factor, num_colinearity_tests, public_key=None, secret_key=None, early_stop=-1):
        self.field = field
        self.early_stop = early_stop
        self.offset = field(1)
        self.omega = omega
        self.inv_omega = omega.inverse()
        self.omega_powers = self.get_powers(omega, initial_domain_length)
        self.inv_omega_powers = self.get_powers(self.inv_omega, initial_domain_length)
        self.domain_length = initial_domain_length
        self.expansion_factor = expansion_factor
        self.num_colinearity_tests = num_colinearity_tests
        
        self.pk = public_key
        self.sk = secret_key
        print("Setting up FRI with %d rounds" % self.num_rounds())
        assert(self.num_rounds() >= 1), "cannot do FRI with less than one round"

    def num_rounds( self ):
        codeword_length = self.domain_length
        num_rounds = 0
        while codeword_length > self.expansion_factor and 2*self.num_colinearity_tests < codeword_length:
            codeword_length /= 2
            num_rounds += 1
            if(self.early_stop != -1 and codeword_length < self.early_stop):
                break
        return num_rounds + 1

    def sample_index( byte_array, size ):
        acc = 0
        for b in byte_array:
            acc = (acc << 8) ^ int(b)
        return acc % size
    
    def sample_array(self, byte_array, size):
        byte_len = self.field.byte_len
        return [self.field.sample(byte_array[byte_len*i:byte_len*(i+1)]) for i in range(size)]

    def sample_indices( self, seed, size, reduced_size, number ):
        assert(number <= reduced_size), f"cannot sample more indices than available in last codeword; requested: {number}, available: {reduced_size}"
        assert(number <= 2*reduced_size), "not enough entropy in indices wrt last codeword"

        indices = []
        reduced_indices = []
        counter = 0
        while len(indices) < number:
            index = Fri.sample_index(blake2b(seed + bytes(counter)).digest(), size)
            reduced_index = index % reduced_size
            counter += 1
            if reduced_index not in reduced_indices:
                indices += [index]
                reduced_indices += [reduced_index]

        return indices

    def eval_domain( self ):
        return [self.offset * (self.omega_powers[i]) for i in range(self.domain_length)]
    
    def prepare_const_array(self, const_array):
        d = self.field.D
        point_type = c_uint64 * d
        res_type = ctypes.POINTER(c_uint64) * len(const_array)
        res = res_type()
        for i in range(len(const_array)):
            point = self.field.to_list(const_array[i])
            point_pointer = ctypes.cast(point_type(*point), ctypes.POINTER(ctypes.c_ulong)) 
            res[i] = point_pointer
        return res
    
    def get_powers(self, x, size):
        return self.field.get_powers(x, size)
        res = [0]*size
        res[0] = self.field(1)
        for i in range(1, size):
            res[i] = res[i - 1]*x
        return res

    def commit( self, codewords, size, num_poly, proof_stream):
        one = self.field.one()
        two = self.field(2)
        inv_two = two**-1
        omega_index = 1
        D = self.field.D

        consts_array = [one]*size

        # repack and hash first codeword
        start = time.time()
        codeword = friC.lib.repack_init_codeword(codewords, size)
        enc_cw_hashes = friC.hash_codeword(codeword, size, return_pointer=True)
        end = time.time()
        print("Commit to the input. Time: %0.3fs" % (end - start))

        root = Merkle.commit(enc_cw_hashes)
        proof_stream.push(root, root)

        # get batching challenge 
        beta_list = friC.sample_field_array(proof_stream.prover_fiat_shamir(), num_poly*self.field.D, self.field.p)
        
        
        # batching
        print("Batching polynomials... ", end="")
        sys.stdout.flush()
        start = time.time()
        encrypted_first_codeword = friC.batch_codewords(codewords, size, beta_list, D, num_poly)
        end = time.time()
        print("Time: %0.3fs" % (end - start))

        N = size

        codewords = []
        hash_list = []
        
        # for each round
        for r in range(self.num_rounds()):
            # prepare next round, but only if necessary
            if r == self.num_rounds() - 1:
                break

            # get challenge
            alpha = self.field.sample(proof_stream.prover_fiat_shamir())

            # collect codeword and hashs
            codewords += [codeword]
            hash_list += [enc_cw_hashes]

            # split and (shallow) fold
            # constant processing
            for i in range(N//2):
                inv_omega_power = self.inv_omega_powers[omega_index*i]*alpha
                for j in range(2**(r + 1)):
                    if(j&1): consts_array[i + j*N//2] *= inv_two*(one - inv_omega_power)
                    else: consts_array[i + j*N//2] *= inv_two*(one + inv_omega_power)

            consts_array_pointer = self.prepare_const_array(consts_array)
            print("Folding round %d." % r, end=" ")
            sys.stdout.flush()
            start = time.time()
            codeword = friC.folding_level(encrypted_first_codeword, consts_array_pointer, N, r, self.pk)
            end = time.time()
            print("Time: %0.3fs" % (end - start))
            if(r == self.num_rounds() - 2):
                enc_cw_hashes = friC.hash_full_codeword(codeword)
                proof_stream.push(enc_cw_hashes, enc_cw_hashes)
            else:
                enc_cw_hashes = friC.hash_codeword(codeword, N//4, return_pointer=True)
                root = Merkle.commit(enc_cw_hashes)
                proof_stream.push(root, root)

            omega_index *= 2
            N = N//2

        # send last codeword
        proof_stream.push(codeword, enc_cw_hashes)

        # collect last codeword too
        codewords = codewords + [codeword]
        hash_list = hash_list + [enc_cw_hashes]

        return codewords, hash_list

    def query( self, current_codeword, current_hash_list, selected_indices, proof_stream):
        vector_type = c_uint64 * len(selected_indices)
        points = friC.lib.codeword_get_packed_points(current_codeword, vector_type(*selected_indices), len(selected_indices))
        
        proof_stream.push(points)

        # reveal authentication paths
        for s in range(len(selected_indices)):
            mk_open =  Merkle.open(selected_indices[s], current_hash_list)
            proof_stream.push(mk_open)
        

    def prove( self, codeword, size, num_poly, proof_stream):
        assert(self.domain_length == size), "initial codeword length does not match length of initial codeword"

        # commit phase
        codewords, hash_list = self.commit(codeword, size, num_poly, proof_stream)

        # get indices
        top_level_indices = friC.sample_indices(proof_stream.prover_fiat_shamir(), size, size>>(self.num_rounds() - 1), self.num_colinearity_tests)
        indices = [index for index in top_level_indices]

        self.selected_points = []

        # query phase
        indices = [index % (size//2) for index in indices] # fold
        selected_indices = indices + [idx + (size//2) for idx in indices]
        self.query(codewords[0], hash_list[0], selected_indices, proof_stream)

        for i in range(1, len(codewords)-1):
            indices = [index % (size>>(i+1)) for index in indices] # fold
            self.query(codewords[i], hash_list[i], indices, proof_stream)

        return top_level_indices
    
    def verify_optimized(self, size, num_poly, proof_stream, sk):
        omega_int = int(self.field.to_list(self.omega)[0])
        return friC.verify_native(proof_stream.ps_obj, size, self.num_colinearity_tests, self.num_rounds(), omega_int, self.expansion_factor, sk)

    def verify( self, size, num_poly, proof_stream, sk):
        he_operations_time = 0
        omega_index = 1
        field = self.field
        num_rounds = self.num_rounds()
        inv_two = field(2)**-1

        # extract all roots, alphas, and beta
        roots = [proof_stream.pull()]
        alphas = [self.field.sample(proof_stream.verifier_fiat_shamir())]

        beta_list = friC.sample_field_array(proof_stream.verifier_fiat_shamir(), num_poly*self.field.D, self.field.p)
        for r in range(self.num_rounds() - 1):
            roots += [proof_stream.pull()]
            alphas += [self.field.sample(proof_stream.verifier_fiat_shamir())]
        
        # extract last codeword
        last_codeword_enc = proof_stream.pull()
        last_codeword_size = size >> (num_rounds - 1)

        start = time.time()
        last_codeword_hash = friC.hash_full_codeword(last_codeword_enc)
        end = time.time()
        he_operations_time += (end - start)

        # check if it matches the given root
        if not friC.compare_hash(roots[-1], last_codeword_hash):
            print("last codeword is not well formed")
            return False
        

        # Decrypt and check if it is low degree
        degree = (last_codeword_size // self.expansion_factor) - 1
        last_omega = self.omega_powers[2**(self.num_rounds()-1)]
        start = time.time()
        last_codeword = friC.decrypt_codeword(last_codeword_enc, last_codeword_size, self.field, sk)
        end = time.time()
        he_operations_time += (end - start)
        
        poly = Polynomial(intt(last_omega, last_codeword))
        if poly.degree() > degree:
            print("last codeword does not correspond to polynomial of low enough degree")
            print("observed degree:", poly.degree())
            print("but should be:", degree)
            print("Polynomial:", poly)
            return False

        # get indices
        top_level_indices = friC.sample_indices(proof_stream.verifier_fiat_shamir(), self.domain_length, self.domain_length >> (self.num_rounds()-1), self.num_colinearity_tests)

        # for every round, check consistency of subsequent layers
        for r in range(0, self.num_rounds()-1):

            # fold c indices
            c_indices = [index % (self.domain_length >> (r+1)) for index in top_level_indices]

            # infer a and b indices
            a_indices = [index for index in c_indices]
            b_indices = [index + (self.domain_length >> (r+1)) for index in a_indices]

            # get eval points
            points = proof_stream.pull()

            if r == 0: ## batching

                start = time.time()
                # decrypt and batch first round points
                dec_points = friC.decrypt_and_batch_points(points, self.num_colinearity_tests*2, beta_list, field, sk, num_poly)

                # get the hash to verify later
                points_hash = friC.hash_packed_points(points, self.num_colinearity_tests*2)
                end = time.time()
                he_operations_time += (end - start)

                prev_index = a_indices.copy()
                dec_points_prev = dec_points
            else:
                start = time.time()
                dec_points = friC.decrypt_packed_points(points, self.num_colinearity_tests, field, sk)
                points_hash = friC.hash_packed_points(points, self.num_colinearity_tests)
                end = time.time()
                he_operations_time += (end - start)

                
                for s in range(self.num_colinearity_tests):
                    ay = dec_points_prev[s][0]
                    by = dec_points_prev[s][1]
                    cay = dec_points[s][0]
                    cby = dec_points[s][1]

                    inv_omega_power = self.inv_omega_powers[omega_index*prev_index[s]]

                    prev_layer_value = inv_two*(ay + by) + alphas[r - 1]*inv_two*(ay - by)*inv_omega_power
                    if(prev_layer_value != cay and prev_layer_value != cby):
                        print("colinearity check failure")
                        print("ay: ", self.field.to_list(ay))
                        print("by: ", self.field.to_list(by))
                        print("inv omega: ", inv_omega_power, "inv_two:", inv_two, "alpha[r-1]", self.field.to_list(alphas[r - 1]))
                        print("s: ", s, "r:", r, "idx", prev_index[s], c_indices[s])
                        print("omega idx: ", omega_index)
                        print("prev_layer_value: ", self.field.to_list(prev_layer_value))
                        print("cay: ", self.field.to_list(cay))
                        print("cay: ", self.field.to_list(cay))
                        sys.stdout.flush()
                        return False
                    
                dec_points_prev = dec_points
                prev_index = a_indices.copy()

                # square omega and offset to prepare for next round
                omega_index *= 2
                # offset = offset**2


            # verify authentication paths
            for i in range(self.num_colinearity_tests):
                path = proof_stream.pull()
                if Merkle.verify(roots[r], a_indices[i], path, points_hash[i]) == False:
                    print("merkle authentication path verification fails for points_hash a")
                    return False
                
            if r == 0:
                for i in range(self.num_colinearity_tests):
                    path = proof_stream.pull()
                    if Merkle.verify(roots[r], b_indices[i], path, points_hash[self.num_colinearity_tests + i]) == False:
                        print("merkle authentication path verification fails for points_hash b")
                        return False
        
        # colin check with the last codeword 
        # (which was already decrypted and verified at the beginning)
        c_indices = [index % (last_codeword_size//2) for index in top_level_indices]
        for s in range(self.num_colinearity_tests):
            ay = dec_points_prev[s][0]
            by = dec_points_prev[s][1]

            inv_omega_power = self.inv_omega_powers[omega_index*prev_index[s]]

            prev_layer_value = inv_two*(ay + by) + alphas[r]*inv_two*(ay - by)*inv_omega_power
            
            if(last_codeword[c_indices[s]] != prev_layer_value and
               last_codeword[c_indices[s]+last_codeword_size//2] != prev_layer_value):
                print("last colinearity check failure")
                print("ay: ", self.field.to_list(ay))
                print("by: ", self.field.to_list(by))
                print("inv_omega_power: ", inv_omega_power, "inv_two:", inv_two, "alpha[r]", self.field.to_list(alphas[r]))
                print("idx", prev_index[s], c_indices[s])
                print("omega idx: ", omega_index)
                print("prev_layer_value: ", self.field.to_list(prev_layer_value))
                print("last_codeword a: ", last_codeword[c_indices[s]])
                print("last_codeword b: ", last_codeword[c_indices[s]+last_codeword_size//2])
                sys.stdout.flush()
                return False               

        # all checks passed
        print("Cummulative time of HE operations: %0.3fs" % he_operations_time)
        sys.stdout.flush()
        return True

