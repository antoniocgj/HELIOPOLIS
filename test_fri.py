from algebra import *
from fri import *
import time, random, argparse, sys
BENCH = False

if(BENCH):
    import yep

DEBUG = False

def test_fri(input_size, expansion_factor, num_colinearity_tests, num_threads):
    field = Field.optFP16()
    degree = 2**input_size - 1
    N = 4096
    # N = 8192
    batching_size = N
    print("Testing FRI with %d polynomials of degree %d with RS code rate 1/%d and %d colinearity tests per round." % (batching_size, degree, expansion_factor, num_colinearity_tests))

    friC.set_multi_threaded(num_threads)
    

    initial_codeword_length = (degree + 1) * expansion_factor
    log_codeword_length = 0
    codeword_length = initial_codeword_length
    while codeword_length > 1:
        codeword_length //= 2
        log_codeword_length += 1

    omega = FieldElement(field.primitive_nth_root(initial_codeword_length))

    print("\nGenerating FHE keys")
    sk = friC.gen_private_key(N, initial_codeword_length)
    pk = friC.gen_public_key(field, sk)

    fri = Fri(field, omega, initial_codeword_length, expansion_factor, num_colinearity_tests, public_key=pk)


    print("\nGenerating %d random polynomials and encrypting" % batching_size)
    enc_polynomials = [random.randint(0, 2**16+1) for i in range(N*(degree + 1))] 
    # for i in range((degree + 1)):
    #     # enc_polynomials += [random.randint(0, 2**16+1) for i in range(batching_size)] + ([0]*(N-batching_size))
    #     enc_polynomials += [i] + ([0]*(N-1))
    enc_input = friC.encrypt_input(enc_polynomials, expansion_factor, sk)

    print("\nDepth-2 homomorphic NTT of size %d for %d encrypted polynomials." % (initial_codeword_length, batching_size), end=" ")
    startNTT = time.time()
    friC.shallow_ntt(enc_input, field, initial_codeword_length)
    endNTT = time.time()
    print("Time: %0.3fs" % (endNTT - startNTT))
    if(False):
        polynomial = Polynomial([field(i) for i in range(degree+1)])
        codeword_dec = friC.decrypt_poly(enc_input, initial_codeword_length, 0, sk, field)
        codeword = shallow_expand( polynomial, omega, initial_codeword_length )
        print("NTT test: ",  codeword_dec == codeword)
        # print("input:", codeword_dec)
        if(codeword_dec != codeword):
            print("check:", codeword[:20])
            exit(0)

    # test valid codeword
    print("\nStarting FRI proof using (up to) %d threads" % num_threads)
    proof_stream = ProofStream(4*num_colinearity_tests*log2(initial_codeword_length))
    startProver = time.time()
    fri.prove(enc_input, initial_codeword_length, batching_size, proof_stream)
    endProver = time.time()
    print("Finish FRI Proof. Total time %0.3fs\n" % (endProver - startProver))

    print("\nStarting FRI verification")
    startVerify = time.time()
    if(BENCH): yep.start('verifier.prof')
    verdict = fri.verify_optimized(initial_codeword_length, batching_size, proof_stream, sk)
    if(BENCH): yep.stop()
    endVerify = time.time()
    print("Finish FRI verification. Total time: %0.3fs\n" % (endVerify - startVerify))


    if verdict == False:
        print("rejecting proof, but proof should be valid!")
        return
    
    print("success! \\o/")
    print("Proof time: %0.3fs" % (endNTT - startNTT + endProver - startProver))
    print("Verify time: %0.3fs" % (endVerify - startVerify))

    import resource
    print("Maximum memory usage (ru_maxrss): %d KB " % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluates FRI over encrypted polynomials")
    parser.add_argument('-i', '--input_size', dest='input_size', type=int, default=7,
                        help='Log_2 of the degree bound of input polynomials (default=7)')
    parser.add_argument('-e', '--expansion', dest='expansion_factor', type=int, action='store', default=2, help='Expansion factor (1 / \\rho) (default=2)')
    parser.add_argument('-m', '--colin_tests', dest='num_colinearity_tests', type=int, action='store', default=102, help='Number of colinearity tests per round of FRI (default=102)')
    parser.add_argument('-t', '--threads', dest='threads', type=int, action='store', default=2, help='Number of threads for the prover (default=2)')
    args = parser.parse_args()
    
    test_fri(args.input_size, args.expansion_factor, args.num_colinearity_tests, args.threads)