# HELIOPOLIS

HELIOPOLIS: Verifiable Computation over Homomorphically Encrypted Data from Interactive Oracle Proofs is Practical

### Paper

[eprint 2023/1949](https://eprint.iacr.org/2023/1949)

```
@misc{cryptoeprint:2023/1949,
      author = {Diego F. Aranha and Anamaria Costache and Antonio Guimarães and Eduardo Soria-Vazquez},
      title = {HELIOPOLIS: Verifiable Computation over Homomorphically Encrypted Data from Interactive Oracle Proofs is Practical},
      howpublished = {Cryptology ePrint Archive, Paper 2023/1949},
      year = {2023},
      note = {\url{https://eprint.iacr.org/2023/1949}},
      url = {https://eprint.iacr.org/2023/1949}
}
```



## Requirements
- Basic development packages (On Debian-based systems: `sudo apt install make gcc cmake`)
- Python 3 (Tested on Python 3.10.12)
- GCC (Tested on GCC 11.4.0)
- Strongly recommended: A processor with [AVX-512 IFMA](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#AVX-IFMA) (Tested on Intel i7-1165G7 and Intel Xeon 8375C). We have added a compilation option that does not require it, but please be aware that it is mostly unoptimized and may be more unstable. 

## Running

1) Clone this repository

```
git clone https://github.com/antoniocgj/HELIOPOLIS.git
```

2) Run `python3 test_fri.py -h` to see the following help text:

```
usage: test_fri.py [-h] [-i INPUT_SIZE] [-e EXPANSION_FACTOR] [-m NUM_COLINEARITY_TESTS] [-t THREADS]

Evaluates FRI over encrypted polynomials

options:
  -h, --help            show this help message and exit
  -i INPUT_SIZE, --input_size INPUT_SIZE
                        Log_2 of the degree bound of input polynomials (default=7)
  -e EXPANSION_FACTOR, --expansion EXPANSION_FACTOR
                        Expansion factor (1 / \rho) (default=2)
  -m NUM_COLINEARITY_TESTS, --colin_tests NUM_COLINEARITY_TESTS
                        Number of colinearity tests per round of FRI (default=102)
  -t THREADS, --threads THREADS
                        Number of threads for the prover (default=2)
```

The Python code will compile the C library during the first execution. 

To recompile the code (in case you want to run it in another machine, for example), delete the build:
```
rm -rf c_lib/lib/libfric.so
rm -rf c_lib/src/third-party/hexl/build/
```

### Execution examples:

1. Parameter set FRI<sub>0</sub> with 4 threads for input of size 2<sup>11</sup>:

```python3 test_fri.py -i 11 -e 2 -m 102 -t 4```



2. Parameter set FRI<sub>3</sub> with 4 threads for input of size 2<sup>7</sup>:
   
```python3 test_fri.py -i 7 -e 16 -m 26 -t 4```


> [!WARNING]
> The optimized version of this implementation requires [AVX-512 IFMA](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#AVX-IFMA). It should be significantly slower and possibly unstable without it (but it should still work).

## Code structure
- `test_fry.py`: main script, it instantiates and runs FRI.
- `friC.py`: ctypes wrapper to call optimized (C/CPP-implemented) functions
- other Python scripts: high-level implementations of HE-FRI and all subprocedures required to implement it (Merkle commitments, NTTs, Fiat-Shamir...).
- `c_lib/src`: optimized C/CPP subroutines used by the Python scripts. It also contains a fully optimized version of HE-FRI's verification (`c_lib/src/verifier.cpp`).

## Implementation changelog

Compared to the first version of the paper:

- Removed dependency on HELib (Replaced with HEXL)
- Removed dependency on Sagemath
- Readjusted HE parameters to make the prover faster at a small cost for the verifier.
- Added a multi-threaded version for the folding and batching operations. 

## License

[Apache License Version 2.0](./LICENSE)

Adapted from the FRI implementation of [Stark Anatomy](https://github.com/aszepieniec/stark-anatomy) - [Apache License Version 2.0](https://github.com/aszepieniec/stark-anatomy?tab=Apache-2.0-1-ov-file#readme) -  Copyright Alan Szepieniec

Additionally, this repository contains code from:

- [MOSFHET](https://github.com/antoniocgj/MOSFHET): [Apache License Version 2.0](https://github.com/antoniocgj/MOSFHET/blob/main/LICENSE) - Copyright Antonio Guimarães et al. - See their [detailed copyright information](https://github.com/antoniocgj/MOSFHET/tree/main?tab=readme-ov-file#license).
- [Intel HEXL](https://github.com/intel/hexl): [Apache License 2.0](https://github.com/intel/hexl/blob/development/LICENSE) - Copyright 2020 Intel Corporation
- [BLAKE3](https://github.com/BLAKE3-team/BLAKE3): [Apache License 2.0](https://github.com/BLAKE3-team/BLAKE3/blob/master/LICENSE_A2) - Copyright 2019 Jack O'Connor and Samuel Neves
