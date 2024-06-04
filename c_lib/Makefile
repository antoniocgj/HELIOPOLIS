MOSFHET_DIR = ./mosfhet/
FAST_DEBUG = 1
include $(MOSFHET_DIR)/Makefile.def
CC = g++
INCLUDE_DIRS += ./include/ ./src/third-party/hexl/hexl/include/
OPT_FLAGS += -Wno-sign-compare -Wno-write-strings 
DEBUG_FLAGS += -Wno-sign-compare -Wno-write-strings 
LIB_FLAGS += -Wno-sign-compare -Wno-write-strings 
LIB_DIRS = ./src/third-party/hexl/build/hexl/lib 
LD_LIBS = hexl
LIBS += $(addprefix -L, $(LIB_DIRS)) $(addprefix -l, $(LD_LIBS)) 

SRC=folding.cpp encrypt.cpp  homNTT.cpp field.cpp verifier.cpp merkle.cpp ip.cpp hash.cpp
BLAKE_SRC=blake3.c blake3_dispatch.c blake3_portable.c blake3_sse2_x86-64_unix.S blake3_sse41_x86-64_unix.S blake3_avx2_x86-64_unix.S blake3_avx512_x86-64_unix.S
RNS_SRC=polynomial.cpp rlwe.cpp misc.cpp
MP_RLWE_SRC=mpLWE.c

ALL_SRC = $(addprefix ./src/, $(SRC)) $(addprefix ./blake3/, $(BLAKE_SRC)) $(addprefix ./src/rns-rlwe/, $(RNS_SRC)) 

all: lib lib/libfric

mainLocalVer: $(SRC_MOSFHET) $(ALL_SRC) mainLocalVer.c
	$(CC) -g -o mainLocalVer $^ $(OPT_FLAGS) $(LIBS)

mainFolding: $(SRC_MOSFHET) $(ALL_SRC) mainFolding.c
	$(CC) -g -o mainFolding $^ $(OPT_FLAGS) $(LIBS)

testMP: $(SRC_MOSFHET) $(ALL_SRC) testMP.c
	$(CC) -g -o testMP $^ $(OPT_FLAGS) $(LIBS)

testIOP: $(SRC_MOSFHET) $(ALL_SRC) testIOP.cpp
	$(CC) -g -o testIOP $^ $(DEBUG_FLAGS) $(LIBS)

testRNS: $(SRC_MOSFHET) $(ALL_SRC) testRNS.cpp
	$(CC) -g -o testRNS $^ $(DEBUG_FLAGS) $(LIBS)

lib/libfric: $(SRC_MOSFHET) $(ALL_SRC) 
	$(CC) -g -fPIC -shared -o lib/libfric.so $^ $(LIB_FLAGS) $(LIBS)


hexl: hexl/build
	cmake --build ./src/third-party/hexl/build

hexl/build:
	cmake -S ./src/third-party/hexl/ -B ./src/third-party/hexl/build

ntt:
	make -f Makefile2

lib:
	mkdir -p lib

