#include <mosfhet.h>

// extern TRLWE_Key _glb_debug_trlwe_key;
// extern TRLWE_Key _glb_debug_trlwe_key2;
// extern TLWE_Key _glb_debug_tlwe_key;

// static void _debug_print128(char * msg, __uint128_t x){
//   uint64_t * x64 = (uint64_t *) &x;
//   printf("%s: 0x%016lx%016lx\n", msg, x64[1], x64[0]);
// }

// static void _debug_print_tlwe(TLWE c, uint64_t scale){
//   Torus v = tlwe_phase(c, _glb_debug_tlwe_key);
//   printf("%lu\n", torus2int(v, scale));
// }

static void _debug_print_trlwe(TRLWE c, uint64_t scale, TRLWE_Key key){
  TorusPolynomial p = polynomial_new_torus_polynomial(c->b->N);
  trlwe_phase(p, c, key);
  for (size_t i = 0; i < c->b->N; i++){
    printf("%lu, ", torus2int(p->coeffs[i], scale));
  }
  printf("\n");
}

static void _debug_print_trlwe2(TRLWE c, uint64_t mod, TRLWE_Key key){
  TorusPolynomial p = polynomial_new_torus_polynomial(c->b->N);
  trlwe_phase(p, c, key);
  for (size_t i = 0; i < c->b->N; i++){
    const uint64_t val = mod_switch(torus2int(p->coeffs[i], 52), 1ULL<<52, mod)%mod;
    printf("%lu, ", val);
  }
  printf("\n");
}