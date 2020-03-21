#ifndef zaddu_h
#define zaddu_h

#include <gmp.h>

extern 
void zaddu(mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2, mpz_t z, mpz_t A,  const mpz_t p);
void eac_end_256_smult(mpz_t x1, mpz_t y1, const unsigned char select, const unsigned char *eac, const mpz_t p, const mpz_t beta);
#endif
