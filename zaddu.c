#include "zaddu.h"
#include <stdio.h>

/*
  Computes (p1, p1+p2) in co-Z arithmetic
  from p1=(X1,Y1,Z)  and p2=(X2,Y2,Z)
*/

void zaddu(mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2,  mpz_t z, mpz_t A,  const mpz_t p){
  //~ A <- X2-X1 
  mpz_sub (A, x2, x1);
  mpz_mod (A, A, p); 
  //~ z <- z(X2-X1) 
  mpz_mul (z, z, A);
  mpz_mod (z, z, p); 
  //~ A <- (X2-X1)^2 
  mpz_mul (A, A, A);
  mpz_mod (A, A, p); 
  //~ X1 <- X1*(X2-X1)^2 
  mpz_mul (x1, x1, A);
  mpz_mod (x1, x1, p); 
  //~ A <- X2*(X2-X1)^2 
  mpz_mul (A, x2, A);
  mpz_mod (A, A, p);
  //~ Y2 <- (Y2-Y1)
  mpz_sub (y2, y2, y1);
  mpz_mod (y2, y2, p);	
  //~ X2 <- (Y2-Y1)^2  - X1 - A 
  mpz_mul (x2, y2, y2);		
  mpz_mod (x2, x2, p);
  mpz_sub (x2, x2, x1);
  mpz_mod (x2, x2, p);
  mpz_sub (x2, x2, A);
  mpz_mod (x2, x2, p); 
  //~ Y1 <- Y1*(A - X1) = Y1*(X2-X1)^3
  mpz_sub (A, A, x1);
  mpz_mod (A, A, p);  		
  mpz_mul (y1, y1, A);
  mpz_mod (y1, y1, p);
  //~ Y2 <- Y2*(X1-X2)-Y1
  mpz_sub (A, x1, x2);
  mpz_mod (A, A, p);  	 	
  mpz_mul (y2, y2, A);
  mpz_mod (y2, y2, p);
  mpz_sub (y2, y2, y1);
  mpz_mod (y2, y2, p);
}

/*
  Scalar multiplication using Euclidean Addition Chain <eac> of length
  256, starting from points (P, phi(P) ) where

  phi: (x,y) --> (beta*x,y)

  <select> is 0 or 1 and represent a "big step"

  B_EAC and P_EAC must be define in main program 
*/
void eac_end_256_smult(mpz_t x1, mpz_t y1, mpz_t z,
		       const unsigned char select, const unsigned char *eac,
		       const mpz_t p, const mpz_t beta)
{
  int i,j,bit;
  mpz_t x2, y2, A;

  mpz_init(A);

  /** START SCALAR MULT **/
  
  /*apply endomorphism*/
  mpz_init_set(x2,x1);
  mpz_mul(x2,x2,beta);
  mpz_mod(x2,x2,p);
  mpz_init_set(y2,y1);

  /*use EAC*/

  for (i=0; i<32; i++){
    for (j=0; j<8; j++){

      bit = (eac[i] & (1<<j))>>j;

      mpz_xor(A,x1,x2);
      mpz_mul_ui(A, A, bit^select);
      mpz_xor(x1,x1,A);
      mpz_xor(x2,x2,A);
      mpz_xor(A,y1,y2);
      mpz_mul_ui(A, A, bit^select);
      mpz_xor(y1,y1,A);
      mpz_xor(y2,y2,A);
       
      zaddu(x1,y1,x2,y2,z,A,p);
    }
  }
  zaddu(x2,y2,x1,y1,z,A,p);

  mpz_clears(x2, y2, A, NULL);
}

void eac_to_affine(mpz_t x1, mpz_t y1, mpz_t z, const mpz_t p){
 
  mpz_invert(z,z,p);
  
  mpz_mul(x1,x1,z);
  mpz_mod(x1,x1,p);
  mpz_mul(x1,x1,z);
  mpz_mod(x1,x1,p);
  
  mpz_mul(y1,y1,z);
  mpz_mod(y1,y1,p);
  mpz_mul(y1,y1,z);
  mpz_mod(y1,y1,p);
  mpz_mul(y1,y1,z);
  mpz_mod(y1,y1,p);
 
}

void zaddu_x_only(mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2, mpz_t A,  const mpz_t p){
  //~ A <- X2-X1 
  mpz_sub (A, x2, x1);
  mpz_mod (A, A, p); 
  //~ z <- z(X2-X1) 
  //mpz_mul (z, z, A);
  //mpz_mod (z, z, p); 
  //~ A <- (X2-X1)^2 
  mpz_mul (A, A, A);
  mpz_mod (A, A, p); 
  //~ X1 <- X1*(X2-X1)^2 
  mpz_mul (x1, x1, A);
  mpz_mod (x1, x1, p); 
  //~ A <- X2*(X2-X1)^2 
  mpz_mul (A, x2, A);
  mpz_mod (A, A, p);
  //~ Y2 <- (Y2-Y1)
  mpz_sub (y2, y2, y1);
  mpz_mod (y2, y2, p);	
  //~ X2 <- (Y2-Y1)^2  - X1 - A 
  mpz_mul (x2, y2, y2);		
  mpz_mod (x2, x2, p);
  mpz_sub (x2, x2, x1);
  mpz_mod (x2, x2, p);
  mpz_sub (x2, x2, A);
  mpz_mod (x2, x2, p); 
  //~ Y1 <- Y1*(A - X1) = Y1*(X2-X1)^3
  mpz_sub (A, A, x1);
  mpz_mod (A, A, p);  		
  mpz_mul (y1, y1, A);
  mpz_mod (y1, y1, p);
  //~ Y2 <- Y2*(X1-X2)-Y1
  mpz_sub (A, x1, x2);
  mpz_mod (A, A, p);  	 	
  mpz_mul (y2, y2, A);
  mpz_mod (y2, y2, p);
  mpz_sub (y2, y2, y1);
  mpz_mod (y2, y2, p);
}

void eac_end_256_smult_x_only(mpz_t x1, mpz_t y1, 
		       const unsigned char select, const unsigned char *eac,
		       const mpz_t p, const mpz_t beta)
{
  int i,j,bit;
  mpz_t x2, y2, A;

  mpz_init(A);

  /** START SCALAR MULT **/
  
  /*apply endomorphism*/
  mpz_init_set(x2,x1);
  mpz_mul(x2,x2,beta);
  mpz_mod(x2,x2,p);
  mpz_init_set(y2,y1);

  /*use EAC*/

  for (i=0; i<32; i++){
    for (j=0; j<8; j++){

      bit = (eac[i] & (1<<j))>>j;

      mpz_xor(A,x1,x2);
      mpz_mul_ui(A, A, bit^select);
      mpz_xor(x1,x1,A);
      mpz_xor(x2,x2,A);
      mpz_xor(A,y1,y2);
      mpz_mul_ui(A, A, bit^select);
      mpz_xor(y1,y1,A);
      mpz_xor(y2,y2,A);
       
      zaddu_x_only(x1,y1,x2,y2,A,p);
    }
  }
  zaddu_x_only(x2,y2,x1,y1,A,p);

  mpz_clears(x2, y2, A, NULL);
}
