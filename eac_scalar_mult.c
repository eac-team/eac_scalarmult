/*
  eac_scalar_mult.c
  -----------------

  scalar multiplication using Euclidean addition chains of length 256
  on a curve with an efficient endomorphism

  E(F_p): y^2 = y^3 + 3
  p = 2^331 - 36301
  phi(x,y) = (beta*x,y)
  beta =  2^( (p-1)/3 ) mod p

  Basic GF(p) arithmetic using de Gnu MP library #include*/

#include <gmp.h>
#define P_331 "4374501449566023848745004454235242730706338861786424872851541212819905998398751846447026354046071347"
#define BETA "1547958758032931671736470639872108762198831514795455684792747161501219944666816365729042179453125869"
/* 
   Computes (p1, p1+p2) in co-Z arithmetic
   from p1=(X1,Y1,Z)  and p2=(X2,Y2,Z)
   
 */

void zaddu(mpz_t x1, mpz_t y1, mpz_t x2, mpz_t y2,  mpz_t z, mpz_t A,  mpz_t p){
  //~ A <- X2-X1 
  mpz_sub (A, x2, x1);
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
  //~ Y2 <- (Y2-Y1)^2 
  mpz_sub (y2, y2, y1);  	// Y2 <- (Y2-Y1)
  //~ X2 <- (Y2-Y1)^2  - X1 - A 
  mpz_mul (x2, y2, y2);		// X2 <- (Y2-Y1)^2
  mpz_mod (x2, x2, p);
  mpz_sub (x2, x2, x1);
  mpz_sub (x2, x2, A);
  mpz_mod (x2, x2, p); 
  //~ Y1 <- Y1*(A - X1) = Y1*(X2-X1)^3
  mpz_sub (A, A, x1);  		// A <- (A-X1)
  mpz_mul (y1, y1, A);
  mpz_mod (y1, y1, p);
  //~ Y2 <- Y2*(X1-X2)-Y1
  mpz_sub (A, x1, x2);  	// A <- (X1-X2)
  mpz_mul (y2, y2, A);
  mpz_sub (y2, y2, y1);
  mpz_mod (y2, y2, p);
}

void eac_scalarmult(mpz_t x1, mpz_t y1, const unsigned char select, const unsigned char *eac)
{
  int i,j,bit;
  mpz_t x2, y2, A, z, beta, p;

  mpz_init(A);
  mpz_init_set_ui(z,1);
  mpz_init_set_str(beta, BETA, 10);
  mpz_init_set_str(p, P_331, 10);

  /** START SCALAR MULT **/
  
  /*apply endomorphism*/
  mpz_init_set(x2,x1);
  mpz_mul(x2,x2,beta);
  mpz_mod(x2,x2,p);
  mpz_init_set(y2,y1);

  /*use EAC*/

  for (i=0; i<1; i++){
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
      gmp_printf("P=[%Zd, %Zd, %Zd]\n\n", x2,y2,z);
      gmp_printf("P=[%Zd, %Zd, %Zd]\n\n", x1,y1,z);
    }
  }
  zaddu(x2,y2,x1,y1,z,A,p);

  /* convert to affine*/
  gmp_printf("[%Zd, %Zd, %Zd]\n\n", x2,y2,z);
  mpz_invert(z,z,p);
  mpz_mul(A,z,z);
  mpz_mod(A,A,p);
  mpz_mul(x1,x1,A);
  mpz_mod(x1,x1,p);

  mpz_mul(A,A,z);
  mpz_mod(A,A,p);
  mpz_mul(y1,y1,A);
  mpz_mod(y1,y1,p);
  mpz_clears(x2, y2, A, z, beta, p, NULL);
}

int main()
{
  mpz_t x1, y1;
  unsigned char select = 1, eac[32] = {0,0,0,0,0,0,0,0,
				       0,0,0,0,0,0,0,0,
				       0,0,0,0,0,0,0,0,
				       0,0,0,0,0,0,0,0};
  
  mpz_init_set_str(x1,"2374220700659496913734031351069009197194002133691926985677385343504881882775867573371443296038832054",10);
  mpz_init_set_str(y1, "1160290385298866558973601951739419258017198715789204242005041424378088369659007871723823978946235971", 10);
  eac_scalarmult(x1,y1,select,eac);

  gmp_printf("[%Zd, %Zd]\n", x1,y1);
  mpz_clears(x1, y1, NULL);
  return 0;
}
