#include <stdio.h>
#include "zaddu.h"
#include "eac_param_358.h"


int main() {

  FILE * fp = fopen(DATA, "r");
  mpz_t x,y,x1, y1, x2, y2, z, p, beta;
  char sx1[FIELD_BASE_10_SIZE+1], sy1[FIELD_BASE_10_SIZE+1],
    sx2[FIELD_BASE_10_SIZE+1], sy2[FIELD_BASE_10_SIZE+1],
    s3[SCALAR_BIT_SIZE+1]; 
  unsigned char eac[32];
  unsigned int select;
  int j,k;
  
  mpz_init(x);
  mpz_init(y);
  mpz_init(x1);
  mpz_init(y1);
  mpz_init(x2);
  mpz_init(y2);
  
  mpz_init_set_str(p,P_EAC,10);
  mpz_init_set_str(beta,B_EAC,10);

  for(int i=0;i<NSAMPLES;i++){
    if ( fscanf(fp,"%s %s %s %s %s %u", sx1,sy1,sx2,sy2,s3,&select)!= 6 ) {
      printf("Read error in file %s\n",DATA);
      return -1;
    }
    mpz_set_str(x,sx1,10);
    mpz_set_str(y,sy1,10);
    mpz_set_str(x1,sx1,10);
    mpz_set_str(y1,sy1,10);
    mpz_set_str(x2,sx2,10);
    mpz_set_str(y2,sy2,10);
    select ^= 1;
    
    mpz_init_set_ui(z,1);
    for(j=0; j<32; j++){
      eac[j] = 0;
      for (k=0; k<8; k++)
	eac[j] |= ((s3[8*j+k]-'0')<<k);
    }
    eac_end_256_smult(x1,y1,z,select,eac,p,beta);
    eac_to_affine(x1,y1,z,p);
    
    if (mpz_cmp(x1,x2)||mpz_cmp(y1,y2)){
      gmp_printf("[error sample %d] P = [%Zd, %Zd]\ngot [%Zd, %Zd]\nexpected [%Zd, %Zd]\n",i,x,y,x1,y1,x2,y2);
      break;
      }
  }
  printf("%d samples tested - no errors\n", NSAMPLES);
  mpz_clears(x,y,x1,y1,x2,y2,z,p,beta,NULL);
  fclose(fp);
} 
    
