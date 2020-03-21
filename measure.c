#define _GNU_SOURCE

#include <unistd.h>
#include <stdio.h>
#include <sys/syscall.h>
#include <inttypes.h>
#include <gmp.h>

#include "zaddu.h"
#include "eac_param_331.h"


#define NTEST 10000 // nombre de fois ou on repete le meme jeu de donnees
#define NSAMPLES 1000 // nombre differents de jeu de donnees


/**** Measurements procedures according to INTEL white paper

      "How to benchmark code execution times on INTEL IA-32 and IA-64" 
 
*****/

// ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Ne pas oublier de desactiver le turbo boost
// /bin/sh -c "/usr/bin/echo 1 > /sys/devices/system/cpu/intel_pstate/no_turbo
// pour fiabiliser la mesure

inline static uint64_t cpucyclesStart (void) {

  unsigned hi, lo;
  __asm__ __volatile__ (	"CPUID\n\t"
				"RDTSC\n\t"
				"mov %%edx, %0\n\t"
				"mov %%eax, %1\n\t"
				: "=r" (hi), "=r" (lo)
				:
				: "%rax", "%rbx", "%rcx", "%rdx");

  return ((uint64_t)lo)^(((uint64_t)hi)<<32);


}

inline static uint64_t cpucyclesStop (void) {

  unsigned hi, lo;
  __asm__ __volatile__(	"RDTSCP\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t"
			"CPUID\n\t"
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
				
  return ((uint64_t)lo)^(((uint64_t)hi)<<32);


}




int main() {
  
  unsigned long long timer , meanTimer =0, t1, t2;

  FILE * fp = fopen(DATA, "r");
  
  // On chauffe les caches
  mpz_t x1, y1, p, beta;

  char s1[FIELD_BASE_10_SIZE+1], s2[FIELD_BASE_10_SIZE+1], s3[SCALAR_BIT_SIZE+1]; 
  unsigned char eac[32];
  unsigned int select;
  int j,k;
  mpz_init_set_str(p,P_EAC,10);
  mpz_init_set_str(beta,B_EAC,10);

  mpz_init_set_str(x1, X1, 10);
  mpz_init_set_str(y1,Y1, 10);
  for (int i=0; i<32; i++) eac[i]=0;
  select = 1;
 
  
  for(int i=0;i<NTEST;i++)
    {
      // appel de la fonction a mesurer a mettre ici
      // juste pour chauffer les caches
      eac_end_256_smult(x1,y1,select,eac,p,beta);
    }
  
  // timing
  for(int i=0;i<NSAMPLES;i++)
    {

      // initialiser un nouveau jeu de donnees a tester
      if ( fscanf(fp,"%s %s %s %u", s1,s2, s3, &select)!= 4 ) {
	printf("Read error in file %s\n",DATA);
	return -1;
      }
      mpz_init_set_str(x1,s1,10);
      mpz_init_set_str(y1,s2,10);
      for(j=0; j<32; j++){
	eac[j] = 0;
	for (k=0; k<8; k++) eac[j] |= (s3[8*j+k]<<k);
      }
      timer = (unsigned long long int)0x1<<63;
      for(int j=0;j<NTEST;j++)
	{
	  t1 = cpucyclesStart();
	  // appeler la fonction ici avec toujours le meme jeu de donnees
	  eac_end_256_smult(x1,y1,select,eac,p,beta);
	  t2 = cpucyclesStop();
	  if(timer>t2-t1) timer = t2-t1;
	}
		
      meanTimer += timer;
    }
  
  printf("\n %lld CPU cycles", meanTimer/NSAMPLES);
  mpz_clears(x1,y1,p,beta,NULL);
  fclose(fp);
  return 0;
}



