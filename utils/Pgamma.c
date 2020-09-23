#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
/* #include "extra-math.h"  */
#include "utils.h"
#include "alloc.h" 

#define LINLEN 1024

/* function to read pairs of a (shape) and x from stdin, and print
   incomplete gamma to stdout, following Num. Rec. in C. It relies on 
   math.h / libm.a defining the log(gamma(x)) function a lgamma
*/

#define ITMAX 200
#define EPS 1.0e-7

/* incomplete gamma function using series represention, for x > a+1 */
void gser(float* gamser, float a, float x, float*gln) {
  int n;
  float sum, del, ap;
  *gln=(float)lgamma(a);

  if (x <= 0.0) {
    if (x< 0.0)
      DIE("x < 0  in GSER,");
    *gamser=0.0;
    return;
  } else { 
    ap=a;
    del=sum=1.0/a;
    for(n=1; n<= ITMAX; n++) {
      ap += 1.0;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS ) {
        *gamser= sum*exp( -x+a*log(x)-(*gln) );
        return;
      }
    }
    DIE("a too loarge, or ITMAX too small in GSER");

  }
} /* gser */

/* incomplete gamma function using continued fractions, for x < a+1 */
void gcf(float*gammcf, float a, float x, float*gln) { 
  int n;
  float gold=0.0, g, fac=1.0, b1=1.0;
  float b0=0.0, anf, ana, an, a1, a0=1.0;
  
  *gln=(float)lgamma(a);
  a1=x;
  for (n=1; n<=ITMAX; n++) {
    an=(float)n;
    ana=an-a;
    a0=(a1+a0*ana)*fac;
    b0=(b1+b0*ana)*fac;
    anf=an*fac;
    a1=x*a0+anf*a1;
    b1=x*b0+anf*b1;
    if (a1) {
      fac=1.0/a1;
      g=b1* fac;
      if ( fabs((g - gold)/g) < EPS)  {
        float e = exp(-x + a*log(x)-(*gln))*g;
        *gammcf = e;
        return;
      }
      gold=g;
    }
  }
  DIE("a too large, or ITMAX too small in GCF,");
} /* gcf */


/* incomplete gamma distribution, for shape factor a */
float gammap(float a, float x) { 
  float gamser, gammcf, gln;
  if (x < 0.0 || a <= 0.0) 
    DIE("Invalid args to gammap,");
  if ( x < (a+1.0)) {
    gser(&gamser, a, x, &gln);
    return gamser;
  } else { 
    gcf(&gammcf, a, x, &gln);
    return 1.0 - gammcf;
  }
}


int main(int argc, char**argv) {
  char line[LINLEN];
  /*  char function[LINLEN]; */

  float x, a;

  while(fgets(line, LINLEN, stdin)) {
    if (strlen(line) > LINLEN)
      DIE("line too long,");

    if(sscanf(line, "%g %g", &a, &x ) != 2)
      DIE("exspected <a> <x>,");
    

    printf("%12.8g\n", gammap(a, x));
  }
  return 0;
}






