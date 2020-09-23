#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

/* #include "nrutil.h" */

#include "utils.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

/* modified from Numerical Recipes */
#define TOL 0.0000001			/* not used; rely on machine precision */
#define NMAX 50				/* max nr of hyperplane rotations  */

static void eigsrt(int n, float *d, float **v) {
  /* sort eigenvalues and -vectors to descending order */
  int k,j,i;
  float p;

  for (i=0;i<n-1; i++) {
    p=d[k=i];
    for (j=i+1; j<n; j++)
      if (d[j] >= p) 
	p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=n-1; j>=0; j--) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
} /* eigsrt */

int jacobi(int n, float **a, float *d, float **v) {
  /*  computes eigenvalues and -vectors from matrix a of dimension
   *  n. Elements of a above diagonal are destroyed. Eigenvectors are
   *  returned in array d, and eigenvectors are the columns of
   *  v. Eigenvalues and -vectors are sorted to descending order. The
   *  number of rotations needed is returned 
   */ 
  int j,iq,ip,i, nrot;
  float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  /* 	b=vector(1,n);
   * 	z=vector(1,n);
   */
  b=ALLOCA(n*sizeof(float));	
  z=ALLOCA(n*sizeof(float));

  for (ip= n-1; ip>=0; ip--) {		/* build identity matrix for v*/
    for (iq= n-1; iq>=0; iq--) 
      v[ip][iq]=0.0;
    v[ip][ip]=1.0;
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }

  nrot=0;
  for (i=0; i < NMAX; i++) {
    sm=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++)
	sm += fabs(a[ip][iq]);
    }
/*     if (sm <= TOL ) { */
    if (sm == 0.0 ) {			/* not very decent ... */
      /* 			free_vector(z,1,n);
       * 			free_vector(b,1,n);
       */
      AFREEA(z);
      AFREEA(b);
      eigsrt(n, d, v);
      return nrot;			/* only exit point */
    }
    if (i < 3)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    
    for (ip = n-2; ip>=0; ip--) {	/* process all above-diagonal elts. */
      for (iq=ip+1; iq < n; iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 3			/* skip rotation if elt. small */
	    && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
	    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]) )
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((float)(fabs(h)+g) == (float)fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;		/* the elt. p,q that was zeroed */
	  for (j=ip-1; j>=0; j--) { ROTATE(a,j,ip,j,iq) } /* 0 < j < p */
	  for (j=ip+1; j<iq; j++) { ROTATE(a,ip,j,j,iq) } /* p < j < q */
	  for (j=iq+1; j<n; j++)  { ROTATE(a,ip,j,iq,j) } /*  q < j < n */
	  for (j=n-1; j>=0; j--)  { ROTATE(v,j,ip,j,iq) }
	  nrot++;
	} /* if i > 3 */
      }	/* for iq */
    } /* for ip */
    for (ip= n-1; ip>=0; ip--) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  } /* for i */
  die("Too many iterations in routine jacobi\n");
  return 0;
} /* jacobi */
#undef ROTATE
