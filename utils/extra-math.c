/* ------------------------------------------------------------------ *
 * Permission to copy and/or modify all or part of this work is       *
 *   granted, provided that this NO WARRANTY and this copyright	      *
 *   notice are retained verbatim and are displayed conspicuously. If *
 *   anyone needs other permissions that aren't covered by the above, *
 *   please contact the author					      *
 *      							      *
 * NO WARRANTY: THIS WORK IS PROVIDED ON AN "AS IS" BASIS.  THE       *
 *   AUTHOR PROVIDES NO WARRANTY WHATSOEVER, EITHER EXPRESS OR        *
 *   IMPLIED, REGARDING THE WORK, INCLUDING WARRANTIES WITH THE WORK, *
 *   INCLUDING WARRANTIES WITH RESPECT TO ITS MERCHANTABILITY OR      *
 *   FITNESS FOR ANY PARTICULAR PURPOSE.			      *
 *								      *
 * Philip Lijnzaad, lijnzaad@embl-heidelberg.de			      *
 * ------------------------------------------------------------------ */

/* some additional math stuff */

#include <signal.h>
#include <stdlib.h>
#include <string.h> 

#include "utils.h"
#include "vecmat.h"
#include "extra-math.h"
#include "jacobi.h"

#ifdef ultrix
#include <mips/fpu.h>
#endif 

#ifdef sgi 
#include <sys/fpu.h>
#endif

#if defined(ultrix) || defined(sgi)
/* provide decent floating point exception handling */

static void fpe_handler(int dummy);

void catch_fpe(int flags) {
  union fpc_csr fpc;
#define FPC fpc.fc_struct
  
  fpc.fc_word=0;			/* since it may contain garbage */
  if( flags & FPE_INEXACT) FPC.en_inexact =1;
  if( flags & FPE_UNDERFLOW) FPC.en_underflow =1;
  if( flags & FPE_OVERFLOW) FPC.en_overflow =1;
  if( flags & FPE_DIVIDE0) FPC.en_divide0 =1;
  if( flags & FPE_INVALID) FPC.en_invalid =1;

  set_fpc_csr(fpc.fc_word);		/* load bits into registers */

  signal(SIGFPE, fpe_handler);		/* install the handler */
}

static void fpe_handler(int dummy) {	/* should be void func(int) */
  union fpc_csr fpc;
#define FPC fpc.fc_struct

  fpc.fc_word=get_fpc_csr();

#define ERR(X) fprintf(stderr, "Floating point exception: %s\n",X)
  if( FPC.se_inexact && FPC.en_inexact ) { ERR("inexact result");
			 FPC.se_inexact=0; set_fpc_csr(fpc.fc_word); }
  if( FPC.se_underflow && FPC.en_underflow ) { ERR("underflow");
			 FPC.se_underflow=0; set_fpc_csr(fpc.fc_word); }
  if( FPC.se_overflow && FPC.en_overflow ) { ERR("overflow");
			 FPC.se_overflow=0; set_fpc_csr(fpc.fc_word); }
  if( FPC.se_divide0 && FPC.en_divide0 ) { ERR("division by zero");
			 FPC.se_divide0=0; set_fpc_csr(fpc.fc_word); }
  if( FPC.se_invalid && FPC.en_invalid ) { ERR("invalid operation");
			 FPC.se_invalid=0; set_fpc_csr(fpc.fc_word); }
  abort();
} /* fpe_handler */

#endif /* if defined ultrix || sgi */

#define FLOAT_TOL 1e-6
#define DOUBLE_TOL 1e-8

double safe_acos(const char *filename, int linenr, double f) {
  if ( (fabs(f) < 1.00) )
    return acos(f);
  if ( (fabs(f) - 1.00)  < DOUBLE_TOL ) { 
    warn("dubious argument %f passed to acos,",f);
    return 0.0; 
  }
  die("invalid argument %f passed to acos,",f);
  return -1.0;
}

double safe_asin(const char *filename, int linenr, double f) {
  if ( (fabs(f) < 1.00) )
    return( asin(f) );
  if ( (fabs(f) - 1.00)  < DOUBLE_TOL ) {
    warn("dubious argument %f passed to asin,",f);
    return(M_PI_2); 
  }
  die("invalid argument %f passed to asin,",f);
  return -1.0;
}


double safe_sqrt(const char *filename, int linenr, double f) {
  if ( f >= 0.00)
    return( sqrt(f) );
  die("invalid argument %f passed to sqrt,", f);
  return -1.0;
}

#ifdef SINGLE_PRECISION_AVAIL
float safe_facos(const char *filename, int linenr, float  f) {
  if ( (fabs(f) < 1.00) )
    return( acos(f) );
  if ( (fabs(f) - 1.00)  < FLOAT_TOL ) {
    warn("dubious argument %f passed to facos", f);
    return 0.0; 
  }
  die("invalid argument %f passed to facos", f);
  return -1.0;
}

float safe_fasin(const char *filename, int linenr, float f) {
  if ( (fabs(f) < 1.00) )
    return( asin(f) );
  if ( (fabs(f) - 1.00)  < FLOAT_TOL )
    return(M_PI_2); 
  die("invalid argument %f passed to fasin, at %s line %d\n",
      f, filename, linenr);
  return -1.0;
}

float safe_fsqrt(const char *filename, int linenr, float f) {
  if ( f >= 0.00)
    return( fsqrt(f) );
  die("invalid argument %f passed to fsqrt,", f); 
  return -1.0;
}
#endif /* SINGLE_PRECISION_AVAIL */

#define RETURN(S) { warn(S); return M_2PI; }
double dihedral(double *a, double *b, double *c, double *d) {
  /*  */
  int i;
  double ab[3], bc[3], cd[3], p[3], dot, lp, lbc, e;

  VEC3MIN(ab, b, a);
  VEC3MIN(bc, c, b);			/* normal of the plane to project on */
  VEC3MIN(cd, d, c);

  dot=DOT3P(ab, bc);
  lbc = DOT3P(bc, bc);
  if (lbc <= 0.0001 )
    RETURN("no dihedral: middle two points coincident\n");

  lp=0.0;
  dot /= lbc;
  for (i=2; i>=0; i--){			/* project onto plane */
    e=p[i]= dot*bc[i] - ab[i];		/* projection of ab */
    lp += SQR(e);
  }

  if (lp <= 0.0001)
    RETURN("no dihedral: first three points colinear\n");

  lp=sqrt(lp);				/* length of p */
  lbc=sqrt(lbc);			/* length of bc */
  for (i=2; i>=0; i--){			/* project onto plane */
    p[i] /= lp;				/* now local x-axis */
    bc[i] /= lbc;			/* now local z-axis */
  }
  CROSS3P(ab, bc, p);			/* ab now local y axis */

  lp = DOT3P(cd, p);			/* projection onto x-axis */
  lbc = DOT3P(cd, ab);			/* projection onto y-axis */

  if ( fabs(lp) <= 0.0001 && fabs(lbc) <= 0.0001 )
    RETURN("no dihedral: latter three points colinear\n");
  return atan2(lbc, lp);
} /* dihedral */

double fdihedral(float *a, float *b, float *c, float *d) {
  /*  */
  int i;
  double dot, lp, lbc, e;
  float ab[3], bc[3], cd[3], p[3];

  VEC3MIN(ab, b, a);
  VEC3MIN(bc, c, b);			/* normal of the plane to project on */
  VEC3MIN(cd, d, c);

  dot=DOT3P(ab, bc);
  lbc = DOT3P(bc, bc);
  if (lbc <= 0.0001 )
    RETURN("no dihedral: middle two points coincident\n");

  lp=0.0;
  dot /= lbc;
  for (i=2; i>=0; i--){			/* project onto plane */
    e=p[i]= dot*bc[i] - ab[i];		/* projection of ab */
    lp += SQR(e);
  }

  if (lp <= 0.0001)
    RETURN("no dihedral: first three points colinear\n");

  lp=sqrt(lp);				/* length of p */
  lbc=sqrt(lbc);			/* length of bc */
  for (i=2; i>=0; i--){			/* project onto plane */
    p[i] /= lp;				/* now local x-axis */
    bc[i] /= lbc;			/* now local z-axis */
  }
  CROSS3P(ab, bc, p);			/* ab now local y axis */

  lp = DOT3P(cd, p);			/* projection onto x-axis */
  lbc = DOT3P(cd, ab);			/* projection onto y-axis */

  if ( fabs(lp) <= 0.0001 && fabs(lbc) <= 0.0001 )
    RETURN("no dihedral: latter three points colinear\n");
  return atan2(lbc, lp);
} /* fdihedral */

#undef RETURN

int root2(double a, double b, double c, double *sol) {
  /* solve quadratic at^2 + bt + c = 0, and write roots to sol[],
   *   sorted. returns nr of roots */
  double det=SQR(b)- 4.0*a*c, q; 

  if (det <= - DOUBLE_TOL)
    return(0);
  if ( fabs(det) < DOUBLE_TOL ) {		/* call it zero */
    sol[0]= -0.5*b/a;
    return(1); 
  }
  if (b<0) 
    q= -0.5*(b-sqrt(det)); 
  else
    q= -0.5*(b+sqrt(det));

  sol[0]=q/a;
  sol[1]=c/q;

  if ( sol[0] > sol[1] )		/* sort it */
    dSWAP(sol[0], sol[1]); 

  return(2); 
} /* root2 */

int root3(double x, double a, double b, double c, double *sol) {
  /* find solution of xt^3 + at^2 + bt +c; write roots (sorted) to sol[],
   *   and return nr of roots. */
  double q, r, q3, r2, d, theta, f, a_3;
  int i,j;

  if (fabs(x) < DOUBLE_TOL)			/* not a cubic equation */
    return( root2(a, b, c, sol));
  
  if (fabs(c) < DOUBLE_TOL) {			/* x*(quadratic equation) */
    i=root2(x, a, b, sol);
    for (j=0; j<i; j++) 
      if ( fabs(sol[j]) < DOUBLE_TOL )		/* x^2*(linear equation): */
	return(i);			/* already 0 solution present */

    sol[i]=0.0;
    if (i==0)
      return(1);
    
    qsort(sol, 3, sizeof(double), double_cmp);
    return(3);
  }

  a /=x;				/* scale to x^3 + ax^2 etc. */
  b /=x;
  c /=x;
  a_3= a/3.0;


  q=( SQR(a)-3.0*b )/9.0;
  r=( 2.0*CUBE(a) - 9.0*a*b + 27.0*c )/54.0;

  if ( fabs(q) < DOUBLE_TOL && fabs(r) < DOUBLE_TOL) { /* of form (x-a)^3: simple */
    sol[0]= a/3.0; 
    return(1); 
  }
  
  q3=CUBE(q);
  r2=SQR(r); 

  d=q3-r2;

  if ( d >= -DOUBLE_TOL) {			/* casus irreducibilis */
    theta=acos(r/sqrt(q3));
    f= -2.0*sqrt(q);

    if ( fabs(theta) < DOUBLE_TOL ) {		/* one extreme 'touches' x axis  */
      sol[0]= f - a_3; 
      sol[1]= -0.5*f -a_3;
      if (sol[0] > sol[1])
	dSWAP(sol[0], sol[1]); 
      return(2);			/* since next solution is the same */
    }

    sol[0]= f*cos(theta/3.0) - a_3;
    sol[1]= f*cos( (theta + M_2PI) /3.0) - a_3;
    sol[2]= f*cos( (theta + M_4PI) /3.0) - a_3;
    qsort(sol, 3, sizeof(double), double_cmp);
    return(3);
  }

  if ( d < -DOUBLE_TOL ) {
    f=sqrt( -d ); 
    
    if (r > 0.0)  {
      f=pow(f + r, 0.33333333333333);
      sol[0]= - ( f + q/f ) - a_3; 
    } else { 
      f=cbrt(f - r);
      sol[0]= pow( f + q/f, 0.33333333333333 ) - a_3; 
    }
    return(1); 
  }
  return -1;
} /* root3 */


static double scale_ellipsoid(int n, float ** xyz, float * center, 
			   double fraction, float ** eigvecs, double *axes) {
  /*
   * does    : adjusts the scale of AXES such that FRACTION part of
   *   the N points XYZ lie within the ellipsoid defined by them 
   *
   * gets    : nr of points, coordinates, 
   *
   * affects : axes[0,1,2]
   *
   * returns : the scale factor
   *
   * warns   : 
   *
   * comment : 
   */
  int i,j, k;
  float d, d2, *dists, scale, u[3];
  
  dists = (float*)ALLOCA(n*sizeof(float));

  LOOPDOWN (j, n) {			/* set up scaled dists. from center */
    d2=0.0;
    VEC3MIN(u, xyz[j], center);
    LOOPDOWN (i,3) {
      d=0.0;
      LOOPDOWN(k,3) 
	d += u[k] * eigvecs[k][i];	/* eigvecs are COLUMNS, not rows !!! */
      d/=axes[i];			/* scale to unit sphere */
      d2 += SQR(d);
    }
    dists[j]=d2;
  }
  qsort(dists, n, sizeof(float), float_cmp);
  k=(int)(fraction*n);
  scale = sqrt(dists[k]);
  LOOPDOWN (i,3)
    axes[i] *= scale;
  AFREEA(dists);
  return scale;				/* not of much use ... */
} /* scale_ellipsoid */

void print_ellipsoid(double factor, float *center, 
			    double * axes, float **eigvecs) {
  /* for debbugging purposes */
  int i,j,k;
  double r,z, x[3], y[3], phi;
# include "vecmat.h"
# define N 10

#define FILENAME  "el"
  FILE * file =fopen(FILENAME, "w");
  if (!file)
    perror(FILENAME),exit(1);

  /* construct unit sphere around origin, scale to ellipsoid with semiaxes
     along xyz, do an additional scaling by factor, rotate, and translate: */
  for (i=0; i<N; i++) {			/* 'circles of latitude' */
    z= 2.0* i/(double)N -1;
    r= sqrt(1.00001 - SQR(z));		/* radius of this circle */
    for (j=0; j<=N; j++) {
      x[0]= z * axes[0];
      phi= M_2PI*j/(double)N;
      x[1]= r*sin(phi)*axes[1];
      x[2]= r*cos(phi)*axes[2];
      MAT3VEC(y, eigvecs, x);		/* rotate to proper orientation */
      /* (don't transpose, as the eig. vec.'s are columns of eigvecs !!! */
      /*	VEC3MAT(y, x, eigvecs); */
      LOOPDOWN(k,3)
	y[k]=factor*y[k] + center[k];
      fprintf(file, "%s %7.3f %7.3f %7.3f\n", j==0 ? ".m" : ".d",
	      y[0], y[1], y[2]);
    }
  }
  fclose(file);
  return;
} /* print_ellipsoid */

double xydev(int n, float **xyz, float *center, 
		    float ** eigvecs, double *axes) { 
  int i,j;
  float *x;
  double d, d2;

  d2=0.0;
  LOOPDOWN(j,n) { 
    x=xyz[j];
    d=0.0;
    LOOPDOWN(i, 3)			/* projection onto local z-axis */
      d += (x[i]-center[i])*eigvecs[i][2];
    d *= axes[2];			/* 'cos eigvecs are unit vecs */
    d2 += SQR(d); 
  } /* j */
  
  return sqrt(d2/(n-1));
} /* xydev */

int ellipsaxes_3D(int n, float ** xyz, double fraction, double *axes) { 
  /*
   * does    : calculates the axes of the EQUIVALENT (see Taylor
   *   &Thornton, J.Mol.Graph. ???) ellipsoid that contains FRACTION points
   *   of the N coordinates XYZ.N
   *
   * gets    : nr of points, coordinates, fraction of points within ellipsoid
   *
   * affects : axes[*] is set to the respective axes, sorted UPWARD
   *
   * returns : 1 if couldn't fit ellips, 0 otherwise.
   *
   * warns   : if axis lengths are nearly zero
   *
   * comment : 
   */
  
  int i, j, k, err;
  double div, factor;
  float *ptr, x;
  float *mat[3], eigvals[3], *eigvecs[3], center[3];

  for (i=2; i>=0; i--) {
    mat[i]=ALLOCA(sizeof(float[4]));
    bzero((char*)mat[i], sizeof(float[4]));
    eigvecs[i]=ALLOCA(sizeof(float[3])); /* eigenvecs are columns, not rows */
  }

  for (k= n-1; k>=0; k--) {
    ptr= xyz[k];
    for (i=2; i>=0; i--)  { 
      x=ptr[i];
      mat[i][3] += x;			/* keep average */
      for (j=2; j>=0; j--)
	mat[i][j] += x*ptr[j];
    } /* for i */
  }
  div= 1.0/(float)(n);
  for (i=2; i>=0; i--) { 
    center[i]=mat[i][3]*div;		/* calc center of points */
    for (j=2; j>=0; j--)
      mat[i][j] -= (mat[i][3]*mat[j][3]*div);
  }

/* eigenvalues of symmetric 3D matrix can be found directly from
 * characteristic polynomial:
 * 
 * det( [ [ a - x,   d,     e   ]
 *        [   d,   b - x,   f   ]
 *        [   e,     f,   c - x ] ] ) = -1 *
 * x^3 - (a + b + c) x^2 + x*(b c + a c + a b - e^2 - f^2 - d^2) + c d^2 
 *       + a f^2 + b e^2 - a b c - 2 d e f == 0
 */    
/* 
 *   a=mat[0][0];
 *   b=mat[1][1];
 *   c=mat[2][2];
 *   d=mat[1][0]; /# assert(fabs(d-mat[1][0] < 0.001)); #/
 *   e=mat[2][0]; /# assert(fabs(e-mat[2][0] < 0.001)); #/
 *   f=mat[2][1]; /# assert(fabs(f-mat[2][1] < 0.001)); #/
 *   
 *   d2=SQR(d);
 *   e2=SQR(e);
 *   f2=SQR(f);
 *   rg = a+b+c;				/# radius of gyration #/
 *   nroots=root3(1.0,			/# x^3 #/
 * 	       -rg,			/# x^2 #/
 * 	       (a*b + b*c + c*a - (d2 + e2 + f2)), /# x #/
 * 	       (c*d2 + a*f2 + b*e2 -(a*b*c + 2*d*e*f)),	/# - #/
 * 	       roots);
 *   if (nroots != 3) { 
 *     warn("found only %d eigenvalues\n", nroots);
 *     return 1;
 *   }
 */

  jacobi(3, mat, eigvals, eigvecs);
  
  err=0;
  for (i=2; i>=0; i--) {			/* convert to float */
/*    axes[i]=sqrt(roots[i]); */
    axes[i]=sqrt(eigvals[i]);
    err += (axes[i] <= 0.0001);
  }
  if(err)
    return warn("no scaling done; some axes too small");

  factor=scale_ellipsoid(n, xyz, center, fraction - 0.0000001, eigvecs, axes);

/*  print_ellipsoid(1.0, center, axes, eigvecs);
 * To make a pdb file out of this, try:
 * awk '{printf "ATOM  %5d  CX  UNK Z%4d     %7.3f %7.3f %7.3f\n", 
 *           NR, NR/10, $2, $3, $4}'
 *
 * printf("XYDEV: %.3f %.3f\n", axes[2], xydev(n, xyz, center, eigvecs, axes) );
 */

  for (i=2; i>=0; i--) {
    AFREEA(mat[i]);
    AFREEA(eigvecs[i]);
  }
  return 0;

} /* ellipsaxes */

/* eigenvalues of symmetric 3D matrix can be found directly from
 * characteristic polynomial:
 * 
 * det( [ [ a - x,   d,     e   ]
 *        [   d,   b - x,   f   ]
 *        [   e,     f,   c - x ] ] ) = -1 *
 * x^3 - (a + b + c) x^2 + x*(b c + a c + a b - e^2 - f^2 - d^2) + c d^2 
 *       + a f^2 + b e^2 - a b c - 2 d e f == 0
 */    

/* int root4(double a, double b, double c, double d, double e, double *sol) {
 *   /# doesn't work yet !!!! #/
 *   double s, sol2[3], a1, a2, a3, r1, r2;
 *   int i,j,n; 
 * 
 *   ABORT("routine root4() doesn't work yet !!!!"); 
 * 
 *   if ( fabs(a) < DOUBLE_TOL )			/# not a quartic equation #/
 *     return( root3(b, c, d, e, sol) );
 *   if (fabs (e)  < DOUBLE_TOL )  {		/# x*(ax^3 +bx^2 +cx +d) #/
 *     i = root3(a, b, c, d, sol);
 *     for (j=0; j<i; j++) 
 *       if (fabs(sol[j]) < DOUBLE_TOL)
 * 	return(i);			/# 0-solution already present #/
 * 
 *     sol[i]=0.0;				/# add 0-solution #/
 *     qsort(sol, i+1, sizeof(double), double_cmp);
 *     return(i+1);
 *   }
 * 
 *   i=root3( 1.0,  -c, (b*d - 4*a*e), (a*SQR(d) + SQR(b)*e - 4*a*c*e), sol2); 
 *   s= sol2[i-1];				/# the most positive solution #/
 *   
 *   r1= SQR(b) - 4*a*(c -s);
 *   r2= SQR(s) - 4*a*e;
 *   
 *   if (r1 < -DOUBLE_TOL || r2 < -DOUBLE_TOL)		/# no real roots #/
 *     return(0); 
 * 
 *   n=0;
 *   if (r1 < DOUBLE_TOL)				/# consider it zero #/
 *     r1=0.0, n++;
 *   else
 *     r1=sqrt(r1); 
 * 
 *   if (r2 < DOUBLE_TOL)				/# consider it zero #/
 *     r2=0.0, n++;
 *   else
 *     r2=sqrt(r2);
 * 
 *   c = b*s - 2*a*d; 
 *   
 *   if (n==2) {
 *     i=root2( 2*a, b, s, sol );
 *     return(i); 
 *   }
 *     
 *   if ( c * r1 * r2 < 0.0) {
 *     i=root2( 2*a, b - r1, s + r2, sol );
 *     j=root2( 2*a, b + r1, s - r2, sol+i );
 *   } else { 
 *     i=root2( 2*a, b - r1, s - r2, sol );
 *     j=root2( 2*a, b + r1, s + r2, sol+i );
 *   }    
 *   n=i+j;				/# total nr of roots found #/
 *   qsort(sol, n, sizeof(double), double_cmp); /# no double roots, i think ! #/
 *   
 *   return(n); 
 * } /# root4 #/
 */
