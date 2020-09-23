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

/* some additional math stuff. To speed up, get rid of checking code by */
/* #define'ing NOSAFE */

#ifndef _EXTRA_MATH_H_

#define _EXTRA_MATH_H_

#include <math.h>
#include "utils.h"

/* float constants */

#ifndef PI
#  define PI          M_PI		/* #defined in math.h */
#endif

#ifndef M_2PI
#  define M_2PI  6.28318530717958647692 /* two-pi */
#endif

#ifndef M_4PI
#  define M_4PI 12.56637061435917295384   /* four pi */
#endif

#ifndef M_1_2PI
#  define M_1_2PI 0.15915494309189533576	/* two-pi reciprocal */
#endif

#ifndef SIN60
#  define SIN60 0.86602541		/* sin(PI/3) */
#endif

#ifndef RSIN60
#  define RSIN60 1.15470053		/* 1.0/sin(PI/3) */
#endif

#ifndef DEGPERRAD
#  define DEGPERRAD 57.295796
#endif

#ifndef RADPERDEG
#  define RADPERDEG 0.017453293
#endif

#ifndef TODEGR
#  define TODEGR(A)    ((A)*180.0*(M_1_PI)) /* likewise */
#endif

#ifndef TORAD
#  define TORAD(A)     ((A)*(PI)/180.0)
#endif

#ifndef M_GOLRAT
#  define M_GOLRAT     0.618033988749895
#endif

#ifndef M_1_GOLRAT
#  define M_1_GOLRAT   1.61803398874989
#endif

#ifndef M_SQRT2
#  define M_SQRT2	     1.414213562373095
#endif

#ifndef M_1SQRT2
#  define M_1SQRT2     0.707106781186548
#endif

#ifndef M_SQRT3
#  define M_SQRT3	 1.7320508075688772
#endif

#ifndef M_1SQRT3
#  define M_1SQRT3 0.57735026918962573
#endif



/* single precision versions of math functions: NOT IN math.h  @#??#?? */

#ifdef SINGLE_PRECISION_AVAIL

float fsqrt(float);
float fpow(float);

float flog(float);
float flog10(float);
float flog1p(float);
float fexp(float);
float fexpm1(float);

float fsin(float);
float fcos(float);
float ftan(float);

float fasin(float);
float facos(float);
float fatan(float);
float fatan2(float, float y);

float fceil(float);
float ffloor(float);
float ftrunc(float);
double trunc(double x);			/* forgotten */

float fsinh(float);
float fcosh(float);
float ftanh(float);

float safe_fsqrt(const char *filename, int linenr, float f);
float safe_facos(const char *filename, int linenr, float);
float safe_fasin(const char *filename, int linenr, float);

#endif /* SINGLE_PRECISION_AVAIL */
/* double cbrt(double x);	        DONT'T USE:  BUGGED !!!! */


#ifndef NOSAFE

#define SQRT(X) safe_sqrt(__FILE__, __LINE__, (X))
#define FSQRT(X) safe_fsqrt(__FILE__, __LINE__, (X))

#define ACOS(X) safe_acos(__FILE__, __LINE__, (X))
#define ASIN(X) safe_asin(__FILE__, __LINE__, (X))
#define FACOS(X) safe_facos(__FILE__, __LINE__, (X))
#define FASIN(X) safe_fasin(__FILE__, __LINE__, (X))

#else  /* if NOSAFE */

#define SQRT sqrt
#define FSQRT fsqrt
#define ACOS acos
#define ASIN asin
#define FACOS facos
#define FASIN fasin

#endif /* NOSAFE */


double safe_sqrt(const char *filename, int linenr, double f);
double safe_acos(const char *filename, int linenr, double x);
double safe_asin(const char *filename, int linenr, double x);

#if  defined(__alpha__) || defined(__osf__) /* don't need it */
#define catch_fpe(DUMMY)		/* mask calls to this function */
#define FPE_INEXACT     0		/* and also use of these constants */
#define FPE_UNDERFLOW   0
#define FPE_OVERFLOW    0
#define FPE_DIVIDE0     0
#define FPE_INVALID     0
#define FPE_DEFAULT	0
#else

void catch_fpe(int fpe_flags);

/* The ultrix and sgi mips machines leave the handling of IEEE floating */
/* point exceptions (eg. sqrt(-1.0), 1.0/0.0) to the programmer. To force */
/* program termination upon occurrence of such an exception, include a */
/* call to the above function, and specify the exceptions to catch by */
/* one or more of the following flags; often catch_fpe(FPE_DEFAULT) will do */

#define FPE_INEXACT    01
#define FPE_UNDERFLOW  02
#define FPE_OVERFLOW   04
#define FPE_DIVIDE0   010
#define FPE_INVALID   020
#define FPE_DEFAULT (FPE_OVERFLOW | FPE_DIVIDE0 | FPE_INVALID | FPE_UNDERFLOW)
/* also underflow, since that's usually the result of a bug; seldom will */
/* you let an algorithm converge by relying on underflow */

#endif /* #ifdef __alpha__ */


/* signal raised after handling the exception */
/* #define CATCH_FPE_SIG DIE_SIG */
/* now just DIE_SIG */

/* float fcabs(struct { float x, float y } z ); 
 * float fhypot(float, float);
 */

/* miscellaneous: */

#define MATRIX(TYPE,mat, nrow, ncol) /* unlike NumRec. matrix() */\
{ int i; mat = (TYPE **) MALLOC(nrow*sizeof(TYPE*));\
    for(i=0;i<nrow;i++)mat[i]=(TYPE *) MALLOC(ncol*sizeof(TYPE));}

#define MATRIX_NZ(TYPE,mat,rl,ru,cl,cu) /* NON-ZERO OFFSET matrix */\
{ int i; mat = (TYPE **) MALLOC((ru-rl+1)*sizeof(TYPE*)) - rl;\
    for(i=rl;i<=ru;i++)mat[i]=(TYPE *) MALLOC((cu-cl+1)*sizeof(TYPE)) -cl;}


int root2(double a, double b, double c, double *sol); /* roots of */
  /* quadratic equation; returns nr of solutions, and write into *sol */
int root3(double a, double b, double c, double d, double *sol);
					/* likewise for cubic equation */
/* int root4(double a, double b, double c, double d, double e, double *sol); 
 *  doesn't work yet			/# likewise for quartic equation #/
 */
int ellipsaxes_3D(int n, float ** xyz, double fraction, double *lengths);
/* calcs. LENGTHS of ellipsoid that contains FRACTION of the N points XYZ */
double dihedral(double *a, double *b, double *c, double *d);
double fdihedral(float *a, float *b, float *c, float *d);


#endif /* #ifndef _EXTRA_MATH_H_ */
