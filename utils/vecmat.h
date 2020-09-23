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

/* vector and matrix stuff */

#include "math.h"

/* initializing */ 
#define  VEC3INIT(A,B) (A)[0]=(A)[1]=(A)[2]=(B)
#define  VEC4INIT(A,B) (A)[0]=(A)[1]=(A)[2]=(A)[3]=(B)
#define VECnINIT(N,A,B) { int i; for (i=0; i<(N); i++) A[i]=(B); }

#define MAT3INIT(A,B){ VEC3INIT((A)[0],B);\
			 VEC3INIT((A)[1],B);VEC3INIT((A)[2],B);}
#define MAT4INIT(A,B){ VEC4INIT((A)[0],B);VEC4INIT((A)[1],B);\
		       VEC4INIT((A)[2],B); VEC4INIT((A)[3],B);}

#define KR_DELTA(I,J) (I++J?1:0)	/* Kronecker delta */

#define ID3INIT(A) { /* 3D-identity matrix */\
		      (A)[1][1]=(A)[2][2]=(A)[0][0]=1.00;\
                    (A)[1][2]=(A)[1][0]=(A)[2][0]=(A)[2][1]=(A)[0][1]=\
		      (A)[0][2]=0.0;}

#define ID4INIT(A) IDnINIT(4,A)		/* should suffice */

#define IDnINIT(N,A)   { /* n-D identity matrix */\
		     int i, j;\
		     for (i=0; i<N; i++)\
		       for (j=0; j<N; j++)\
			   (A)[i][j]=(i==j)?1.00:0.00; }
/* assigning */
#define  VEC3ASS(A,B)  { (A)[0]=(B)[0]; (A)[1]=(B)[1]; (A)[2]=(B)[2]; }
#define  VEC4ASS(A,B)  { (A)[0]=(B)[0]; (A)[1]=(B)[1]; \
			   (A)[2]=(B)[2]; (A)[3]=(B)[3];}
#define VECnASS(N,A,B)  { int i; for (i=0; i<(N); i++) A[i]=B[i]; }

#define MAT3ASS(A,B)  { VEC3ASS((A)[0],(A)[0]);\
			VEC3ASS((A)[1],(B)[1]);\
		        VEC3ASS((A)[2],(B)[2]); }

#define MAT4ASS(A,B)  { VEC4ASS((A)[0],(A)[0]);\
			VEC4ASS((A)[1],(B)[1]);\
		        VEC4ASS((A)[2],(B)[2]);\
		        VEC4ASS((A)[3],(B)[3]); }
#define MATnASS(N,A,B) { int i; for (i=0; i<N; i++) VECnASS(N, A,B); }

/* vector operations yielding scalar */
#define DOT3P(A,B)    ( (A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2] )

#define NORM3_2(A)     (DOT3P((A),(A)))

#define COS(A,B)   ( DOT3P(A,B)/(fsqrt(DOT3P(A,A) * DOT3P(B,B))))

#define ANGLE(A,B)   (ACOS(DOT3P(A,B)/(sqrt(DOT3P(A,A) * DOT3P(B,B)))))

#define DIST3_2(A,B)   ( ((A)[1]-(B)[1])*((A)[1]-(B)[1])\
                          +((A)[2]-(B)[2])*((A)[2]-(B)[2])+\
			((A)[0]-(B)[0])*((A)[0]-(B)[0]))
/* ^^^ this is for obtaining distances, not for finding dist's < something ! */
/* finding distances smaller than or bigger than can be speeded up by */
/* searching structured (i.e. first look if center of residue is close */
/* enough), then by in turn considering 
	distance12_squared = 0.0;
	for (k = 1; k <=3; k++) {
		distance12 = *atm2_center++ - *atm1_center++;
		if (distance12 >= radius_sum) return (0);
		if (distance12 <=  (-radius_sum)) return (0);
		distance12_squared += distance12 * distance12;
	} additional test on distance12_squared every k ?
	This code was nicked from Michael Connolly 
*/

#define _DIFF(A,B) ((A)<(B)?((B)-(A)):((A)-(B))) /* no function calls !! */
#define _DISTCOMPARE(X,Y, RELOP, D) \
((_DIFF((X)[1],(Y)[1]) RELOP D)?\
 1:\
 ((_DIFF((X)[2],(Y)[2]) RELOP D)?\
  1:\
  ((_DIFF((X)[0],(Y)[0]) RELOP D )?\
   1:\
   (( SQR((X)[1]-(Y)[1]) + SQR((X)[2]-(Y)[2]) RELOP SQR(D))?\
    1:\
    ((DIST3_2(X,Y) RELOP SQR(D))?\
     1:0)))))

/* in fact this should be a C++ inline function !! Hopefully the optimizer */
/* will recognize the part common to 4th & 5th clause ... But apart from */
/* this i think this is formally provable the fastest way to get the */
/* answer. But of course this depends on the likelyhoods of hits over */
/* misses; if that is large, this may be slower */

#define DIST_GT(X,Y,D) (_DISTCOMPARE(X,Y, >, D))
#define DIST_GE(X,Y,D) (_DISTCOMPARE(X,Y, >=, D))
#define DIST_LT(X,Y,D) (!(_DISTCOMPARE(X,Y, >=, D)))
#define DIST_LE(X,Y,D) (!(_DISTCOMPARE(X,Y, >, D)))

/* matrix operations yielding scalar */

#define VEC3DET(A,B, C) ( (A)[0]*((B)[1]*(C)[2]-(B)[2]*(C)[1]) \
                -(B)[0]*((A)[1]*(C)[2]-(C)[1]*(A)[2]) \
	        +(C)[0]*((A)[1]*(B)[2]-(B)[1]*(A)[2]) )

#define DET3(A) (VEC3DET( (A)[0], (A)[1], (A)[2]))


#define VECnTIM(N, A, B) { int i; for (i=0; i<(N); i++)  A[i]*=(1.0*B); }

#define  VEC3DIV(A, B) { float rec=(1/(1.00*(B)));  VEC3TIM(A, rec); }

#define VECnDIV(N, A, B) { float rec=(1/(1.0*(B))); VECnTIM(N, A, rec); }

#define  VEC3NEG(A)    { (A)[1]= -(A)[1]; (A)[2]= -(A)[2]; (A)[0]= -(A)[0];  }

#define NRMIZE3(A) { double r; r=sqrt(DOT3P(A,A)); r=1/r; VEC3TIM(A, r); }


/* vector-scalar operations yielding vector */
#define  VEC3TIM(A, B) do {(A)[0]*=(1.0*B); (A)[1]*=(1.0*B);\
			     (A)[2]*=(1.0*B); } while(0)

#define  VEC3TIMASS(C, B, A) do {(C)[0]=(A)[0]*(B); (C)[1]=(A)[1]*(B);\
			     (C)[2]=(A)[2]*(B); } while(0)
/* vector-scalar operations yielding vector, plus assignment */

#define VECnTIMASS(N, A, B, C) { int i; for (i=0; i<(N); i++)\
				   A[i]=B[i]*1.0*(C); }

#define  VEC3DIVASS(A, B, C) { float rec=(1/(1.00*(C))); \
				 VEC3TIMASS(A, B, rec);}

#define VECnDIVASS(N, A, B, C) { float rec=(1/(1.0*(B)));\
  VECnTIMASS(N, A,B, rec); }

/* vector-vector operations yielding vector */
#define  VEC3INC(A, B)    { (A)[0]+=(B)[0]; \
			    (A)[1]+=(B)[1]; \
			    (A)[2]+=(B)[2]; }

#define  VEC3INCTIM(A, d, B)    { (A)[0] += (d)*(B)[0]; \
				  (A)[1] += (d)*(B)[1]; \
				  (A)[2] += (d)*(B)[2]; }

#define  VEC3DEC(A, B)    { (A)[0] -= (B)[0]; \
			    (A)[1] -= (B)[1]; \
			    (A)[2] -= (B)[2]; }

#define  VEC3DECTIM(A, d, B)    { (A)[0] -= (d)*(B)[0]; \
				  (A)[1] -= (d)*(B)[1]; \
				  (A)[2] -= (d)*(B)[2]; }

/* vector-vector operations yielding vector, plus assignment */
#define  VEC3PLUS(A, B, C) { (A)[0] = (B)[0] + (C)[0];\
			     (A)[1] = (B)[1] + (C)[1];\
			     (A)[2] = (B)[2] + (C)[2]; }

#define  VEC3PLUSTIM(A, B, d, C) { (A)[0] = (B)[0] + (d)*(C)[0];\
				   (A)[1] = (B)[1] + (d)*(C)[1];\
				   (A)[2] = (B)[2] + (d)*(C)[2];}

#define  VEC3MIN(A, B, C) { (A)[0] = (B)[0] - (C)[0];\
			    (A)[1] = (B)[1] - (C)[1];\
			    (A)[2] = (B)[2] - (C)[2]; } 

#define  VEC3MINTIM(A, B, d, C) { (A)[0] = (B)[0] - (d)*(C)[0];\
				  (A)[1] = (B)[1] - (d)*(C)[1];\
				  (A)[2] = (B)[2] - (d)*(C)[2];}

#define CROSS3P(A,B,C) { (A)[0]=(B)[1]*(C)[2]-(B)[2]*(C)[1]; \
                        (A)[1]=(B)[2]*(C)[0]-(B)[0]*(C)[2]; \
                        (A)[2]=(B)[0]*(C)[1]-(B)[1]*(C)[0]; }

/* scalar matrix operations */

#define MAT3TIM(A,B) { VECTIM((A)[1],B); VECTIM((A)[2], B); VECTIM((A)[0],B); }

#define MAT3DIV(A,B) { MAT3TIM(A, (1/B)) }

#define MAT3NEG(A)   { VEC3NEG((A)[1]); VEC3NEG((A)[2]); VEC3NEG((A)[0]); }

#define MAT3INC(A,B) { VEC3INC((A)[1],(B)[1]); VEC3INC((A)[2],(B)[2]); \
			 VEC3INC((A)[0],(B)[0]);}

#define MAT3DEC(A,B) { VEC3DEC((A)[1],(B)[1]); VEC3DEC((A)[2],(B)[2]); \
			 VEC3DEC((A)[0],(B)[0]);}

#define MAT3PLUS(A,B,C){ VEC3PLUS((A)[1],(B)[1],(C)[1]);VEC3PLUS((A)[2],(B)[2],(C)[2]);\
			    VEC3PLUS((A)[0],(B)[0],(C)[0]); }

#define MAT3MIN(A,B,C) { VEC3MIN((A)[1],(B)[1],(C)[1]);VEC3MIN((A)[2],(B)[2],(C)[2]); \
			    VEC3MIN((A)[0],(B)[0],(C)[0]); }

#define TRANSP3(A) { float h01,h02,h12;\
		        h01=(A)[0][1];     h02=(A)[0][2];     h12=(A)[1][2]; \
		    (A)[0][1]=(A)[1][0]; (A)[0][2]=(A)[2][0]; (A)[1][2]=(A)[2][1]; \
	            (A)[1][0]=h01;     (A)[2][0]=h02;     (A)[2][1]=h12; }
#define TRPASS(A,B)  { (A)[1][1]=(B)[1][1]; (A)[1][2]=(B)[2][1]; (A)[1][3]=(B)[3][1]; \
		       (A)[2][1]=(B)[1][2]; (A)[2][2]=(B)[2][2]; (A)[2][3]=(B)[3][3]; \
	               (A)[3][1]=(B)[1][3]; (A)[3][2]=(B)[2][3]; (A)[3][3]=(B)[3][3]; }

/* matrix multiplication */
/* Beware: don't do MATxVEC(a,b,a), since there is no local copy */

#define MAT3VEC(v,M,u) {\
     (v)[0]=(M)[0][0]*(u)[0] + (M)[0][1]*(u)[1] + (M)[0][2]*(u)[2];\
     (v)[1]=(M)[1][0]*(u)[0] + (M)[1][1]*(u)[1] + (M)[1][2]*(u)[2];\
     (v)[2]=(M)[2][0]*(u)[0] + (M)[2][1]*(u)[1] + (M)[2][2]*(u)[2];\
     }					/* v = Mu, column vectors */

#define VEC3MAT(v,u,M) {\
     (v)[0]=(M)[0][0]*(u)[0] + (M)[1][0]*(u)[1] + (M)[2][0]*(u)[2];\
     (v)[1]=(M)[0][1]*(u)[0] + (M)[1][1]*(u)[1] + (M)[2][1]*(u)[2];\
     (v)[2]=(M)[0][2]*(u)[0] + (M)[1][2]*(u)[1] + (M)[2][2]*(u)[2];\
     }					/* v = uM, row vectors */

/* Beware: don't do MATxMAT(a,b,a), since there is no local copy */


#define MAT3MAT(A,B,C) { int I,J,K;\
			   for(I=0; I<3; I++)\
			     for(J=0; J<3; J++)\
			       { A[I][J]=0.00;\
				 for (K=0; K<3; K++)\
				   A[I][J]+=B[I][K]*C[K][J];\
			       }\
		       }

/* following is for rotation of T about an arbitrary axis N (should be */
/* unit length vector. Result is a row-vector matrix (y=xR) */
#define ARBANGROTMAT(R,N,T)\
{ double C=cos((T)), S=sin((T));\
       (R)[0][0]=SQR((N)[0])+(1.0-SQR((N)[0]))*C;\
       (R)[1][1]=SQR((N)[1])+(1.0-SQR((N)[1]))*C;\
       (R)[2][2]=SQR((N)[2])+(1.0-SQR((N)[2]))*C;\
       (R)[0][1]=(N)[0]*(N)[1]*(1.0 - C)+(N)[2]*S;\
       (R)[0][2]=(N)[0]*(N)[2]*(1.0 - C)-(N)[1]*S;\
       (R)[1][2]=(N)[1]*(N)[2]*(1.0 - C)+(N)[0]*S;\
       (R)[1][0]=(N)[1]*(N)[0]*(1.0 - C)-(N)[2]*S;\
       (R)[2][0]=(N)[0]*(N)[2]*(1.0 - C)+(N)[1]*S;\
       (R)[2][1]=(N)[1]*(N)[2]*(1.0 - C)-(N)[0]*S;\
     }
