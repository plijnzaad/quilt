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

/* Some general purpose utilities */

#ifndef _UTILS_H_			/* avoid multiple inclusion */
#define _UTILS_H_ 1

/* #include "compat.h" */
/* types: */
#include <stdio.h>			/* needed for some typedefs */

/* debugging version of fopen/fclose: writes stuff to stderr:
 */
FILE* _myfopen(const char *filename, const char *mode);
#define FOPEN(name,mode) (UPDATE_FILE_AND_LINE,_myfopen(name,mode))

void _myfclose(FILE * file);
#define FCLOSE UPDATE_FILE_AND_LINE,_myfclose


/* #if defined sgi || defined SUN4 /# SGI, sun4 may have defined uint: #/
 * #include <sys/types.h>
 * #else
 *   typedef unsigned int uint;
 * #endif
 */

#include <sys/types.h>			/* for uint defintion */
#ifdef  __ultrix__
typedef unsigned long ulong;		/* they forgot that */
#endif
#ifndef __alpha__
typedef unsigned char uchar;
#endif

#ifdef TRUST_ALLOCA			/* if alloca (stack allocation) */
#  define ALLOCA alloca
#  define AFREEA(X)
#else 
#include <alloca.h>
#  define ALLOCA malloc			/* not available or trusted */
#  define AFREEA free
#endif 

/* #define AFREEA                          /# @@@debugging */


#ifndef while_not
#  define while_not(x) while(!(x))
#endif

#ifndef if_not
#  define if_not(x) if(!(x))
#endif

#define BOOLIFY(X)  ((X) != 0)		/* make boolean out of expr */
#define XOR(A,B)    (BOOLIFY(A) != BOOLIFY(B)) /* exclusive or */
#define LOGEQV(A,B) (BOOLIFY(A) == BOOLIFY(B)) /* log. equivalence  */

#define LOOPUP(I, N) for((I)=0; (I)<(N); (I)++)
#define LOOPDOWN(I, N) for((I)=(N)-1; (I)>=0;  (I)--)

/* general warning and debugging stuff: */

extern const char * __file__;		/* declared versions of macros */
extern int  __line__;			/* __FILE__  and __LINE__ */
extern char message_message[];		/* written to by error() */
#define UPDATE_FILE_AND_LINE ((__file__=__FILE__),(__line__=__LINE__))

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define R_ok 0				/* general return value if all's ok */
#define R_err 1				/*  ..     ..     .. if not all's ok */


#if defined __STDC__ || defined __ANSI__ 
#  define QUOTE(X) #X
#else
#  define QUOTE(X) "X"
#endif 

#if defined (ultrix) && !defined (__GNUC__)
#  define CONCAT(x,y) x/**/y
#elif defined (VMS) && !defined (__ALPHA)
#  define CONCAT(x,y) x/**/y
#elif defined (CCKR)			/* K&R C  */
#  define CONCAT(x,y) x/**/y
#else  				
#  define CONCAT(x,y) x##y		/* ANSI C */
#endif

#ifdef __STDC__
#include <assert.h>
#else
extern int __assert(const char * assertion, const char* file, const int line);
  /* in case abort() doesn't work  */
#  ifndef NDEBUG
#    define assert(X) if (X); else\
	__assert(QUOTE(X), __FILE__, __LINE__)
#  else
#    define assert(x) ((void)0) 
#  endif
#endif

#define assert_not(X) assert(! (X))


#define BUG die("Bug at %s line %d\n", __FILE__, __LINE__)
/* following functions are much like the perl ones.  */
#define WARN UPDATE_FILE_AND_LINE,warn
#define DIE UPDATE_FILE_AND_LINE,die
#define PERROR UPDATE_FILE_AND_LINE,myperror

#define ABORT UPDATE_FILE_AND_LINE,die	/* old style */
#define WARNING UPDATE_FILE_AND_LINE,warn /* old style */

					/* i don't see a way around this ...*/
#ifdef SUN4				/* SUN does not know raise() ! */
extern int raise(int __sig);		/* send __sig to your self */
#endif

/* following is a bit like perl functions warn and die; they will put file */
/* and linenr after printf-style formatted message, if it lacks a trailing */
/* newline or space character */

#define WARN_TERMINATOR ":, "		/* if message ends with this, FILE & */
					/* LINE will be appended  */
int set_err_file(const char *filename);
  /* set file for error output; if NULL or "", reset to stderr */

FILE* get_err_file(void);
  /* find out what the err file is; will return stderr + warning if it's
     not was set explicitly with set_err_file() */

extern int die(const char * fmt, ...);	
  /* writes to err_file and dies by raising DIE_SIG. If last character of */
  /* formatted string does equal any of WARN_TERMINATOR, does */
  /* fprintf(stderr, " at line %s line %d\n", __file__, __line__) too */ 

extern int myperror(const char * fmt, ...);
/* like WARN, but does a perror("") first */

extern int warn(const char * fmt, ...);
  /*  same as die, but doesn't ... */

extern int message(int status, const char *fmt, ...); 
  /* writes stuff to extern char message_message[], and returns status */



/* STRING HANDLING: */

#define MAXLINLEN 512

char * findfile(const char * Path, const char * File, const char * Mode,
		FILE ** filep); 
/* searches File in Path to be opened in Mode and returns (staticly stored)
 * full pathname or NULL. Path can be a colon separated list of
 * directories, or the name of an environment variable with such a list as
 * value. If Path is NULL or "", the path defaults to '.:$HOME'. If
 * filep!=NULL, *filep is set to the opened file
 */

void * _memsav(const void * data, int length, const char *, int); 
#define memsav(D, L) _memsav( (D), (L), __FILE__, __LINE__)
  /* allocates memory and */
  /* copies data of length; if len==0, treats it as a 0-terminated string */
#define MEMSAV(TO, FROM, L) \
	do { int __L = ((L)?(L):strlen((FROM))) + 1; \
	     TO=(char*)ALLOCA(__L); \
	     memcpy(TO,FROM,__L); } while(0)
  /* like memsav but with alloca */

char * _strsav(const char * s, const char *, int);
#define strsav(S) _strsav((S), __FILE__, __LINE__)
  /* strcpy(malloc(strlen(s)+1), s); */
  /* make dynamically allocated copy of string */

void * _readfile(FILE * file, const char *, int);
#define readfile(F) _readfile((F), __FILE__, __LINE__)
  /* reads complete contents of FILE into dynamically allocated memory,
     which is returned */ 
void * _readfilename(const char *path, char * filename, const char *, int);
#define readfilename(P,N) _readfilename((P), (N), __FILE__, __LINE__)
  /* reads complete contents of FILENAME (to be found in PATH)into
     dynamically allocated memory, which is returned */

extern char * charbits(uchar c); 
  /* returns c as a string containing binary numbers */
extern char * intbits(uint i);  
  /* returns i as a string containing binary numbers */

extern char * strrev(const char * s);	/* malloc'ed reversed string  */
extern const char * strrstr(const char * search, const char * string);
  /* reverse strstr: strrstr("bc", "abcbcd") returns address of 2nd b */
extern int strncpy0(char * s1, const char * s2, int l);
  /* like strncpy, but always stores a '\0' at end, and returns length of */
  /* resulting string*/
extern int lstrcpy(char *  s1, const char * s2);
  /*like strcpy, but returns result length*/
extern int lstrcat(char *  s1, const char * s2);
  /*like strcat, but returns result length*/
extern int rplcpy(char * dest, const char * src, char c, char t);
  /* copies string src to dest, while changing c to t; if t=='\0' c is deleted
   * returns number of occurences of c */
extern int leastrp(char * string);
  /* strips leading whitespace off string; returns length of result */
extern int trastrp( char * string);
  /* strips trailing whitespace off string; returns length of result */
extern int allstrp( char * string);
  /* strips leading and trailing whitespace; returns length of result */
extern int chardel(char *string, const char* delete);
  /* deletes all characters in the set delete from string; returns new length */

extern int readpast(FILE * file, char * line, const char * match);
  /* does fgetl(line, file), and reads until it found match; returns length */

int fill(char *text, const char * pre, const char * breakon, int width);
  /* fills TEXT to WIDTH, line-breaking after characters in BREAKON (on
   *  " \t,.;:-" if BREAK == "") prefixing the lines with PRE. TEXT is
   *  modified (will grow, generally) and the nr of lines is
   *  returned. 
   */

extern int _i_arg1; 
extern float _f_arg1, _f_arg2;			
extern double _d_arg1, _d_arg2;			

#ifndef SQR
#define SQR(A) ((A)*(A))
#define _SQR(A)   (_arg1=(A), SQR(_arg1)) /* to avoid double function evals */
#define CUBE(A) ((A)*(A)*(A))
#define _CUBE(A) (_arg1=(A_), CUBE(_arg1))
#endif

#define ABS_ORD_MAGN(A) (fabs(floor(log10(fabs(A))))) /* when bypassing %g */

#define ABS(I) abs((I))
#define FABS(I) fabs((I))		/* gcc's fabs faster than macro! */
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define _MINf(a,b) (_f_arg1=(a),_f_arg2=(b),(_f_arg1<_f_arg2?_f_arg1:_f_arg2))
#define _MAXf(a,b) (_f_arg1=(a),_f_arg2=(b),(_f_arg1>_f_arg2?_f_arg1:_f_arg2))
#endif

#define iSWAP(A,B) do { int i=(A);(A)=(B);(B)=i; } while (0)
#define fSWAP(A,B) do { float i=(A);(A)=(B);(B)=i; } while (0)
#define dSWAP(A,B) do { double i=(A);(A)=(B);(B)=i; } while (0)

/**********************************/

#ifdef __ultrix__
/* prototypes that were forgotten by ultrix: */
void     bcopy(const char *, char *, int);
int      bcmp(const char *, const char *, int);
void     bzero(char *, int);
void     blkclr(char *, int);
#endif /* __ultrix__  */

#ifndef _BSD
#  define _BSD 1
#endif

#ifndef RAND_MAX
#  define RAND_MAX 2147483647	/* OSF/1 kludge: check this !!#!@ */
#endif

int random_range(int n, int *array);
/* randomly fill array with the range 0 .. N-1 */
int randomize(int n, int size, void *array);
/* randomizes, in place, an ARRAY of N things having SIZE */

void qrank(const void *array, int nelts, int size, 
	   int (*compare)(const void *, const void *),
	   int * indices, int * ranks);   
  /*
   * does    : Sorts indices and determines the ranks of elements of an array
   *
   * gets    : the ARRAY to be looked at, which has NELTS elements of size
   *   SIZE, to be compared using the function COMPARE (exactly as in 
   *   qsort). INDICES or RANKS are not modified if they are NULL.
   *
   * affects : The INDICES, if not NULL, are sorted such that 
   *   ARRAY [ INDICES[ 0 <= i < NELTS ] ] is ordered. The RANKS[i], if
   *   not NULL, are set to the rank of element ARRAY[i] 
   *
   * returns : void
   *
   * warns   : if both INDICES and RANKS are NULL
   *
   * comment : See qsort() for details. To sort ARRAY (to a copy),
   *   consider one of the following:
   *
   * for (i=0; i<nelts; i++) { 
   *   sorted_array_1 [ i ] = array[ indices[i] ];
   *   sorted_array_2 [ ranks[i] ] = array[i];
   * }
   */

/* compare functions for sorting: */
int double_cmp(const void * a, const void *b); 
int float_cmp(const void * a, const void *b);
int int_cmp(const void * a, const void *b);
void debug_pause(void);			/* writes pid to stdout, and waits */
					/* 2 secs to be grabbed */

/* some macro tricks: */
#define BIT(N) (1<<(N))			/* BIT(0)=1, BIT(1)=2, BIT(3)=4,etc. */
/* to store/retrieve bit-item I of size Q (bits) in/from an uint array A :
 */

#define SETBITS(A,I,Q,D) A[(I/(32/Q))] |= \
  ( (D&(BIT(Q)-1))<<((I-(I/(32/Q))*(32/Q))*Q) )
#define GETBITS(A,I,Q) ((A[(I/(32/Q))])>>((I-(I/(32/Q))*(32/Q))*Q)&(BIT(Q)-1))


#endif /* _UTILS_H_ */

