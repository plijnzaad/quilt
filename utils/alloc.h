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

/* alloc/dynamic memory watch: */

/* allocating functions,  give decent error message and quit if out of */
/* memory. Also, they write stuff to the file DMW_LOG ,  if (1) dmw_on is */
/* non-zero and (2) DOdmw(1) has been called. */ 

#ifndef _ALLOC_H_			/* no collision with system files ?? */

#define _ALLOC_H_			/* avoid multiple inclusion */

#include "compat.h"

#define DMW_DIR  "/junk/lijnzaad"	/* where to leave logfile */
#define DMW_FILE "dmw.log"		/* log file name. Gets big */
#define DMW_LOG  DMW_DIR"/"DMW_FILE	/* full pathname */

#ifndef GC				/* if no garbage collection */

#define DMW_MODE "w"			/* opening mode for log file; can */
					/* also be "a", to append */

extern int dmw_on;			/* controls overal dynamic memory */
					/* watching. Set it from debuyger */

/* allocating safely: following functions give decent error message and */
/* quit if out of memory. Also, they write stuff to the file */
/* DMW_LOG if (1) dmw_on is non-zero and (2) DOdmw(1) has been called. */

void * _malloc(int bytes, const char * filename, int linenr); 
#define MALLOC(n) _malloc((n), __FILE__,__LINE__)

void * _calloc(int nelem, int elsize, const char *, int); 
#define CALLOC(n, size) _calloc((n), (size), __FILE__,__LINE__)

void * _realloc( void * ptr, int nbytes, const char *, int); 
#define REALLOC(ptr, size) _realloc((ptr), (size), __FILE__,__LINE__)

void * arealloc(void*ptr, int nbytes, const char *, int); 
#define AREALLOC(ptr, size) arealloc((ptr), (size),__FILE__,__LINE__)
    /* does calloc(1, size) if ptr==NULL,  REALLOC() otherwise */

void _free(void * ptr, const char *, int);
#define FREE(ptr) _free((ptr), __FILE__,__LINE__),(ptr)=NULL

void _dodmw(int onoff, const char * filename, int lin);
#define DOdmw(onoff) _dodmw((onoff), __FILE__,__LINE__)
/* Controls watching of dynamic memory in THIS piece of program text. */
/* Overall switching is done from debugger by setting dwm_on to non-zero */
/* analyzing the resulting file "dmw_trail" done by program dmw_ana */


#else  /* #defined GC: want Boehms garbage collection */

/* #include <gc.h> ??? */

#ifdef SUNOS
#define strtok dangerous !		/* messy strtok for shared lib's  */
#endif

#include "gc.h"

#define FREE(X) GC_free((X))		/* not necessary, but may speed up */
#define _free(F,L,P) GC_free((P))	/* forget about file/lines for now */

#define MALLOC GC_malloc
#define CALLOC(A,B) GC_malloc((A)*(B))
#define REALLOC(ptr, size) GC_realloc((ptr), (size))

#endif /* #ifdef GC */

#endif /* #ifndef _ALLOC_H_ */
