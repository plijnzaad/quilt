/* ------------------------------------------------------------------ *
 * Permission to copy and/or modify all or part of this work is, for  *
 *   non-profit purposes, granted, provided that this NO WARRANTY     *
 *   and this copyright notice are retained verbatim and are          *
 *   displayed conspicuously. If anyone needs other permissions that  *
 *   are not covered by the above, please contact the author          *
 *      							      *
 * NO WARRANTY: THIS WORK IS PROVIDED ON AN "AS IS" BASIS.  THE       *
 *   AUTHOR PROVIDES NO WARRANTY WHATSOEVER, EITHER EXPRESS OR        *
 *   IMPLIED, REGARDING THE WORK, INCLUDING WARRANTIES WITH THE WORK, *
 *   INCLUDING WARRANTIES WITH RESPECT TO ITS MERCHANTABILITY OR      *
 *   FITNESS FOR ANY PARTICULAR PURPOSE.			      *
 *								      *
 * Philip Lijnzaad, lijnzaad@embl-heidelberg.de			      *
 * ------------------------------------------------------------------ */

/* module that is meant to contain all compatibility stuff */
#ifndef _COMPAT_H_ 
#define _COMPAT_H_ 


#ifdef SunOS
  /* exspect LOADS of warnings about unknown functions ! */
#  define  ulong unsigned long
#  define  FILENAME_MAX 255		/* from OSF1-stdio.h */
  /* as yet, no work around for uniform use of stdarg.h <-> varargs.h. So */
  /* forget SunOS for now ... */
#endif /* SunOS */


#ifdef Darwin                           /* aka MacOSX */
#  define  ulong unsigned long
#  define  FILENAME_MAX 255		/* from OSF1-stdio.h */
#endif /* Darwin */

#ifdef __APPLE__
#  define  ulong unsigned long
#  define  FILENAME_MAX 255		/* from OSF1-stdio.h */
#endif /* Darwin */



#endif /* _COMPAT_H_ */

