/*      getargs - command line argument processor for C programs
 * 
 * See getargs.doc and header file for full explanation
 * (C)  Copyright 1985, Allen I. Holub.  All rights reserved.
 * This program may be copied for personal, non-profit use only
 * 
 * history...
 * 
 * May 1985     published in Dr. Dobb's Journal #103.
 * 19 May 85    Transcribed by James R. Van Zandt
 * Oct 1993     Extended/modified by Philip lijnzaad@embl-heidelberg.de
 *
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

/*	getargs.h - typedefs and defines needed for getargs	*/

#ifndef _GETARGS_H_			/* avoid multiple inclusion */
#define _GETARGS_H_ 1			/* (unlikely, but anyway) */

#include "compat.h"
#include <stdio.h> 			/* for prototypes */
#include <stdlib.h>  
#include <string.h>

typedef struct arg_ {
  char * name;				/* name of command line switch 	  */
  int type;				/* variable type (see below)	  */
  void *address;			/* pointer to storage of variable */
  char  *comment;			/* pointer to error message 	  */
  int (*check)(void*);			/* function to be called to check */
					/* values (not called if NULL)    */
} arg_t, *arg_p;

/* prototype: */
int getargs(int argc, char **argv, arg_p tabp);

/* following may be set by the user; it's printed by pr_usage. If it's not */
/* set, a default is provided. If it's not set, the full table will be */
/* printed both upon error and with type HELP. If it's set, the full table */
/* wil only be printed with type HELP; upon error, only usage_message will */
/* be printed */

extern char * usage_message;

/* the option table should be something like 
 *
 *   arg_t Argtab[]={  
 *  /# name   type  address  comment (starts at 21st column) check_function #/
 *     {"-q", BOOL, &opt_q,  "do it quick'ndirty", NULL },
 *     {"-fast", BOOL, &opt_q,  "same as -q", NULL },
 *     {"-n", INT, &ntimes,  "do it N times" , check_that_nr},
 *     {"-r", DOUBLE, &radius,  "use radius <float>" , check_the_radius},
 *     {"-o", STRING, &filename, "output goes to STRING", &check_itslength},
 *     {"-p", PROC, &foobar, "invoke foobar (no arguments)", NULL },
 *     {"-format ", PROC1, &do_format, "format with STRING", NULL },
 *     END_ARGTAB
 *   };
 *
 *   argc = getargs(argc, argv, Argtab); /# invocation #/
 */

/* option types: */
#define HELP    1			/* only prints table and exits */
#define BOOL	2			/* loads 1 (or 0) in int*      */
#define INT	3			/* loads int in int*. Reads 0x<val> */
					/* as hexadecimal, and 0<val> as */
					/* octal */
#define DOUBLE	4			/* loads double in double*	*/
#define CHAR	5			/* loads char in char*		*/
#define STRING	6			/* loads char* in char** 	*/
#define PROC	7			/* calls void (*func)(void)     */
#define PROC1	8			/* calls (*func)(char *nextarg)	*/
#define R_FILE  9			/* loads FILE* in **FILE, opened */
					/* for reading; stdin if "-"  */
#define W_FILE  10			/* loads FILE* in **FILE, opened */
					/* for writing; stdout if "-"  */
#define A_FILE  11			/* for appending; stdout if "-" */
#define COMMENT 12			/* for introducing a comment in */
  /* the table of options. The name of the option itself is printed as a */
  /* comment; the rest of the slots are not used, so can be 0 */

#define END_ARGTAB { NULL, 0, NULL, NULL, NULL } /* use to terminate table */


/* CHECK FUNCTION-DEFINING MACROS:  */


#define DEF_NUM_LIM(FUNCNAME, TYPE, LOWERTEST, UPPERTEST, MESSAGE)\
  int FUNCNAME(void *param) { if( *(TYPE*)param LOWERTEST && \
				 *(TYPE*)param UPPERTEST )\
				 return(0);\
			       fprintf(stderr, "%s", MESSAGE);\
			       return(1); }
/* defines a numerical limits check function: */
/*     	(eg. DEF.(check_radius, double, >0.0, <=10.0, "invalid radius") */

#define DEF_STR_LIM(FUNCNAME, LOWERTEST, UPPERTEST, MESSAGE)\
  int FUNCNAME(void * charpp) { int i = strlen(*(char**)charpp);\
				if( i LOWERTEST && i UPPERTEST )return(0);\
				fprintf(stderr, "%s", MESSAGE);return(1); }
/* defines a string length check function: */
/*     	(eg. DEF.(check_filename, 1, FILENAME_MAX, "filename too long") */


#define DEF_CHCK_FUN(FUNCNAME, PARAM, TEST, MESSAGE)\
  int FUNCNAME(void * PARAM) { if((TEST))return(0); \
				fprintf(stderr, "%s", MESSAGE);return(1); }
/* general simple check functions. Be sure to use
 * proper casts inside TEST-expression, as PARAM has type void*
 * DEF_CHCK_FUN(check_opt_abc,		/# name of function to be defined #/
 * 	     dummy,			/# argument is not used #/
 * 	     opt_a + opt_b + opt_c <= 1, /# boolean expression. If 0, error #/
 * 	     "choose only one of -a -b -c" /# message #/
 * 	     )
 * checks if boolean options -a, -b and -c are used mutually exclusively
 */
#endif /* _GETARGS_H_ */
