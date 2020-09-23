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


#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include "getargs.h"

typedef int (*PFV)(void);
typedef int (*PFS)(char*);

#define USAGE_MESSAGE "Usage:"
char * usage_message;			/* printed by pr_usage; user-setable */
static int help=0;

#define OPT_FIRSTCHAR "-+"		/* first chars that imply an option */

#define ERR_EXIT_STATUS 1		/* exit status upon user error */
#define INT_ERR_EXIT_STATUS 2		/* exit status upon internal error */

static int boolval;			/* used in setarg */

static int setarg(arg_p argp, char ** argv) {
  /* set an argument.  argp points at the argument table entry corresponding
   * to argv[0]. Returns number of items that now can be taken out of
   * argument list, ie. 1 for BOOLEAN and PROC, 2 for rest, and 0 for
   * error */
  char *cp;
  PFV pfv;				/* pointer to function of void  */
  PFS pfs;				/* pointer to function of string */
  double g;
#define USER_CHECK if(argp->check&&(*argp->check)(argp->address )){\
  fprintf(stderr, " (option %s)", *argv);return(0); }


  if (argp->type==HELP) {
    help=1;
    return(0);
  }

  /* check address: */
  if (!argp->address) {
    fprintf(stderr, "getargs internal error: bad address for option %s\n",
	    argp->name); 
    exit(INT_ERR_EXIT_STATUS);
  }

  /* following two option types don't take values:  */
  switch (argp->type) {
  case BOOL:				/* set it */
    *(int*)argp->address = boolval;	/* 1 if -option, 0 if +option */
    USER_CHECK;
    return(1);
    
  case PROC:
    pfv = argp->address;		/* call it, without arg */
    if ( (*pfv)() ) {			/* error if non-zero return */
      fprintf(stderr, " (option %s)", *argv);
      return(0); 
    }
    USER_CHECK;				/* also call normal check function */
    return(1);

  default:
    break;
  }

  /* following functions all need a value, so see if it's at all there: */
  if ( (! argv[1]) || (!argv[1][0] ) ) {
    fprintf(stderr, "option %s needs additional argument", argv[0]);
    return(0);
  }

  switch (argp->type) {
  case INT:
    *(int*)argp->address = strtol(argv[1], &cp, 0);
    if (*cp) {
      fprintf(stderr, "funny integer: %s\n", argv[1]); 
      return( 0 ); 
    }
    USER_CHECK;
    return(2);
  
  case CHAR:
    if ( strlen(argv[1]) > 1 ) {
      fprintf(stderr, "funny character: %s\n", argv[1]); 
      return( 0 );
    }
    *(char*)argp->address = argv[1][0];
    USER_CHECK;
    return(2);

  case STRING:				/* just read a string */
    *(char **) argp->address = argv[1];
    USER_CHECK;
    return(2);

  case DOUBLE:
/*     if (sscanf(argv[1], "%lf%.1s", (double*)argp->address, &dummy) != 1 ) {
 *       fprintf(stderr, "funny double: %s\n", argv[1]); 
 *       return(0); 
 *     }
 */
    g=strtod(argv[1], &cp);		/* doesn't work under SunOS !?! */
    *(double*)argp->address = g;
    if (*cp) {
      fprintf(stderr, "funny double: %s\n", argv[1]); 
      return(0); 
    }
    USER_CHECK;
    return(2);

  case PROC1:
    pfs = argp->address;		/* call it with argument argv[1] */
    if( (*pfs)(argv[1]) ) {		/* error if non-zero return */
      fprintf(stderr, " (option %s)", *argv);
      return(0); 
    }
    USER_CHECK;
    return(2);

  case R_FILE:				/* open arg as file for reading */
    if ( ! strcmp(argv[1], "-") )	/* special case: stdin */
      *(FILE**)argp->address = stdin;
    else
      *(FILE**)argp->address = fopen(argv[1], "r");

    if (! *(FILE**)argp->address) {
      fprintf(stderr, "%s ", argv[0]);	/* tell which option errors */
      perror(argv[1]);			/* and use system to tell what error */
      return(0);
    }
    USER_CHECK;
    return(2);

  case W_FILE:				/* open arg as file for writing */
    if ( ! strcmp(argv[1], "-") )	/* special case: stdout */
      *(FILE**)argp->address = stdout;
    else
      *(FILE**)argp->address = fopen(argv[1], "w");

    if (! *(FILE**)argp->address) {
      fprintf(stderr, "%s ", argv[0]); 	/* tell which option errors */
      perror(argv[1]);
      return(0);
    }
    USER_CHECK;
    return(2);

  case A_FILE:				/* open arg as file for appending */
    if ( ! strcmp(argv[1], "-") )	/* special case: stdout */
      *(FILE**)argp->address = stdout;
    else
      *(FILE**)argp->address = fopen(argv[1], "a");

    if (! *(FILE**)argp->address) {
      fprintf(stderr, "%s ", argv[0]); 	/* tell which option errors */
      perror(argv[1]);
      return(0);
    }
    USER_CHECK;
    return(2);

  default:
    break;
  }

  fprintf(stderr, "getargs internal error: bad argument type %d\n",
	  argp->type);
  exit(INT_ERR_EXIT_STATUS);
} /* setarg */


static arg_p findarg(char * name, arg_p tabp) {
  for ( ; tabp->name ; tabp++) {
    if (tabp->type != COMMENT && ! strcmp(tabp->name, name) )
      return (tabp);
  } 
  return (NULL);
} /* findarg */

static void pr_usage(arg_p tabp) {
  /* if this is a HELP request or usage_message is not set, print out the */
  /* table in the form  <optionname> [<arg>] <comment> [value]. Otherwise, */
  /* just print the set usage_message */ 

  if (!help && usage_message) {		/* this is an error, not a HELP type */
    fprintf(stderr, "\n%s\n", usage_message);
    return;
  }

  if(!usage_message)			/* provide some default */
    usage_message = USAGE_MESSAGE;
  fprintf(stderr, "\n%s\n", usage_message);

  fprintf(stderr, "Options:\n");

  for ( ; tabp->name; tabp++) {
    switch (tabp->type) {
    case COMMENT:
      fprintf(stderr, "%s\n", tabp->name);
      break;

    case HELP:				/* has no default */
      fprintf(stderr, "  %8s          %s\n", tabp->name, tabp->comment);
      break;

    case BOOL:
      fprintf(stderr, "  %8s          %s [%s]\n", tabp->name, tabp->comment,
	      (*(int*)(tabp->address) !=0 ) == (tabp->name[0]=='-') ? 
	      "true" : "false");	/* always a default */
      break;

    case PROC:				/* no default */
      fprintf(stderr, "  %8s          %s\n", tabp->name, tabp->comment);
      break;

    case INT:
      fprintf(stderr, "  %8s <n>      %s [%d]\n", /* always a default */
	      tabp->name, tabp->comment, *(int*)(tabp->address));
      break;

    case CHAR:				/* check for default */
      fprintf(stderr, "  %8s <char>   %s ", tabp->name, tabp->comment);
      if ( *(char*)tabp->address)	/* there's a default */
	fprintf(stderr, "['%c']" , *(char*)(tabp->address));
      fprintf(stderr, "\n");
      break;

    case STRING:			/* check for default (and for NULL) */
      fprintf(stderr, "  %8s <string> %s ", tabp->name, tabp->comment);
      if (*(char**)(tabp->address)	/* treat NULL also as empty string */
	  && **(char**)(tabp->address)) /* not an empty string */
	fprintf(stderr, "[\"%s\"]", *(char**)(tabp->address));
      fprintf(stderr, "\n");
      break;

    case DOUBLE:
      fprintf(stderr, "  %8s <float>  %s [%.3g]\n", /* always a default */
	      tabp->name, tabp->comment, *(double *)(tabp->address));
      break;

    case PROC1:
      fprintf(stderr, "  %8s <string> %s\n", tabp->name, tabp->comment);
      break;

    case R_FILE:
    case W_FILE:
      fprintf(stderr, "  %8s <file>   %s\n", tabp->name, tabp->comment);
      break;
    }
  }
} /* pr_usage */

int getargs(int argc, char **argv, arg_p tabp) {

/* Process command line arguments, stripping all command line switches out */
/* of argv.  Return a new argc. If an error is found, exit(ERR_EXIT_STATUS) */ 
/* is called (getargs won't return) and a usage message is printed, */
/* showing all arguments in the table */


  int nargc, i;
  char **nargv, *cp, oneletter[3], *savarg;
  arg_p argp;

  nargc = 1;				/* argv[0] is program name  */
  memset(&oneletter[0], 0, sizeof(char[3])); /* make zero */
/* OR #include <bstring.h> :  bzero(oneletter, 3); */

/* #define USAGE_EXIT do { pr_usage(tabp); exit(ERR_EXIT_STATUS); } while(0)
 */

#define USAGE_EXIT pr_usage(tabp),exit(ERR_EXIT_STATUS)

  for (nargv = ++argv; --argc > 0; argv++) { /* walk all arguments */
    /* for now no provision for '\-'; do it with '--' */
    if ( strchr( OPT_FIRSTCHAR, (int)*argv[0]) ) { /* an option (or more) */
      if (! strcmp("--", *argv)) {	/* signals end of options */
	argv++;				/* skip this one */
	nargc += (argc-1);		/* all the rest but '--' */
	for ( ; --argc > 0 ; argv++)	/* shift the rest */
	  *nargv++ = *argv;
	return(nargc); 
      }
      if (! argv[0][1] ) {		/* "-" or "+": never an option, */
	*nargv++ = *argv;		/* always an argument: shift it  */
	nargc++;
	continue;
      }
      boolval = (*argv[0] == '-');	/* 1 if -option, 0 if +option */
      if ( (argp = findarg(*argv, tabp)) ) { /* first try find this arg as */
					  /* a complete word in the table */ 
	if ( (i=setarg(argp, argv)) ) {
	  i--;
	  argv += i;
	  argc -= i;
	  continue;
	} else
	  USAGE_EXIT;			/* wrong or no argument */
      } else {				/* try finding it as a letter arg */
	oneletter[0]= argv[0][0];	/* any of OPT_FIRSTCHAR (typically -)*/
	for ( cp = &argv[0][1]; *cp; cp++ ) {
	  oneletter[1] = *cp;
	  if ( (argp=findarg(oneletter, tabp)) ) { /* find this arg in table */
	    /* we seem to have found a matching one-letter option. */
	    /* Temporarily change this *argv, so as not to confuse */
	    /* the error messages: */
	    savarg = *argv;		/* save current arg */
	    *argv = &oneletter[0];	/* argv[1] still 'normal' */
	    if ( (i = setarg(argp, argv)) ) { /* not-nil: succes !  */
	      *argv = savarg;		/* restore original */
	      if (i==2) {		/* option that takes a value: */
		argv++;			/* must have been eaten */
		argc--;
		break;			/* continues at upper for-loop */
	      }	/* else: continues at cp++ */
	    } else			/* 0: failure */
	      USAGE_EXIT;		/* wrong or no argument */
	  } else {
	    fprintf(stderr, "Unrecognized option \'%s\'\n", *argv);
	    USAGE_EXIT;			/* option not found */
	  }
	}
      }
    } else {				/* not an option: shift arguments */
      *nargv++ = *argv;
      nargc++;
    }
  } 
  return(nargc);
} /* getargs */
