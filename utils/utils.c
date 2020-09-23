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

/* #include <stdarg.h> */
#include "compat.h"
#include <stdarg.h>
#ifndef _VA_LIST_
#define _VA_LIST_			/* for avoiding errors with gcc ?!!? */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined ultrix && ! defined  __GNUC__ /* annoying ultrix (infrix ...) */
#undef  __STDC__			/* forces ansi_compat.h to be read */
#include <signal.h>			/* avoid prototype of signal(); */
#include <limits.h>			/* error in <limits.h> !!!! */
#define __STDC__ 
#else
#include <signal.h>
#include <limits.h>
#endif

#include <time.h>
#include <ctype.h>

#include "utils.h" 
#include "alloc.h"

#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>

static struct { 
  /* variables/data local to this module */
  int (*qrank_cmp)(const void*, const void*); /* qsort-type compare function */
  const void * qrank_array;
  int qrank_size;
} local;

/* general arguments for macros in utils.h (bloody awful) */
extern int _i_arg1; 
extern float _f_arg1, _f_arg2;			
extern double _d_arg1, _d_arg2;			

/* local prototype */
static void _warn(const char * default_message, const char *fmt, va_list args);

char * charbits(uchar C) {
  char  *cp;
  static char  word[32+1];
  int i;

  word[32]=0; cp = &(word[32]);
  for (i=1; i<=8; i++) 
    *cp-- = (C & BIT(i))?'|':'.' ;
  return (cp+1);
}

char * intbits( uint I) {
  int i; 
  char *cp;
  static char  word[32+1];

  word[32]=0; cp = &(word[32]);
  for (i=0; i<=31; i++)
    *cp-- = (I & BIT(i) )?'|':'.' ;
  return (cp+1);
}

int strncpy0(char * s1, const char * s2, int l) {
  int i=0;

  while( (s1[i]=s2[i]) && (i <l))
    i++;
  s1[i]='\0';
  return(i);
}

int lstrcpy(char *  s1, const char * s2) {
  int i=0;
  while ( (s1[i]=s2[i]) )
    i++;
  s1[i]='\0';
  return(i);
}

int lstrcat(char * s1,  const char * s2) {
  int i=0;
  while(s1[++i])
    ;
  while( (s1[i]=s2[0]) ) {
    i++; s2++;
  }
  s1[i]='\0';
  return (i);
}

int rplcpy(char * dest, const char * src, char c, char t) {
  int i=0, j=0, k=0;
  while (src[i]) {
    if (src[i]==c) {
      if (t)
	dest[j++]=t;
      i++;
      k++;
    }
    else
      dest[j++]=src[i++];
  }
  dest[j]=0;
  return(k);
}

int leastrp(char *  string) {
  char *cp=string, *cp2=string;
  if (! string[0])			/* empty */
    return(0);
  while ( isspace(*cp2) )		/* walk to non-whitespace */
    cp2++;
  while( (*cp = *cp2) )			/* copy */
    cp++, cp2++;
  return(strlen(string));
}

int trastrp(char * string) {
  int i=strlen(string);
  if (!i)				/* empty string */
    return(0); 
  i--;

  while ( isspace(string[i] ))
    i--;
  string[i + 1]='\0';
  return(i);
}

int allstrp(char * string) {
  trastrp(string);
  return leastrp(string) ;
}

int chardel(char *string, const char * delete) {
  int l;
  char *s, *temp, *t, c;

  l=strlen(string);
  MEMSAV(temp, string, l);
  t=temp; s=string;

  l=strlen(delete);
  if (l==1) { 				/* common case: make faster */
    c=delete[0];
    while( *t ) {
     if ( *t != c)
       *s++ = *t;
     t++;
   }
  } else { 
    while( *t ) {
      if ( strchr(delete, *t) )
	*s++ = *t;
      t++;
    }
  }
  l= s - string;			/* new length */
  *s='\0';
  AFREEA(temp);
  return l;				/* length */
} /* chardel */


#ifdef SUN4				/* they don't have ise() */
#include <sys/types.h>			/* for fstat etc. */
int raise(int sig) {
  return( kill(getpid(), sig));
#define _raise raise
}
#endif

/* following are used as run-time versions of the corresponding uppercase */
/* macros, updating done through UPDATE_FILE_AND_LINE in the other macros */

#define ERR_MSG_LENGTH 1024		/* size of error buffer */

const char * __file__ = NULL;         /* run-time version of __FILE__ */
int __line__ = 0;                     /* id. of __LINE__ */

char message_message[ ERR_MSG_LENGTH ];	/* program global variable */

static FILE *err_file = (FILE*) NULL;   /* not to be seen outside */
/* report on this (18 jan 2014)
  gcc -g -Wall -I -DLinux                   -c -o utils.o utils.c
  utils.c:208:1: error: initializer element is not constant
  make: *** [utils.o] Error 1
*/

FILE* get_err_file(void) {
  if (err_file==NULL) {
    fprintf(stderr, "*** err_file not set explicitly using set_err_file(); "
            "using stderr instead\n");
    set_err_file(NULL);
  }
  return err_file;
}

int set_err_file(const char *filename) { 
  /* set file for error output; if NULL or "", reset to stderr */

  if ( filename== NULL || filename[0]=='\0' ) { 
    if (err_file) { 
      /* assume it was a proper FILE*; if not, we're hosed */
      fclose(err_file);
      err_file=NULL;
    }
    err_file=stderr;
    return 0;
  }

  err_file=fopen(filename, "w");
  if (!err_file)
    perror(filename),exit(1);		/* this one does go to stderr ... */
  return 0;
} /* set_err_file */

static void _warn(const char * default_message, 
		  const char *fmt, va_list args) {
  /* used by warn and die */
  int len;
  char warn_message[ ERR_MSG_LENGTH ];	/* local version */

  vsprintf(warn_message, fmt, args);		/* global buffer */

  len = strlen(warn_message);
  if (len > ERR_MSG_LENGTH ) {
    warn_message[ 73 ]='\0';		/* otherwise infinite recursion */
    DIE("internal error: error string too long:\n>%s ...", 
	ERR_MSG_LENGTH-1, warn_message);
  }

  fprintf(get_err_file(), "%s", len==0 ? default_message : warn_message);

  /* append the at <file> <line> bit */
  if (len==0 || strchr(WARN_TERMINATOR , warn_message[len-1]) )
    fprintf(get_err_file(), " at %s line %d\n", __file__, __line__);

  fflush(stdout);			/* so we have output already */
  fflush(stderr);
  fflush(get_err_file());

  message_message[0]='\0';		/* clear global variable */

  return;
} /* _warn */

FILE * _myfopen(const char * name, const char  *mode) {
  FILE *file=NULL;
  fprintf(stderr, "Opening file %s at %s line %d ",
          name, __file__, __line__);
  file=fopen(name, mode);
  fprintf(stderr, "yielded %p\n", file);
  return file;
}

void _myfclose(FILE * file) {
  fprintf(stderr, "Closing FILE* %p at %s line %d\n",
          file, __file__, __line__);
  fclose(file);
  return;
}


int myperror(const char *fmt, ...) {
  va_list args;

  va_start(args, fmt);
  _warn("i/o-error: ",                    /* default message if fmt+args empty */
	fmt, args);
  va_end(args);
  perror("");

  return 1;
} /* die */


int warn(const char *fmt, ...) {
  va_list args;

  va_start(args, fmt);
  _warn("Warning: something wrong",	/* default message if fmt+args empty */
	fmt, args);
  va_end(args);
  
  return 1;
} /* warn */

int die(const char *fmt, ...) {
  va_list args;

  va_start(args, fmt);
  _warn("Died",				/* default message if fmt+args empty */
	fmt, args);
  va_end(args);

  raise(SIGABRT);
  return 1;
} /* die */


int message(int status, const char *fmt, ...) {
  int len;
  va_list args;

  len=strlen(message_message);
  va_start(args, fmt);
  vsprintf(message_message+len, fmt, args); /* append */
  va_end(args);

  return status;
} /* message */

#ifdef OWN_ASSERT			/* avoid multiple definition */
int __assert(const char * assertion, const char* file, const int line) {
  die("%s: %d: assertion failed (programming bug):\n\t%s\n", 
      file, line, assertion);
  return(1);
} /* __assert */
#endif /* OWN_ASSERT */

char * _strsav(const char * s, 
	      const char * file, int line) {		/* save string s in dynamic memory */
  int l=strlen(s);
  char*cp=_malloc(l+1, file, line);
  memcpy(cp, s, l+1);
  return (char*) cp;
}

void * _memsav(const void * data, int len, const char * file, int line) {
  int l = len ? len : (strlen((char*)data)+1);
  char *cp = _malloc(l, file, line);
  return (void*) memcpy((void*)cp, data, l);
}

char * strrev(const char * str) {	/* return mirror image of string */
  int i; 
  int l=strlen(str);
  char * rev;
  
  rev=(char*)MALLOC(l+1);		/* dynamically allocated */
  for (i=0; i<l; i++)
    rev[l-1 -i]=str[i];
  rev[l]='\0';
  return(rev);
}

const char * strrstr(const char * se, const char * str) { /* reverse strstr */
  int i, lse, lstr;
  char * revse, *revstr;
  const char * cp; 

  lse=strlen(se);
  lstr=strlen(str);

#define STRREV(str, len, rev)			\
  rev= ALLOCA(len+1);				\
  for (i=0; i<len; i++)				\
    rev[ len-1-i ]=str[ i ];			\
  rev[ len ]='\0';

  STRREV(str, lstr, revstr);
  STRREV(se, lse, revse);

  cp=strstr(revse, revstr);

  if(cp)				/* NULL can just be passed */
    cp = str + lstr - ( lse + (uint)(cp - revstr) ) ;   

  AFREEA(revse);
  AFREEA(revstr);

  return(cp);
} /* strrstr */

#undef STRREV

#define _ACCESS(TYPE, VAR) (*(TYPE*)(VAR))

#define A _ACCESS(double,a)
#define B _ACCESS(double,b)
int double_cmp(const void* a, const void* b) {
  return ( (2*(A>B)-1)*(A!=B) );
}
#undef A
#undef B

#define A _ACCESS(float,a)
#define B _ACCESS(float,b)
int float_cmp(const void * a, const void *b)  {
  return ( (2*(A>B)-1)*(A!=B)  );
}
#undef A
#undef B

#define A _ACCESS(int,a)
#define B _ACCESS(int,b)
int int_cmp(const void * a, const void *b)  {
  return(A - B); 
}
#undef A
#undef B

int fill(char *text, const char * pre, const char * breakon, int width) { 
  /*  fills TEXT to WIDTH, line-breaking after characters in BREAKON (on
   *   punctuation if BREAK == "")
   *   prefixing the lines with PRE. TEXT is modified, and the nr of lines
   *   is returned. TEXT will be generally be larger after a fill
   */
  int len, prelen, nlines=0, newlen, cl;
  char *new, *np, *tp;
  const char *br;

  len=strlen(text);
  newlen=2*len;				/* first guess */
  new=(char*)ALLOCA(newlen);

  prelen=strlen(pre);

  br= (breakon[0])? breakon : " \t,.;:-";

  cl=0;
  np=new;
  for(tp=text; *tp; tp++) {		/* step through text */
    if(cl==0) {
      strcpy(np, pre);			/* start with pre-string */
      cl=prelen; np += prelen;
    }
    if (*tp=='\n' || *tp == ' ') {
      *np = ' ';			/* replace newlines by ' ' */
      if (np[-1] != ' ')		/* but not if previous also ' ' */
	np++, cl++;
    } else {
      *np=*tp;
      np++, cl++;
    }
    if (cl >= width  && strchr(br, *tp)) { /* time for a \n */
      *np='\n'; np++; cl=0; nlines++;
    }
  } /*  */
  if (np[-1]!='\n')
    *np++ ='\n';
  *np='\0';
  strcpy(text, new);
  AFREEA(new);
  return nlines+1;		/* i think */
} /* fill */

int random_range(int n, int *array) { 
  /* write a random range of 0 .. N-1, into ARRAY of length N */
  int i, j, k, ran, left;

  for (i=0; i<n; i++)
    array[i]= -1;

  for (i=0,left=n; i<n; i++, left--) {
    ran=(int)left*(random()/(RAND_MAX+1.0));	/* between 0 and left */
    /* put number at ran-th empty place: */
    for (j=0; array[j] >= 0; j++)
      ;					/* find index 0: first available */
    for (k=0; k<ran; k+=(array[j] < 0) )  /* go to index `ran' */
      j++;
    array[j]=i;
  }
  return 0;
} /* random_range */

int randomize(int n, int size, void *array) {
  /* randomizes, in place, ARRAY of N things of SIZE */
  int i;
  int *new_order;
  char *ran_array;

  new_order=(int*)ALLOCA(n*sizeof(int));
  random_range(n, new_order);

  ran_array=(char*)ALLOCA(n*size);
  for (i=0; i<n; i++)
    memcpy(ran_array + new_order[i]*size, (char*)array + i*size, size);

  memcpy(array, ran_array, size*n);

  AFREEA(new_order);
  AFREEA(ran_array);
  return 0;
} /* randomize */

#define ARRAY_LOC(X) (void*)local.qrank_array+(*(int*)X)*local.qrank_size
static int qrank_compare(const void *a, const void *b) {
  /* wrapper around user-provided comparison function, so she can use the
   * same comparison-function as in qsort() */
  return (*local.qrank_cmp)(ARRAY_LOC(a), ARRAY_LOC(b)); /* that's all */
}
#undef ARRAY_LOC

void qrank(const void *array, int nelts, int size, 
	   int (*compare)(const void *, const void *),
	   int * indices, int * ranks) {  
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
   *   you might consider one of the following:
   * for (i=0; i<nelts; i++) {
   *   sorted_array[ i ] = array[ indices[i] ];
   *       /# OR #/
   *   sorted_array[ ranks[i] ] = array[i];
   * }
   */
  int i, *inds;

  if (indices)
    inds=indices;
  else {				/* indices not wanted, so we need */
    if (!ranks)				/* a scratch array here */
      return (void)warn("qrank: both INDICES and RANKS are NULL\n");
    inds=(int*)ALLOCA(nelts*sizeof(int));
  }  
  for (i=0; i<nelts; i++)
    inds[i]=i;
  
  local.qrank_cmp=compare;		/* make visible to sort function */
  local.qrank_array=array;
  local.qrank_size=size;

  qsort(inds, nelts, sizeof(int), qrank_compare);

  if(ranks)				/* if at al wanted */
    for (i=0; i<nelts; i++)
      ranks[ inds[i] ] =i;

  if(! indices)				/* clean up scratch array */
    AFREEA(inds);
  return;				/* or return rank correlation ? */
} /* qrank  */


#include <unistd.h>
void debug_pause(void) { 
  fprintf(get_err_file(),"\
debugging ... pid = %d ... suspend me ... you have 2 secs\n",
	      (int)getpid());
  sleep(2); 
  return;
}

char * findfile(const char * path, const char * name, const char * mode, 
       FILE ** filep) { 
  int l=0, c=0;
  char  *dir, *cp, *dirs, *sav;
  static char filename[ FILENAME_MAX ];
  FILE * file=NULL;
  
  /* default path is ".:<home dir>" : */

#define DEFAULTPATH 	cp=getenv("HOME");\
  		    	l = strlen(cp);\
			dirs=(char*)ALLOCA(l+1+2);\
			memcpy(dirs, ".:", 2);\
			memcpy(dirs+2, cp, l+1);

  
  if (!path || !path[0] || strchr(name, '/') ) /* NULL or "" */
    { DEFAULTPATH }			/* $cwd and $HOME */
  else { 
    cp=getenv(path);
    if(cp) {				/* is path an environment var ? */
      l = strlen(cp);
      if(! l)				/* empty env var */
	{ DEFAULTPATH }			/* $cwd and $HOME */
      else				/* yes, not empty */
	MEMSAV(dirs, cp, l);
    } else {				/* not an env var */
      l=strlen(path);
      MEMSAV(dirs, path, l);		/* but itself a list */
    }
  }
  
  /* save a copy of dirs, so we can waste it */
  MEMSAV(sav, dirs, strlen(dirs) );
  cp=dirs;				/* initialize search */
  while( (dir=strtok(cp, ":")) ) {	/* orig destroys 'dirs', but we */
    cp=NULL;				/* still have 'sav' */
    l=strlen(dir);
    if (dir[l-1] == '/')		/* sometimes happens */
      dir[l-1]='\0';			/* delete that */
    if (strcmp(dir, "."))		/* i.e. if not "." */
      sprintf(filename, "%s/%s", dir, name); /* dir/file */
    else
      strcpy(filename, name);		/* ./file: use just file to allow */
    /* normal unix file naming, both absolute and relative */
    
    if (strlen(filename) > FILENAME_MAX )
      DIE("file name too long:%s", filename);
    file=fopen(filename, mode);
    if (!file)
      continue;

    c=fgetc(file);
    if(c == EOF)
      fclose(file);			/* will continue; */
    else
      break;				/* ready: found a file */
  } /* while strtok */

  if(!file)				/* found nothing */
    warn("file %s not found in %s\n", name, sav);

  AFREEA(dirs);
  AFREEA(sav);

  if(!file) {				/* found nothing */
    if (filep)
      *filep=NULL;
    return NULL;			/* "" */
  }

  if (filep) {				/* file pointer wanted */
    ungetc(c, file);
    *filep=file;
  } else
    fclose(file);

  return filename ;
} /* findfile */


void * _readfile(FILE * file, const char *_file, int line) { 
  /* read FILE into dynamically allocated memory, and return it */
  int n;
  size_t s;
  struct stat buf;
  void * vp;

  if ( fstat( fileno(file), &buf) ) { /* also get stat on file */
    perror("");
    return NULL;
  }

  s=(size_t)buf.st_size;
  vp= _malloc(s, _file, line);
  n=fread(vp, 1, s,  file);		/* perform the read */
  if(n != s) {
    warn("read %d characters of exspected %d\n", n, s);
    perror("");
    return NULL;
  }
  return vp;
} /* readfile */

void * _readfilename(const char *path, char * name, 
		     const char *_file, int line) { 
  char * found;
  FILE * file;
  void * vp;

  found=findfile(path, name,  "r", &file);
  if (!found) {
    perror(name);
    return NULL;
  }
  vp = _readfile(file, _file, line);
  if (!vp)				/* something wrong; perror() called */
    warn(name);				/* already, now give name */
  return vp;
} /* readfile */
