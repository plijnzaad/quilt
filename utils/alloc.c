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

/* alloc/dynamic memory watch */

/* module for doing allocation that logs itself if needed, and dies */
/* with decent error message if out of memory */

#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include "utils.h"
#include "alloc.h"

int dmw_on;				/* dynamic control of dmw; set */
					/* from debugger */
static int _dmw;			/* switch; static control: which */
					/* sections of program are to be */
					/* controlled. Not touched by user! */
static FILE  * _dmw_file;		/* also no touched. Usage DOdmw(0 */
					/* or 1) to overall switch on or off */

#define DECEASE(S)  __file__=filename,__line__=linenr,die(S)

void * _malloc(int size, const char *filename, int linenr) {
  int * ip; 
  
  if (size) {
    ip=malloc(size);
    if (_dmw)
      fprintf(_dmw_file,"%s:%d: malloc(%d) returned 0x%lx\n",
	      filename,linenr,size, (ulong)ip);
    if(ip)
      return(ip);
  }
  die("malloc(%d) failed at %s %d\n", size, filename, linenr);
  return NULL;
} /* _malloc */

void * _calloc(int nelem, int elsize, 
	       const char *filename, int linenr)  {
  int * ip;

  if (nelem * elsize) { 
    ip=calloc(nelem, elsize);
    if (_dmw)
      fprintf(_dmw_file,"%s:%d: calloc(%d,%d) returned 0x%lx\n",
	      filename,linenr,nelem, elsize, (ulong)ip);
    if(ip)
      return(ip);
  }
  die("calloc(%d, %d) failed at %s %d\n", nelem, elsize, filename, linenr);
  return NULL;
} /* _calloc */

void * _realloc(void * ptr, int size, const char *filename, int linenr) {
  int * ip;

  if(size) { 
    ip=realloc(ptr, size);
    
    if (_dmw)
      fprintf(_dmw_file,"%s:%d: realloc(0x%lx,%d) returned 0x%lx\n",
	      filename,linenr,(ulong)ptr, size, (ulong)ip);
    if(ip)
      return(ip);
  }
  die("realloc(0x%lx, %d) failed at %s %d\n", ptr, size, filename, linenr);
  return NULL;
} /* _realloc */

void * arealloc( void*ptr, int size, const char *filename, int linenr) { 
  int * ip;

  if(size) { 
    if (ptr==NULL)
      ip=_calloc( 1, size, filename, linenr);
    else
      ip=_realloc(ptr, size, filename, linenr);
    
    if(ip)
      return(ip);
  }
  die("arealloc(0x%x, %d) failed at %s %d\n", ptr, size, filename, linenr);
  return NULL;
} /*  arealloc */

void _free(void * ptr, const char *filename, int linenr) {
  if (_dmw)
    fprintf(_dmw_file,"%s:%d: free() given 0x%lx\n",
	    filename,linenr,(ulong)ptr);
  if (!ptr)
    die("NULL pointer given to free() at %s %d\n", filename, linenr);
  free(ptr);
} /* _free */

void _dodmw(int switch_on, const char *filename, int linenr) { 
  /* used to open/close the log file */
  time_t t;

  if (!dmw_on)
    return;

  _dmw= switch_on;			/* for use in the _XXX functions */

  if (switch_on) {			/* start or resume logging */
    if ( !_dmw_file) {			/* file not yet opened */
      if_not (_dmw_file=fopen(DMW_LOG, "w"))
	ABORT("unable to open file '%s' in current directory", DMW_LOG);
      if( setvbuf(_dmw_file, NULL, _IOLBF, 133) /* line-buffer it */)
	ABORT("unable to setvbuf");
      time(&t);
      fprintf(_dmw_file, "# starting log at %s %d, %s" /* no '\n' ! */
	      , filename, linenr, ctime(&t));
    } else
      fprintf(_dmw_file, "# resuming log at %s %d\n", filename, linenr);
  }
  else {				/* suspend logging */
    if (_dmw_file) {			/* may be switched off several times */
      fprintf(_dmw_file, "# suspending log at %s %d\n", filename, linenr);
      fclose(_dmw_file);
      _dmw_file=NULL;
    }
  }
} /* _dodmw */
