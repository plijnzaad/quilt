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

/* lisp stuff. Allocation error messages will refer to here, */
/* unfortunately. May change this one day. */

#include <stdarg.h>
#include <stdlib.h>

#ifndef _VA_LIST_
#define _VA_LIST_			/* for avoiding errors with gcc ?!!? */
#endif

#include "utils.h"			/* for MALLOC and args */
#include "alloc.h"
#include "list.h"

list_p _Lcons(const char * file, int line,
	      car_p ptr, list_p list) {	/* does (cons car cdr); list returned*/
  list_p l = (list_p)_malloc(sizeof(list_t), file, line);
  l->ptr=ptr;				/* may be NULL */
  l->next=list;				/* may be NULL */
  return l;
}

list_p Llast(list_p list) {		/* returns last elt of list */
  list_p l;

  for (l=list; l->next; l=l->next) 
    ;					/* that's it */
  return l;
}

int Llength(list_p list) {		/* returns length of list */
  int i=0;
  list_p l;

  for (l=list; l; l=l->next) 
    i++;				/* that's it */
  return i;
}

list_p Lextend(list_p list, car_p ptr){	/* 'cons at end'; last elt returned */
  list_p last, newlast;

  last=Llast(list);
  newlast=(list_p)CALLOC(1, sizeof(list_t));
  last->next=newlast;
  newlast->car=ptr;
  return newlast;
}

int _Lfree(list_p list, void(*free_func)(), 
	   const char* file, int line) { 
  /* frees all cons-cells of list, and if ptr!=NULL, also cars, unless */
  /* they're NULL */ 
  int i=0;
  list_p l,n;

  for (l=list; l; l=n) {
    i++;
    n=l->next;
    if (free_func != NULL && l->ptr != NULL)	/* free car */
      (*free_func)(l->ptr, file, line);
    _free(l, file,line);
  }
  return i;
}

void * _Lnlist(int n, int elsize, const char * file, int line) {
  /* returns linked list, contiguously allocated in one go */
  int i;
  list_p list, l;
  ulong next;

  list= (list_p)_calloc(n,elsize, file, line);

  for (i=0,l=list; i<n-1; i++, l=(void*)next) {
    next=((ulong)l)+elsize;
    l->next=(void*)next;
  }

  return list;
} /* Lnlist */

list_p Lreverse(list_p list, int size)  { /* original list destroyed */
  int i, len=Llength(list);
  list_p l, *array;
  
  if (len==0)
    return(NULL);
  if (len==1)
    return(list);

  array=(list_p*)ALLOCA(len*size);

  i=0;
  for (l=list; l; l=l->next)
    array[i++]=l;
  for (i=len-1; i>0; i--) 
    array[i]->next=array[i-1];
  array[0]->next=NULL;

  list=array[len-1];

  AFREEA(array);

  return list;
}

list_p Lsort(list_p list, 
	     int (*compare)(const void * one, const void * other)) { 
  int i, len=Llength(list);
  list_p l, *array;
  
  if (len < 2)
    return(list);

  array=(list_p*)ALLOCA(len*sizeof(list_p));

  i=0;
  for (l=list; l; l=l->next)
    array[i++]=l;
  
  qsort(array, len, sizeof(list_p), compare);

  for (i=0; i<len-1; i++) 
    array[i]->next=array[i+1];
  array[len-1]->next=NULL;
  list=array[0];

  AFREEA(array);

  return list;
} /* Lsort */

list_p * Lpreceed(list_p *listp, list_p l) {
  list_p k;

  assert(listp);			/* calling error otherwise */
  if(*listp == l)			/* searching for first elt */
    return listp;			/* (this can be NULL ...) */

  for (k= *listp;  k; k=k->next)
    if(k->next==l)			/* searching for a NULL list */
      return( & k->next);
  return(NULL);				/* simply not found */
} /* Lpreceed */

list_p Lexcise(list_p *listp, list_p l) {
  list_p *lp;
  
  lp=Lpreceed(listp, l);
  if (!lp)				/* not found */
    return NULL;
  assert( *lp == l);
  *lp=l->next;
  return l;
} /* Lexcise */


list_p Lnconc(int nargs, list_p first, ...) {
  /* puts NARGS lists together like nconc and returns result */
  int i;
  list_p list, last=NULL, l;
  va_list valist;

  list=first;				/* maybe NULL */
  if (list)
    last=Llast(list);			/* for the next time */
  va_start(valist, first); 
  nargs--;
  for (i=0; i<nargs; i++) {
    l = va_arg(valist, list_p);
    if (!l)				/* emtpy list */
      continue;
    if (list) 				/* if not the first one */
      last->next=l;
    else
      list=l;
    last=Llast(l);			/* for the next time */
  }
  va_end(valist);
  return(list);
} /* Lappend */

list_p Lrotate(list_p list, int n) {
  /* Lrotate(2, (a b c d e)) = (c d e a b) */
  list_p l, last;
  int len;

  if (n==0 || list==NULL)
    return list;

  len=0;				/* get length and last elt of list */
  for (l=list; l->next; l=l->next) 
    len++;
  len++;				/* since counted till penultimate ! */
  last=l;

  if (n<0)
    n += -n*len;			/* make positive */
  n %= len;				/* wrap around */
  if (n==0)				/* may still happen */
    return list;

  len=0;				/* use just as counter */
  for(l=list; l; l=l->next) {		/* walk to 'cut-point' */
    len++;
    if(len==n)				/* here it is */
      break;
  }
  last->next=list;			/* circularize */
  list=l->next;				/* new head of list */
  l->next=NULL;				/* put end to it */
  return list;
} /* Lrotate */

list_p Lcopy(list_p list) {
  /* copy of top-level cells */
  list_p l, copy, new, *prevnextp ;

  copy=NULL;
  prevnextp= &copy;
  for(l=list; l; l=l->next) {
    new=(list_p)MALLOC(sizeof(list_t));
    new->next=NULL;			/* 'cos could be last one */
    new->ptr=l->ptr;			/* copy of car */
    assert(*prevnextp==NULL);
    *prevnextp=new;
    prevnextp = &new->next;
  }
  return copy;
} /* Lcopy */

list_p Lrcopy(list_p list) {
  /* makes a reverse copy of list's top level cells */
  list_p copy, l;

  copy=NULL;
  for (l=list; l; l=l->next) 
    copy=Lcons(l->car, copy);
  return copy;
}

list_p Lrandomize(list_p list) {
  int i, len;
  list_p l, *array;

  len=Llength(list);
  if (len<2)
    return list;			/* NULL or single elt. */

  array=(list_p*)ALLOCA(len*sizeof(list_p));

  i=0;
  for (l=list; l; l=l->next)
    array[i++]=l;
  
  randomize(len, sizeof(list_p), array);

  for (i=0; i<len-1; i++) 
    array[i]->next=array[i+1];
  array[len-1]->next=NULL;
  list=array[0];

  AFREEA(array);

  return list;
} /* Lrandomize */

list_p Lnthcdr(int n, list_p list) {
  /* returns list starting at nth elt in list */
  int i;
  list_p l;

  for (l=list,i=0; l; l=l->next,i++)
    if (i==n)
      break;
  return l;
} /* Lnthcdr */

list_p Lsplit(int n, list_p list) {
  /* truncates list a element n, and returns tail */
  int i;
  list_p l, prev;

  prev=NULL;
  for (i=0,l=list; l; i++,prev=l,l=l->next)
    if (i==n)
      break;
  if (prev)
    prev->next=NULL;			/* doesn't hurt if at end of list */
  return l;
} /* Lsplit */
