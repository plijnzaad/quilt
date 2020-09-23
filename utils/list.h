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

/* Lisp stuff. Only cons & extend allocate; others mess with existing lists */

#ifndef _LIST_H_
#define _LIST_H_ 1
#include "compat.h"

#define ADDLINK(ROOT,LINK, NEXT,TYPE) \
((LINK)=(*((ROOT)?&((LINK)->NEXT):&(ROOT)))=(TYPE *)CALLOC(1, sizeof(TYPE)))
/* E.G:
 *   while(i<=10) {
 *     link=ADDLINK(root, link, follow, struct _link);
 *     link->val=i++;			
 */


typedef void * car_p;
typedef struct list_ {			/* generic linked list */
  struct list_ * next;
  car_p ptr;
} list_t, cons_t, 			/* synonyms */
  *list_p, *cons_p;

/* for lisp afficionados: */
#define car 	ptr
#define cdr 	next
#define caar 	car->car 
#define caaar 	caar->car
#define cddr  	cdr->cdr
#define cdddr 	cddr->cdr
#define cadr  	cdr->car		/* == (car (cdr l)) */
#define cdar	cdr->car		/* == (cdr (car l)) */
/*   etc ... */

#define Lget(l, str, f) (((struct str *)(l->ptr))->f)
  /* meaning: consider l's car to be a struct str*, and get its field f  */
  /* may be assigned to, of course  */

#define Lcons(X,Y) _Lcons(__FILE__, __LINE__, (X), (Y))
list_p _Lcons(const char *file, int line, 
	      void * car, list_p cdr);	/* does (cons car cdr); list returned*/
list_p Lextend(list_p cdr, void * pointer); /* cons at end */

void * _Lnlist(int n, int elsize, const char *file, int line);
  /* returns linked list of length N, of elements of size ELSIZE. The
   * elements are bzero()ed */
#define Lnlist(X,Y) _Lnlist((X), (Y), __FILE__, __LINE__)

list_p Lcopy(list_p list);		/* copy of the of top-level cells */
list_p Lrcopy(list_p list);		/* same as Lcopy, but reversed list */

#define Lprepend(car, cdr) (car)->next=(cdr),(cdr)=(car)
					/* (setq a (cons a 'b)) */
#define Lunshift Lprepend
#define Lshift(cdr) (cdr=cdr->next)

list_p Llast(list_p list);		/* returns last elt of list */
int Llength(list_p list);		/* length of list */
  /* 'cons at end'; last elt returned */

#define Lfree(X, F) _Lfree((X), (F), __FILE__, __LINE__ )
int _Lfree(list_p list, void(*free_func)(), 
	   const char * file, int line); /*  */
  /* frees all cons-cells of list, and if free_func!=NULL, also cars, unless */
  /* they're NULL; free_func should be the address of the function to be */
  /* called to free the cars (which might even be Lfree 8-)*/ 

#define FREE_CARCONS(x) do {if((x)->car)FREE(x->car); FREE(x);} while(0)
                                       /* free both car and cons-cell */

list_p Lreverse(list_p list, int size);	/* original list destroyed */

list_p Lsort(list_p list, int (*compare)(const void * one, const void * other));
  /* sorts list in place; original order is lost. The compare function  */
  /* with type 'int (*compare)(list_p * one, list_p * other)', compares */
  /* pointers to list_p's, not list_p's themselves ! */
list_p Lrandomize(list_p);
  /* randomizes list in place, and returns it */

list_p *Lpreceed(list_p *listp, list_p ELT); 
/* returns adress whose contents are ELT list; searches in listp; NULL if */
/* not found. NB: searching for a NULL list will return  a pointer to */
/* '.next' of last element ! */
list_p Lexcise(list_p *list, list_p ELT);
/* excises ELT out of list, (rewires listp), and returns ELT (or NULL if */
/* not found); ELT ->next is not NULLed */

list_p Lrotate(list_p list, int n);	/* Lrotate(1,(a b c d)) = (b c d a) */

list_p Lnthcdr(int N, list_p);
/* returns (nthcdr N LIST), i.e. list starting at nth elt. of list (or NULL)*/
list_p Lsplit(int N, list_p);
/* like (nthcdr N LIST), but LIST itself truncated to length N */

list_p Lnconc(int nlists, list_p list, ...);
/* puts NLISTS lists together like nconc. Terminate list with ENDLIST */

#endif /* _LIST_H_ */
