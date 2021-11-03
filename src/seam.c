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
 * Philip Lijnzaad, plijnzaad@gmail.com   			      *
 * ------------------------------------------------------------------ */

/* routines for identifying a seam, and for triangulating the simple cases */

#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <limits.h>
#include "extra-math.h"

#include "utils.h"
#include "alloc.h" 
#include "list.h"
#include "vecmat.h" 

#include "quilt.h"


static con_p next_con(con_p con, int end, int *lenp) {
  /*
   * does    : returns the next (existing) con along a (clockwise) path
   *   along a seam; NULL if it is the same con again
   *
   * gets    : con to start from, and the id of the bndry at the other end
   *
   * affects : *bndry_p is set to the new bndry; *lenp is set to the
   *   length of the traversed stretch over the bndry. If the same con is
   *   found, the length is obviously that of the complete boundary
   *
   * returns : the newly found con, which may actually be the same one, if
   *   patch is connected by only one boundary
   *
   * warns   : 
   *
   * comment : 
   */
  int i, j, offset, length;
  bndry_p bn;
  list_p * cons, l;

  bn=con->bndries[ end ];			/* the other bndry */
  offset=con->flags[ end ].offset;
  cons=bn->cons;
  length=bn->length;
  
  for(l=cons[offset]; l->next; l=l->next) /* if l->next, more than one */
    if (l->ptr==con) {			/* edge on point, so walk until */
      *lenp=1;				/* current found, and return next */
      return (con_p)l->next->ptr;	/* If fails, walk along bndry */
    }
  if (length==1) {			/* isolated point: apparently  */
    *lenp=1;				/* started at end of list, */
    return cons[offset]->ptr;		/* so circularize */
  }
  /* else: walk along bndry */
  for (i=(offset+1)%length,j=1; j<=length; i=(i+1)%length,j++)
    if(cons[i])	{			/* a (sorted!) list of cons: take */
      *lenp=j+1;				/* the first one */
      return (con_p)cons[i]->ptr;	/* may be equals con ! */
    }
  assert(0);				/* can not happen */
  return NULL;
  /* not found: patch connected to the rest by just one con. return this */
  /* con, to make sure it doesn't count as a cross-over */
  *lenp=length;
  return con;
} /* next_con */

static int seam_loop(con_p firstcon, int firstend, int *nseamsp, 
		     int maxnseams, seam_p seams) {
  int ncross, nseams, end, length;
  con_p con, newcon;
  bndry_p bndry;
  SHORT *first, *last;


#define THING (((shell_p)bndry->patch->atom->pointer)->pointflags + \
  bndry->path[ newcon->flags[ end ].offset ] )
/* #define THING ( bndry->path + newcon->flags[ end ].offset ) */

  end = firstend;
  newcon=firstcon;
  bndry=newcon->bndries[ end ];		/* the bndry following the con */
  first=THING;

  ncross= *nseamsp = nseams=0;

  do {
    con=newcon;
    con->flags[ end ].seen = 1;
    if (con->flags[ ! end ].seen) {	/* taken in opposite direction */
      assert(nseams>0);			/* must have been in this call */
      ncross++;				/* earlier: edges cross */
    }

    newcon=next_con(con, end, &length);	/* find new one, and length of old */
    last= ((shell_p)bndry->patch->atom->pointer)->pointflags +
      bndry->path[ (con->flags[end].offset + length -1 + bndry->length)
		  % bndry->length ];

    /* first save data on previous seam_t: */
    seams[ nseams ].end=end;
    seams[ nseams ].bndry=bndry;
    seams[ nseams ].first=first;
    seams[ nseams ].last=last;
    seams[ nseams ].con=con;
    seams[ nseams ].length=length;
    nseams++;
    if (nseams >= maxnseams/2) {
      if (nseams >= maxnseams)
	  DIE("seam too long; could not accomplish triangulation, ");
      if(ncross) {			/* clean up first !!! */
	*nseamsp = nseams;
	return 4;
      }
    }

    end = (newcon->bndries[0] == bndry); /* new end */

    bndry = newcon->bndries[ end ];	/* new bndry */
    first = THING ;

/*    assert(newcon != con);		/# not allowed anymore */

  } while_not (newcon == firstcon && end==firstend );
  *nseamsp=nseams;

  assert(ncross < maxnseams );
  assert(nseams>1 && nseams < maxnseams );

  if (ncross > 0)
    return 4;
  return 0;				/* normal exit */
} /* seam_loop */

int find_seam(con_p firstcon, int firstend, int *nseamsp, 
	      int maxnseams, seam_p seams) { 
  /*
   * does    : assembles a seam by traveling along cons constituting it.
   *
   * gets    : a starting con + direction, plus place to leave results.
   *
   * affects : writes results (con, end, bndry, from, to) into
   *   seams[*]. Also sets the con->flags[end].seen flags on traversed
   *   cons.  *nseamsp is set to the number of seams found. 
   *
   * returns : 0  seam constitutes complete triangle
   *           1  normal seam, length >= 2
   *           2  abnormal seam consisting of one connection (length == 2)
   *           3  abnormal seam, bndry crossings
   *           4  abnormal seam, con crossings
   * warns   : 
   *
   * comment : 
   */
  int i, length, nseams;

  i=seam_loop(firstcon, firstend, nseamsp, maxnseams, seams);

  if(i) {
    assert( i==4);
    return i;				/* crossing */
  }

  nseams= *nseamsp;

  length=0;
  for (i=0; i<nseams; i++) {
    assert(seams[i].length >= 1);
    length += seams[i].length;
  }
  assert(length>=2);
  if (length==2)			/* abnormal: just one connection */
    return 2;

  if(length==3) {			/* found a triangle */
    assert(nseams <= 3);		/* since length always >= 1 */
    for (i=0; i<nseams; i++) {
      seams[i].con->flags[ seams[i].end ].done =1; 
      seams[i].con->flags[ seams[i].end ].seen =0;
    }
    return 0;				/* complete triangle exit */
  }
  return 1;				/* normal exit: has to be subdiv'd */
} /* find_seam */


