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

/* routines for establishing connectivity of atomic patches, and sorting */
/* used for the full triangulation  (con.c (1 n) is used for patch thing */

#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <limits.h>
#include "extra-math.h"

#include "utils.h"
#include "pdb.h"
#include "alloc.h" 
#include "list.h"
#include "vecmat.h" 

#include "quilt.h"

/* program global datum */

typedef struct geom_ { 
  shell_p shella, shellb;
  float *coorda, *coordb, *xa, *xb;
  double * saxis;
} geom_t, *geom_p;

static double calc_delta(double d2, double ra, double rb, 
			sphere_p spha, sphere_p sphb ) { 
  /* calculates a scaled squared distance that serves as cutoff for */
  /* allowing a connection between boundary points on both edges */
  double a, b, phi, cosphi, c, halfalpha, beta, ra2, psi
    , rb2, x, y, delta1, delta2;

  a=spha->edgelengths[1]*ra;
  b=sphb->edgelengths[1]*rb;

  halfalpha=spha->halfalpha;
  beta=sphb->beta;

  if(spha!=sphb) {			/* choose largest pair */
    x=spha->beta;
    y=sphb->halfalpha;			/* other way around */
    if (halfalpha + beta < x+y) {
      halfalpha=x;
      beta=y;
    }
  }
  
  ra2=SQR(ra);
  rb2=SQR(rb);
  c = (ra2 + rb2 -d2)/(2.0*ra*rb);	/* this is cos(pi-psi), not cos(psi) */
  /* (psi is angle between tangent planes at intersection) */

  /* max angle between two edges: */
  psi= M_PI - ACOS(c);			/* cannot avoid this one */
  /* alternative would be to use a sqrt, and tabulate (halfalpha + beta)'s */
  phi= psi + halfalpha + beta; 

#undef DEFINED
#ifdef DEFINED
  if (phi < M_PI_2) {			/* yo, save a cosine */
    if(a>b)
      delta = SQR(a)+0.25*SQR(b);
    else
      delta = SQR(b)+0.25*SQR(a);
  } else { 
    cosphi = cos(phi);			/* could trade for a sqrt + table */
    delta = SQR(a) + SQR(b) - M_SQRT3 * a*b *cosphi;
  }

  return delta/rb2;
#else

  if(a>b)
    delta1 = SQR(a)+0.25*SQR(b);
  else
    delta1 = SQR(b)+0.25*SQR(a);
  cosphi = cos(phi);			/* could trade for a sqrt + table */
  delta2 = SQR(a) + SQR(b) - M_SQRT3 * a*b *cosphi;
/*   printf("! %f %f %f\n", phi*57.3, delta1, delta2); */
  return MAX(delta1,delta2)/rb2;
#endif
} /* calc_delta */

#ifdef undefined
static void set_midpoint(con_p con, geom_p geom, double *midpoint) {
  int i, p;
  const float *xuspt, *xatom, *dir;
  double r, len;
  bndry_p bn;

  bn=con->bndries[0];
  p=bn->path[ con->flags[0].offset];
  xuspt=geom->coorda + 3*p;		/* unit sphere coords */
  xatom=geom->xa;
  r=xatom[3];				/* radius; */
  dir=con->direction;
  len=dir[3];
  for (i=0; i<3; i++)
    midpoint[i]=  xatom[i] + r*xuspt[i] + 0.5*len*dir[i];
  return;
} /* set_midpoint */

static int is_buried(shell_p shell, double *midpoint) {
  /* returns 1 if midpoint is buried by any of shell's neigbours, else 0 */
  int i, k;
  shell_p sh;
  const float *xyz;
  double d,d2, r2;

  for (i=0; i<shell->nngbs; i++) {
    sh=shell->ngbs[i];
    xyz= sh->atom->xyz;
    r2=SQR(xyz[3]);
    d2=0.0;
    for (k=0; k<3; k++) { 
      d=xyz[k] - midpoint[k];
      d2 += SQR(d);
    }
    if (d2 < r2 )
      return 1;
  } /* for shp */
  return 0;
} /* is_buried */

static int midpoint_buried(con_p con, geom_p geom) { 
  /* returns 1 if midpoint of this connection is buried by any of bndries */
  /* ngbs; otherwise 0 */
  double midpoint[3];

  set_midpoint(con, geom, midpoint);
  return (is_buried(geom->shella, midpoint) || 
	  is_buried(geom->shellb, midpoint) );
} /* midpoint_buried */
#endif /* undefined */

static int keep_con(con_p prelim_cons, 
		bndry_p a, int offsa, bndry_p b, int offsb,
		double *dir, double d, geom_p geom) {
  /* fills in further relevant stuff about the con into prelim_cons */
  int k;
  con_p con;
  double rb;


  con= &prelim_cons[offsa];		/* lvalue */
  con->bndries[0]=a;
  con->flags[0].offset=offsa;		/* will this be truncated ???? */
  con->bndries[1]=b;
  con->flags[1].offset=offsb;

  rb=geom->xb[3];			/* radius of b atom */
  for (k=0; k<3; k++) 
    con->direction[k]= dir[k]*rb/d;	/* unscaled, normalized direction */
  con->direction[3]= d;

#ifdef undefined
  if ( midpoint_buried(con, geom) ) {
    warn("midpoint buried\n");
    return 0;
  }
#endif
  return 1;
} /* keep_con */

static int NConsAllocated, NConsDeleted, NConsTotal;

con_p alloc_con(const con_t *con) {
  con_p newcon;

  NConsTotal++;
  NConsAllocated++;
  newcon=(con_p)malloc(sizeof(con_t));
  if (con)
    *newcon=*con;
  else
    bzero(newcon, sizeof(con_t));
  return newcon;
} /* _alloc_con */

void free_con(con_p con) {
  NConsTotal--;
  NConsDeleted++;
  FREE(con);
  return;
}

con_p save_con(bndry_p a, int offsa, bndry_p b, int offsb,
	       double *dir, double d2, double rb) {
  /* 
   * does    : allocates space to hold a con, and fills in relevant stuff
   *   of one or both paths involved
   *   
   * gets    : one bndry+offset, and other bndry+offset, plus geometric data
   *   dir should be {x,y,z}/rb; d2 is SQR(dist/rb); rb should be 1.0 if
   *   unscaled data is given
   *
   * affects : nothing
   *
   * returns : the new con. 
   *   
   *
   * comment : con->bndries[0]->patch->atom garantueed to be 
   *   smaller than con->bndries[1]->patch->atom
   */
  int k;
  con_p con;
  double d;

  con=alloc_con(NULL);
  con->bndries[0]=a;
  con->flags[0].offset=offsa;		/* will this be truncated ???? */
  con->bndries[1]=b;
  con->flags[1].offset=offsb;

  d=sqrt(d2);				/* this is the scaled distance */
  for (k=0; k<3; k++) 
    con->direction[k]= dir[k]/d;	/* unscaled, normalized direction */
  con->direction[3]= d*rb;		/* unsquared, unscaled distance */

  return con;
} /* save_con */


static int decr_bndry_con(bndry_p bndry) {
  /* decrements nr of cons of bn, and takes action if it was the last one; */
  /* retuns 1 if it was, 0 otherwise */

  bndry->ncons--;
  if ( bndry->ncons > 0 )
    return 0;				/* normal exit: still cons left */

  bndry->flags |= ISWIDOW;
  ((shell_p)bndry->patch->atom->pointer)->flags |= REMOVE_WIDOWS;
  return 1;
} /* decr_bndry_con */


static int delete_conref(con_p con, int end) {
  /* deletes list_p that references con at side 'end', excising it out of */
  /* linked list; returns 1 if remaining list is empty; 0 otherwise */
  int offset;
  bndry_p bn;
  list_p l, prev;

  bn=con->bndries[end];
  offset=con->flags[end].offset;
  prev=NULL;
  for(l=bn->cons[ offset ]; l; prev=l,l=l->next)
    if(l->ptr==con)			/* (fast, and happens seldom) */
      break;

  assert(l);				/* otherwise not found */

  if(prev)				/* middle or end of list */
    prev->next=l->next;
  else					/* begin of list */
    bn->cons[ offset ]=l->next;
  FREE(l);

  return decr_bndry_con(bn);
} /* delete_conref */

int remove_con(con_p con) { 
  /* deletes conrefs and con; returns 0 if nothing happend, 1 if bndry was */
  /* deleted, 2 if patch, 4 if shell; if happend to 2nd, returns patterns */
  /* shifted left by 3 */
  int i;

  i=delete_conref(con, 0) + 2*delete_conref(con,1);
  free_con(con);			/* con itself */
  return i;
} /* remove_con */

int remove_con_halfknown(list_p l, bndry_p bndry) { 
  /* same as remove_con, but faster if conref one side (bndry) already knonw */
  int i,end;
  con_p con;
  
  con=l->ptr;
  end= (con->bndries[0]==bndry);	/* other end */
  FREE(l);
  i= decr_bndry_con(bndry)		/* this end */
    + 2*delete_conref(con, end);	/* other end */
  free_con(con);			/* con itself */
  return i;
} /* remove_con_halfknown */

static void mark_longest(bndry_p bndry, int offset) {
  /* marks longest con at one offset */
  int end, maxend;
  list_p l, cons;
  con_p con, maxcon;
  float max;

  cons=bndry->cons[offset];
  
  if(cons == NULL || cons->next==NULL )
    return;

  maxcon=NULL;
  max= -10.0;
  for (l=cons; l; l=l->next) {
    con=l->ptr;
    end=(con->bndries[0]==bndry);
    if(con->direction[3] > max) {
      max=con->direction[3];
      maxcon=con;
      maxend=end;
     }
  }
  assert(maxcon);
  maxcon->flags[maxend].seen=1;
  return;
} /* mark_longest */

static int purge_longest(bndry_p bndry, int offset) {
  /* deletes longest con at offset, if it's also longest at other bndry */
  int i,j;
  list_p l, next, cons;
  con_p con;

  cons=bndry->cons[offset];
  
  if(cons==NULL)
    return 0;

  for (l=cons,i=j=0; l; l=l->next,i++) {
    con=l->ptr;
    if(con->flags[0].seen ==1 && con->flags[1].seen ==1) {
      j++;
      con->flags[0].done=1;		/* tag it */
    }
  } /* for i */
  assert(j<i);

  if(j)					/* have to delete some */
    for (l=cons; l; l=next) {		/* delete everything marked */
      next=l->next;
      con=l->ptr;
      if(con->flags[0].done==1)
	remove_con(con);
      else {			/* reset the rest */
	con->flags[0].done=0;
	con->flags[1].done=0;
	con->flags[0].seen=0;
	con->flags[1].seen=0;
      }
    }	/* for l */
  else					/* nothing to be deleted: reset */
    for (l=cons; l; l=l->next) {
      con=l->ptr;
      con->flags[0].done=0;
      con->flags[1].done=0;
      con->flags[0].seen=0;
      con->flags[1].seen=0;
    }
  return j;				/* n of cons deleted */
} /* purge_longest */

static int clean_conns(box_p box) {
  /* if both end points of a connection have more than one connection, and */
  /* a connection is the longest among them in BOTH endpoints, then that */
  /* connection is deleted by this routine. Serves to clean up potentially */
  /* crossing edges */
  int i, j, deleted;
  shell_p shell;
  patch_p pa;
  bndry_p bn;

  deleted=0;
  for (i=0,shell=box->shells; i<box->natoms; i++,shell++) {
    if ( remove_widows(shell, &deleted) )
      continue;
    for(pa=shell->patches; pa; pa=pa->next)
      for(bn=pa->bndries; bn; bn=bn->next)
	for (j=0; j<bn->length; j++) 
	    mark_longest(bn, j);
  } /* for i */

  for (i=0,shell=box->shells; i<box->natoms; i++,shell++) { 
    if ( remove_widows(shell, &deleted) )
      continue;
    for(pa=shell->patches; pa; pa=pa->next)
      for(bn=pa->bndries; bn; bn=bn->next)
	for (j=0; j<bn->length; j++)
	  deleted += purge_longest(bn, j);
    remove_widows(shell, &deleted);
  } /* for shell */
  return deleted;
} /* clean_conns */

static double phi_of(bndry_p bndry, con_p con, 
		     const double *x0, const double *y0) {
  /*
   * does    : calculates angle of a con with respect to x0 and y0
   *
   * gets    : a con, a bndry (for finding out sense), and (orthonormal)xy axis
   *
   * affects : *dotp is set to the dotproduct of connection vector with
   *   surface normal (and should in fact be > 0 or so)
   *
   * returns : angle in range 0, 2bpi
   *
   * warns   : about cons with z < 0 and 
   *
   * comment : 
   */
  int k;
  double x, y, phi, sign;
  const float *dir;

  sign=2.0*(con->bndries[0]==bndry) - 1.0;
  dir=con->direction;

/*   dot=DOT3P(n, dir);
 *   for (k=0; k<3; k++ )		/# get projection onto xy plane #/
 *     p[k]=sign*(dir[k]-dot*n0[k]);
 *   *dotp = sign*dot;
 */

  x=y=0.0;
  for (k=0; k<3; k++ ) { 
    x += dir[k]*x0[k];
    y += dir[k]*y0[k];
  }	  

  phi=atan2( sign*y, sign*x);
  if (phi<0.0)
    phi+=M_2PI;

  return phi;
} /* phi_of */

static list_p * con_place_normal(con_p toadd, bndry_p bndry, int offset, 
				 const double * xaxis) { 
  /*
   * does    : finds place for toadd in list bndry->cons[offset], sorted
   *   with respect to normal and xaxis. Called for UTURN's and ISOLATED
   *
   * gets    : list to work from, etc.
   *
   * affects : list is taken apart, and inserted into bndry->cons[offset]*
   *
   * returns : NULL if thing connection already present,
   *   &predecessor->next otherwise
   *
   * warns   : about anlge being too small
   *
   * comment : if an error occurs (connection already present), the
   *   routine returns immediately, discarding the rest of the toadd list
   */
  int i, n, e1, e2, p;
  list_p l, *prevp;
  con_p lcon;
  double yaxis[3], phi, thisphi, dot;
  const float * normal;
  static double threshold=0.8;		/* wild guess */

  assert(toadd);
  assert( bndry->cons[offset] );

  p=bndry->path[offset];
  assert( ((shell_p)bndry->patch->atom->pointer)->pointflags[p]
	 & ( UTURN | ISOLATED) ); 
  normal= &((shell_p)bndry->patch->atom->pointer)->sphere->coordinates[3*p];

  e1=(toadd->bndries[0]==bndry);
  if(bndry->length==1) {		/* otherwise, it was already tested */
    dot= (2*e1 -1)*DOT3P(normal, toadd->direction);
    if(dot > threshold) {
      warn("angle too small in con_place_normal:");
      return NULL;
    }
  }

  CROSS3P(yaxis, normal, xaxis);

  n=0;					/* nr of existing cons of known phi */
  thisphi=phi_of(bndry, toadd, xaxis, yaxis); /* phi of new con */
  
  prevp = bndry->cons+offset;
  for (i=0,l= *prevp; l; prevp= &l->next,l= *prevp, i++) {
    lcon=l->ptr;
    e2=(lcon->bndries[0]==bndry);
    
    if (toadd->bndries[ e1 ] == lcon->bndries[ e2 ]
	&& toadd->flags[ e1 ].offset == lcon->flags[ e2 ].offset ) {
      warn("already present:");		/* occurs ??? */
      return NULL;
    }    
/*     if(i>=n)				/# existing phi not known #/
 *       phis[n++]=phi_of(bndry, lcon, xaxis, yaxis); /# calc, and store #/
 *     phi=phis[i];
 */
    phi=phi_of(bndry, lcon, xaxis, yaxis);

    if ( phi > thisphi )		/* first one that's larger */
      break;

/*       memmove( &phis[i+1], &phis[i], (n-i)*sizeof(double) ); /# shift #/
 *       phis[i]=thisphi;
 */

  } /* for l */

  return prevp;
} /* con_place_normal */

static int local_top(bndry_p bndry, int offset, 
		     double *pray, double *mray, double *nray) {
  /*
   * does    : establishes local topology at bndry->path[offset]; finds
   *   vectors from offset along path in two directions, and a vector
   *   halfway these two, pointing to the 'left' of the path
   *
   * gets    : 
   *
   * affects : *pray, *mray,*nray, which are set to previous, middle and
   *   next vector
   *
   * returns : 
   *
   * warns   : 
   *
   * comment : 
   */
  int i, m, p, pi, ni, li, pe, ne, me, prev, next, length, nngbs;
  SHORT *path, *ngbs;
  sphere_p sphere;
  edge_p pedge, medge, nedge;
  shell_p sh;
  point_p point;
  double sign;

  length=bndry->length;
  path=bndry->path;
  p=path[offset];
  sh=bndry->patch->atom->pointer;
  sphere=sh->sphere;
  point=sphere->Points+p;

  pi=ni= -1;
  prev = path [ (offset + length -1)%length ];
  next = path [ (offset + length +1)%length ];
  
  nngbs=point->nngbs;
  ngbs=point->ngbs;
  
  for (i=0; i<nngbs; i++) {
    if ( ngbs[i] == prev )
      pi=i; 
    if ( ngbs[i] == next )
      ni=i; 
  }
  assert(pi>=0 && ni>=0);
  
  /* find tangent edge in middle */
  li=(pi-ni+nngbs)%nngbs;		/* nr of edges between ne & pe */
  m=(pi - li/2 + nngbs)%nngbs;
  me=point->edges[m];

  medge=sphere->Edges+me;
  sign= 2.0*(medge->ends[0]==p) - 1.0;
  VEC3TIMASS(mray, sign, medge->direction);

  if (li & 1) {				/* this isn't enough yet */
    m=(m-1+nngbs)%nngbs;
    me=point->edges[m];
    medge=sphere->Edges+me;
    sign= 2.0*(medge->ends[0]==p) - 1.0;
    for (i=0; i<3; i++) 
      mray[i]=M_1SQRT3*(mray[i]+ sign*medge->direction[i]);
  }
  
  pe=point->edges[pi];
  ne=point->edges[ni];
  
  /* rest is easy */
  pedge=sphere->Edges+pe;
  sign= 2.0*(pedge->ends[0]==p) - 1.0;
  VEC3TIMASS(pray, sign, pedge->direction); /* edge on bndary */

  nedge=sphere->Edges+ne;
  sign= 2.0*(nedge->ends[0]==p) - 1.0;
  VEC3TIMASS(nray, sign, nedge->direction); /* can be same as pedge ... */
  
  assert( pe!=ne  || sh->pointflags[p] & UTURN );
  return (pe==ne) ;
} /* local_top */

static int con_present(con_p newcon, bndry_p bndry) {
  /*
   * does    : checks if the new con is indeed new to bndry at offset
   *
   * gets    : con to be checked for presence, bndry+offset, and end
   *
   * affects : nothing
   *
   * returns : 1 if present, 0 otherwise
   *
   * warns   : 
   *
   * comment : we don't have to search through all bndries to check if con
   *   is not there, as it will turn up in 2nd invocation anyway.
   *   
   */
  int end, offset;
  SHORT point;
  list_p l;
  con_p con;
  patch_p patch;
  bndry_p bn;

  end=(newcon->bndries[0]==bndry);	/* other end */
  bn=newcon->bndries[ end ];		/* partner bndry of newcon */
  assert(bn);
  offset=newcon->flags[ end ].offset;	/* offset in partner bndry */

  point=bn->path[offset];		/* ID of partner point: this is */
  patch=bn->patch;			/* the point that should be unique */

  offset=newcon->flags[ ! end ].offset;	/* offset in  bndry !!! */
  for(l=bndry->cons[offset]; l; l=l->next) {
    con=l->ptr;
    end=(con->bndries[0]==bndry);
    bn=con->bndries[end];
    assert(bn);

    if ( bn->patch == patch ) {		/* same patch: watch out  */
      offset = con->flags[ end ].offset;
      if (point==bn->path[offset])	/* and same (physical) point */
	return 1;
    }
  } /* for l */
  return 0;
} /* con_present */

list_p * con_place(con_p newcon, bndry_p bndry, int offset) {
  /*
   * does    : finds place for newcon in bndry->cons[offset], sorted so as to
   *   minimize the total curvature of the set of rays emanating out of point
   *
   * gets    : connection to be added
   *
   * affects : newcon is destroyed, bndry->cons[offset]* changes
   *
   * returns : address of predecessor->next; NULL if not found or invalid
   *
   * warns   : 
   *
   * comment : upon error, immediately returns, only excising offending
   *   cons from bndry and partner
   */

  int i, k, n, maxn, e1, e2;
  list_p l, *maxp, *prevp;
  con_p con;
  double sign, rays[3* MAX_NCONS ], max, dot, dot1, dot2, dot3, dot4
    , *a,*b,*c, ab,ad,bc,bd,cd, dray[3], lray[3], mray[3];
  static double threshold=0.8;		/* wild guess */

  assert(newcon);
  assert( bndry->cons[offset] );

  if(con_present(newcon, bndry)) { 
    warn("already present\n");
    return NULL;
  }
  if (bndry->length == 1) {		/* ISOLATED */
/*  assert(  ((shell_p)bndry->patch->atom->pointer)->pointflags[p] &ISOLATED);
 */
    con=bndry->cons[offset]->ptr;	/* first edge */
    sign=2.0*(con->bndries[0]==bndry) - 1.0;
    VEC3TIMASS(&rays[0], sign, con->direction);
    return con_place_normal(newcon, bndry, offset, &rays[0]);
  } /* if length==1 */

  e1=(newcon->bndries[0]==bndry);
  sign=2.0*e1 - 1.0;
  VEC3TIMASS(dray, sign, newcon->direction); /* dray is direction of con */
  
  i=local_top(bndry, offset, &rays[0], mray, lray);

  if(i)					/* a UTURN  point */
    return con_place_normal(newcon, bndry, offset, &rays[0]);

  dot=DOT3P(mray, dray);
  if(dot > threshold) {
    warn("angle too small in con_place:");
/*    pcon(NULL, newcon); */
    return NULL;
  }
  
  /* search closest ray: */    
  maxp= &bndry->cons[offset];		/* parameters of first bndry edge */
  max=DOT3P(dray, &rays[0]);
  maxn=0;

  n=1;					/* index of first con in array */
  prevp= &bndry->cons[offset];
  for(l=bndry->cons[offset]; l; prevp= &l->next,l= *prevp) { 
    con=l->ptr;
    e2=(con->bndries[0]==bndry);
    sign=2.0*e2 - 1.0;
    VEC3TIMASS(&rays[3*n], sign, con->direction);
    dot=DOT3P(dray, &rays[3*n]);
    if (dot > max) {
      maxp=prevp;
      max=dot; 
      maxn=n;
    }
    n++;
  } /* for l */
  assert(maxp);

  dot=DOT3P(dray, lray);		/* don't forget last one */
  if (dot > max) {
    maxp= &(Llast(bndry->cons[offset])->next); /* also easy: at end of list */
    maxn= -1;
  }
  if(maxn == 0 || maxn==-1)		/* ray is closer to one of bndries */
    return maxp;			/* easy; no further work needed */

  memcpy(&rays[3*n], lray, sizeof(double[3])); /* last ray, on bndry */

/* a   b
 * .___.  
 *    . \
 *   d   .c        d is new con to be added; b is the con it is closest to
 */        

  a= &rays[3*(maxn-1)];
  b= &rays[3*(maxn)];
  c= &rays[3*(maxn+1)];
  ab=ad=bc=bd=cd=0.0;
  for (k=0; k<3; k++) {
    ab += a[k]*b[k];
    ad += a[k]*dray[k];
    bc += b[k]*c[k];
    bd += b[k]*dray[k];
    cd += c[k]*dray[k];
  }
/*   if(ab < 0.0 || bc < 0.0) {		/# suspect: might be > 180 angle  #/
 *     CROSS3P(normal??); 
 *   }
 */  /* d between a & b: */
  dot1=ad*bd-ab;			/* face angle adb */
  dot2=bc*bd-cd;			/* face angle dbc */
  
  /* d between b & c: */
  dot3=ab*bd-ad;			/* face angle abd */
  dot4=bd*cd-bc;			/* face angle bdc */
  
  if(dot3+dot4 > dot1+dot2 )		/*  d should be between b & c,  */
    maxp= & (*maxp)->next;		/* so change insertion point */
  
  return maxp;
} /* con_place */

static int sort_connections(box_p box) {
  /*
   * does    : sorts all connections out of all bndries of shell
   *
   * gets    : shell
   *
   * affects : ordering of shell->patches->bndries->cons[i]; may or may
   *   not delete connections deemed invalid
   *
   * returns : nr of connections deleted
   *
   * warns   : 
   *
   * comment : 
   */
  int i, ii, length, deleted;
  patch_p patch;
  bndry_p bndry;
  list_p rest, l, next, *w;
  shell_p shell;

  deleted=0;
  for (ii=0; ii<box->natoms; ii++) { 
    shell=box->shells + ii;
    if ( remove_widows(shell, &deleted) )
      continue;
    for(patch=shell->patches; patch; patch=patch->next) {
      for(bndry=patch->bndries; bndry; bndry=bndry->next) {
	assert( !(bndry->flags & ISWIDOW));
	length=bndry->length;
	for (i=0; i<length; i++) {	/* walk complete path */
	  if (bndry->cons[i]) {
	    rest=bndry->cons[i]->next;
	    bndry->cons[i]->next=NULL;	/* can leave one on it: is sorted */
	    for (l=rest; l; l=next) {	/* add all of them one by one */
	      next=l->next;
	      w=con_place(l->ptr, bndry, i);
	      if ( w==NULL ) {
		warn("rejecting edges during sorting\n"); /* other end: */
		deleted++;
		remove_con_halfknown(l, bndry);
		continue;
	      }	/* w==NULL */
	      l->next= *w;
	      *w=l;
	    } /* for l */
	  } /* if != NULL */
	} /* for i < length */
      } /* for bndry */
    } /* for patch */
    remove_widows(shell, &deleted);
  } /* for ii < box->natoms */
  return (deleted > 0);
} /* sort_connections */

static void accept_cons(con_p prelim_cons, bndry_p bndry, bndry_p bn) { 
  /* inserts all prelim_cons into bndry and bn such, that any point of */
  /* bn(dry) has at most one con to another bndry  */
  int i, length, offset, o, p;
  con_p c;
  const con_t * con;
  list_p *cons, *ngbcons, l;
  float dist;

  length=bndry->length;
  cons=bndry->cons;
  ngbcons=bn->cons;

  for (i=0; i<length; i++) {
    con= &prelim_cons[i];
    if(con->bndries[0]) {		/* if not NULL */
      dist=con->direction[3];
      /* first check the neighbours to see if it already has a con to it */
      offset=con->flags[0].offset;	/* offset in bndry */
      o=con->flags[1].offset;		/* offset in bn */

      l=ngbcons[o];
      if(l &&				/* not empty; see which is better */
	 ((c=l->ptr)->bndries[0])==bndry ) { /* previous also to bndry */
	if(dist < c->direction[3]) {	/* new one better; throw away old */
	  p = c->flags[0].offset;	/* previous offset in bndry */
	  l=cons[p]; cons[p]=l->next;	/* excise ... */ 
	  Lprepend(l, cons[ offset ]);	/* ... and rewire */
	  *c = *con;			/* overwrite ! */
	} /*  else: old one better; forget about it */
      } else { 				/* bndry is not yet a neighbour */
	c=alloc_con(con);		/* so make it one  */
	cons[ offset ]=Lcons(c, cons[ offset ] ); /* new link on bndry */
	ngbcons[ o ]=Lcons(c, ngbcons[ o ] ); /* new link on bn */
	bndry->ncons++;
	bn->ncons++;
      }	/* else */
    } /* if (con->bndries[0]) */
  } /* for i < length */
  return;
} /* accept_cons */

static int connect_2_bndries(bndry_p bndry, bndry_p bn, 
			     geom_p geometry, int mode) {
  /*
   * does    : connects two bndries, and gets rid of crossing edges
   *
   * gets    : two bndries to be connected, geometry, and shells involved
   *
   * affects : bn{dry}->cons[i], which are arrays of linked lists of con's
   *
   * returns : nr of connections made
   *
   * warns   : 
   *
   * comment : 
   */
  int i, j, k, p, q, ncons, nc, length, len, minj;
  SHORT *path, *qath;
  float *coorda, *coordb, *Xp, *xq, *xa, *xb;
  shell_p shell, sh, *np, *nq;
  double xp[3], d2, d, delta, ra_rb, ra, rb, saxis[3], dir[3], mindir[3]
    , mind2;
  list_p conlist;
  con_t prelim_cons[ MAX_PATHLENGTH ];

  if(mode)
    bzero((char*)prelim_cons, sizeof(prelim_cons));

  ncons=0;
  shell=geometry->shella;
  sh=geometry->shellb;
  coorda=geometry->coorda;		/* points' coordinates */
  coordb=geometry->coordb;
  xa=geometry->xa;			/* atom's coordinates */
  xb=geometry->xb;
  ra=geometry->xa[3];
  rb=geometry->xb[3];
  ra_rb=ra/rb;
  
  d2=0.0;
  for (k=0; k<3; k++) {
    d=xa[k] - xb[k];
    saxis[k]= d/rb;			/* scaled inter-atom axis */
    d2 += SQR(d);			/* (unscaled) squared distance */
  }

  delta=calc_delta(d2, ra, rb, shell->sphere, sh->sphere);

  length=bndry->length;
  len=bn->length;
  path=bndry->path;
  qath=bn->path;

  for (i=0; i<length; i++) {
    bzero((char*)(prelim_cons+i), sizeof(con_t));
    for(np=bndry->pt.ngbs[i]; *np; np++) {
      if(*np==sh) {			/* sh is ngb of this point: go ahead */
	p=path[i];			/* (will happen only once) */
	conlist=NULL;			/* list of con's out of p  */
	nc=0;				/* their number */

	Xp=coorda +3*p;			/* point's unit sphere coordinates */
	for (k=0; k<3; k++) 
	  xp[k] = ra_rb*Xp[k] + saxis[k]; /* 'scaled point' */
	
	/* find i's closest neighbour: */
	mind2=HUGE_VAL;
	minj= -1;			/* sentinel */
	for (j=0; j<len; j++) {		/* walk this bndry */
	  for(nq=bn->pt.ngbs[j]; *nq; nq++)
	    if (*nq==shell)  {		/* q does have shell a ngb; proceed */
	      q=qath[j];		/* candidate point on other bndry */
	      xq= coordb+3*q;
	      d2=0.0;
	      for (k=0; k<3; k++) {
		d = dir[k] = xq[k]-xp[k];
		d2 += SQR(d);
	      }

	      if ( d2 < delta		/* close enough to be a connection? */
		  && d2 < mind2) {	/* and closer than prev. */
		mind2=d2;
		minj=j;
		memcpy(mindir, dir, sizeof(dir));
	      }	/* if d2 < delta */
	    } /* if *nq==shell */
	} /* for j < bn->length */
	if(minj>=0 )
	  ncons += keep_con(prelim_cons, bndry, i, bn, minj, 
			    mindir, (rb*sqrt(mind2)), geometry);
      } /* if *np==sh */
    } /* for(np=bndry->pt.ngbs[i]; *np; np++) */
  } /* for i<bndry->length */

  geometry->saxis=saxis;		/* needed for gap_con() */
  
  if(ncons)				/* found connections */
    accept_cons(prelim_cons, bndry, bn); /* & put in bndry & bn */


  return ncons;
} /* connect_2_bndries */

/* #define BLOCK 8				/# by which to extend ALLOCATION #/
 * #define BITS (BLOCK-1)			/# (BLOCK should be a power of 2) #/
 */

static  int unbury(bndry_p bndry) { 
  int i,j, q,b,n, indir, new, k, len;
  shell_p shell;
  point_p Points, pp;
  SHORT *pointflags;

  n=0;
  shell=bndry->patch->atom->pointer;
  pointflags=shell->pointflags;

  Points=shell->sphere->Points;
  pp=Points+bndry->path[0];
  q=bndry->path[1];

  for (i=0; i<6; i++)			/* find starting direction */
    if (pp->ngbs[i]==q)
      break;
  indir=pp->whoamis[i];			/* incoming dir at q */
  pp=Points+q;

  len=bndry->length;
  for (i=1; i<len+1; i++) {
    q=bndry->path[ (i+1)%len ];
    for(j=1; 1; j++ ) {
      assert(j<=3);			/* path should be convex ! */
      k=MOD(indir+j, pp->nngbs);
      b=pp->ngbs[k];
      if (b == q)
	break;
      new= ((pointflags[b] & (BURIED | UNWANTED)) > 0);
      n+=new;
      pointflags[b] &= ~(BURIED | UNWANTED);
    } /* for j */
    indir=pp->whoamis[ k ];		/*  */
    assert(q==pp->ngbs[ k ]);
    pp=Points+b;
    assert(pp->nr==q);
  } /* for i <  length */
  assert(n);				/* nr of quasi-buried points */
  /*   assert(n <= 4);*/               /* (4 was seen once in 10,000 cases) */
  assert(n <= 5);
  return n;
} /* unbury */

static int rescue_bndry(bndry_p bndry) {
  /* stupid function that triangulates the interior of single bndry */
  shell_p shell;

  if (bndry ->length < 5)
    return 0;

  shell=bndry->patch->atom->pointer;
  if(bndry->patch->bndries->next == NULL &&  /* not a multiple bndry */
     bndry->patch->npoints < 2*bndry->length) { /* and not nearly isolated */
    bndry->flags |= ISWIDOW;
    shell->flags |= REMOVE_WIDOWS;
    return 0;
  }
  return unbury(bndry);
} /* rescue_bndry */


static int prelim_connections(box_p box, int mode) {
  /*
   * does    : for shell, connects all its bndries to neighbouring sh's >
   *   shell. This is the place where the preliminary interatom edges are
   *   found. Later on they're cleaned up/supplemented to complete the
   *   triangulation. Only bndry->cons[0..lenght] are generated
   *
   * gets    : a shell
   *
   * affects : both shell->patches[*]->bndries[*]-> ... and those of ngbs
   *
   * returns : nr of bridges found for shell
   *
   * warns   : 
   *
   * comment : mainly a bookkeeping function; real work done by bndry_bridge()
   */
  int i, j, ii, deleted, n;
  shell_p shell, sh;
  patch_p patch, pa;
  bndry_p bndry,bn;
  bndry_p bndry_ngbs[ MAX_NNGBS ];
  geom_t geometryt, *geometry = &geometryt;

  deleted=0;
  for (ii=0; ii<box->natoms; ii++) { 
    shell=box->shells+ii;

    if ( remove_widows(shell, &deleted) )
      continue;
    geometry->coorda=shell->sphere->coordinates;
    geometry->xa=shell->atom->xyz;
    geometry->shella=shell;

    for(patch=shell->patches; patch; patch=patch->next) {
      for(bndry=patch->bndries; bndry; bndry=bndry->next) {
	bzero( (char*)bndry_ngbs, MAX_NNGBS*sizeof(bndry_p));
	for (i=0; i<bndry->nngbs; i++) {
	  sh=bndry->ngbs.shell[i];	/* union of all pt.ngbs[*] */
	  if (sh > shell) {		/* important ! avoid doing it twice! */
	    for(pa=sh->patches; pa; pa=pa->next)
	      for(bn=pa->bndries; bn; bn=bn->next)
		for (j=0; j<bn->nngbs; j++) /* search for mutual neighours: */
		  if(bn->ngbs.shell[j]==shell) { /* (need not be case) */
		    geometry->coordb=sh->sphere->coordinates;
		    geometry->xb=sh->atom->xyz;
		    geometry->shellb=sh;
		    connect_2_bndries(bndry, bn, geometry, mode);
		  } /* for i< bndry->nngbs */
	  } /* if sh > shell */
	} /* for i < bndry->nngbs */
	if (bndry->ncons == 0)	{	/* no connections found for bndry: */
	  if ( (n=rescue_bndry(bndry)) )
	      warn("\
un-buried %d points in interior of UNFAIR boundary of length %d %s\n",
		   n, bndry->length, atom_spec(box, ii));
	  bndry->flags |= ISWIDOW;
	  shell->flags |= REMOVE_WIDOWS;
	  bndry->length=0;		/* otherwise, its points will be */
					/* set to UNWANTED */
	  shell->flags &= ~MBNDRIES;	/* should do find_patches(),in fact */
	} /* if no connections */
      } /* for bndry */
    } /* for patch */
    remove_widows(shell, &deleted);
  } /* for ii<box->natoms */
  return 0;
} /* prelim_connections */

#define MAXPATCHNEIBORS 32		/* max nr of neighbouring patches */
static int connect_patches(box_p box) {
  int i, j, nngbs, end, err, ii, deleted;
  con_p con;
  patch_p patch, pa, ngbs[ MAXPATCHNEIBORS ];
  bndry_p bndry, bn;
  list_p l;
  shell_p shell;

  deleted=0;
  for (ii=0; ii<box->natoms; ii++) { 
    shell=box->shells + ii;
    if ( remove_widows(shell, &deleted ) ) 
      continue;
    err=0;
    for(patch=shell->patches; patch; patch=patch->next) {
      nngbs=0;
      for(bndry=patch->bndries; bndry; bndry=bndry->next) { 
	assert( ! ( bndry->flags & ISWIDOW ) );	/* most important */
	assert(bndry->ncons);
	for (i=0; i<bndry->length; i++) 
	  for (l=bndry->cons[i]; l; l=l->next) {
	    con=l->ptr;
	    end=(con->bndries[0] == bndry);
	    bn=con->bndries[end];
	    pa=bn->patch;		/* candidate for being new ngbr */
	    for (j= 0; j < nngbs ; j++)
	      if( pa == ngbs[j])	/* already present */
		goto NEXT_CON;
	    ngbs[ nngbs++ ]=pa;
	      assert( nngbs <= MAXPATCHNEIBORS );
	  NEXT_CON:
	    ;
	  } /* for l */
      }	/* for bndry */
      assert(nngbs);			/* otherwise a WIDOW */
      if(nngbs <= 2 && patch->npoints > 1) /* somewhat unusual */
	warn("patch has only %d neighbouring patches %s\n", 
	     nngbs, atom_spec(box, ii));
      patch->nngbs=nngbs;
      if(nngbs)
	patch->ngbs=memsav(ngbs, nngbs*sizeof(patch_p));
    } /* for patch */
  } /* for ii < natoms */
  return deleted;
} /* connect_patches */

int connections(box_p box, int mode) {
  int i;
  /* make connections between all bndries; if mode==1, do it in */
  /* patch-mode: don't reduce delta, and 1 connection per bndry suffices */

  i=prelim_connections(box, mode);

  if(mode) {				/* patch mode */
    connect_patches(box);
    return i;
  }

  free_all_neibordata(box);		/* has to be done after prelim ! */
  clean_conns(box);			/* get rid of 'longest' ones */
  sort_connections(box);
  stitch(box);				/* rest of triangulation */
  connect_patches(box);
  return NConsTotal;			/* exporting file global variable */
} /* connections */

void free_connections(box_p box) {	/* frees just cons, not refs ! */
  int i, j, natoms;			/* (that is done by free_bndries) */
  shell_p shell;
  patch_p p;
  bndry_p bndry;
  con_p con;
  list_p l, lnext;

  natoms=box->natoms;

  for (i=0,shell=box->shells; i<natoms; i++, shell++) {
    for(p=shell->patches; p; p=p->next) {
      for(bndry=p->bndries; bndry; bndry=bndry->next) {
	if (bndry->cons) {
	  for (j=0; j<bndry->length; j++) { /* walk bndry, freeing all lists */
	    for (l=bndry->cons[j]; l; l=lnext) {
	      lnext=l->next;
	      con=l->ptr;
	      if (con->bndries[0])	/* not NULL: first time we see it */
		con->bndries[0]=NULL;	/* mark for next time */
	      else
		free_con(con);		/* next time: free it */
	      FREE(l);
	    } /* for l */
	  } /* for j */
	} /* if bndry->cons,  */
      }/* for bndries */
    } /* for patches */	    
  } /* for i<natoms */
  return;
} /* free_connections */
