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

/* routines for establishing topology of patches on one atom */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "extra-math.h"

#include "utils.h"
#include "alloc.h" 
#include "list.h"
#include "pdb.h"
#include "vecmat.h" 

#include "quilt.h"

struct {
  box_p box;
  int firsttime;
} fileglobal;

typedef struct {
  double axes[ 3*MAX_NNGBS ];
  double dotps[ MAX_NNGBS ];
} geom_t, *geom_p;

#define MAX_PT_NNGBS MAX_NNGBS

static bndry_p save_bndry(shell_p shell, const bndry_t *bndry) { 
  /*
   * does    : saves (allocates & copies) a just found bndry, and returns it
   *
   */
  int length;
  bndry_p bn;

  bn = (bndry_p)memsav(bndry, sizeof(bndry_t));
  length=bn->length;
  assert(length>0);
  bn->path = memsav(bndry->path, length*sizeof(SHORT));
  bn->cons=(list_p*)CALLOC(1, length*sizeof(list_p*));

  if (fileglobal.firsttime) {
    if(bn->nngbs > 0) { 
/*       if (bn->pt.ngbs)
 * 	FREE(bn->pt.ngbs);
 *       if (bn->ngbs.shell)
 * 	FREE(bn->ngbs.shell);
 */
      bn->pt.ngbs = memsav( bndry->pt.ngbs, length*sizeof(shell_p*) );
      bn->ngbs.shell = memsav(bndry->ngbs.shell, bn->nngbs * sizeof(shell_p));
    } else { /* a widow */
      assert(bndry->nngbs == 0);
      bn->ngbs.bndry=NULL;		/* otherwise it might get FREEd */
      bn->pt.ngbs=NULL;			/* likewise */
      bn->flags |= ISWIDOW;
      shell->flags |= REMOVE_WIDOWS;
    }
  } else {
    /* WIDOW_HOOD cannot be decided until after transfer_connections */
    bn->nngbs=0;
    bn->pt.ngbs=NULL;
    bn->ngbs.shell=NULL;
  }
  return bn;
} /* save_bndry */

static void ngb_geometry(shell_p shell, geom_p geom) {
  /*
   * does    : rebuilds axes to, and dot products with neighbours
   *
   * gets    : shell to apply it on; places to leave the axes and dotps
   *
   * affects : axes[*] and dotps[*]
   *
   * returns : nothing
   *
   * warns   : 
   *
   * comment : 
   */
  int j, natngbs;
  shell_p * atngbs;
  float *xa, *xb;
  double  *axes, *dotps, *axis, d2, ra, rb, sdotp;

  axes=geom->axes;			/* an lvalue */
  dotps=geom->dotps;			/* an lvalue */

  atngbs=shell->ngbs;
  natngbs=shell->nngbs;
  for (j=0; j<natngbs; j++) {
    axis=axes+3*j;			/* lvalue */
    xa=shell->atom->xyz;
    xb=atngbs[j]->atom->xyz;
    d2=DIST3_2(xa,xb);
    ra=xa[3]; rb=xb[3];
    sdotp = ((d2 - SQR(rb))/ra + ra)*0.5;
    VEC3MIN(axis, xb, xa);		/* assign */
    dotps[j] = sdotp;
  } /* for j */
  return;
} /* ngb_geometry */

static int find_buriors(shell_p shell, bndry_p bndry, 
			point_p pp, int indir, int outdir, 
			const geom_t *geom, shell_p *ptngbs) {
  /* finds which shells bury neigbour points of pp; writes into ptngbs */
  int i, j, k, nngbs, nptngbs;
  int done[ MAX_PT_NNGBS ];
  const double  *axes, *dotps, *axis;
  double dot;
  shell_p sh;

  axes=geom->axes; 
  dotps=geom->dotps;

  nptngbs=0;
  assert(shell->nngbs > 0 && shell->nngbs < MAX_PT_NNGBS);
  bzero((char*)done, shell->nngbs*sizeof(int));
  nngbs=pp->nngbs;
  for (i=(indir+1)%nngbs; i != outdir; i = (i+1)%nngbs) { 
    assert(shell->pointflags[ pp->ngbs[i] ] & (BURIED | UNWANTED ));
    for (j=0; j<shell->nngbs; j++) {
      axis=axes + 3*j;			/* axis towards ngb j */
      dot=DOT3P(pp->ngbs_xyz[i], axis);
      if (dot >= dotps[j] 	 	/* this pt lies buried within atom j */
	  && !done[j]) {		/* and was not yet found */
	sh=shell->ngbs[j];
	ptngbs[nptngbs++]=sh;		/* add shell to list of pt ngbs */
	done[j]=1;
	/* keep bndry neighbours; avoid doubles: */
	for (k=bndry->nngbs-1; k>=0; k--) 
	  if(bndry->ngbs.shell[k] == sh)
	    break;
	if(k<0)
	  bndry->ngbs.shell[ bndry->nngbs++ ] = sh;
      }
    } /* for j */
  } /* for i */
  return nptngbs;
} /* find_buriors */

static void save_ptngbs(shell_p shell, bndry_p bndry, point_p pp, 
			int indir, int outdir, const geom_t *geom) {
  int size, nptngbs;
  shell_p * pn;
  shell_p ptngbs[ MAX_PT_NNGBS ];	/* find_buriors will write into this*/


  if (fileglobal.firsttime==0)		/* this is messy ... */
    return;

  nptngbs=find_buriors(shell, bndry, pp, indir, outdir, geom, ptngbs);
  /* (also updates the bndry->ngbs.shell thing) */
  assert(nptngbs <= MAX_PT_NNGBS);

  size=(nptngbs+1)*sizeof(shell_p);	/* 1 extra for terminating NULL */
  pn=(shell_p*)MALLOC(size);
  memcpy(pn, ptngbs, size);
  pn[ nptngbs ]=NULL;			/* NULL terminate it */

  bndry->pt.ngbs[ bndry->length ]=pn;

  return;
} /* save_ptngbs */

static point_p find_start_pt(shell_p shell, int start, int*dirp) {
  int i, w, dir, maxlen, nngbs;
  point_p pp, qp;
  SHORT  *pointflags;


  if ( start >= 0 ) {			/* MBNDRY, doing it a 2nd time */
    assert(*dirp>=0);
    return  &shell->sphere->Points[start]; /* dirp should be OK */
  }

  pointflags=shell->pointflags;
  
  /* Have to find a point + direction on the boundary */
  /* First find an exposed point that's not yet part of a BNDRY or PATCH */
  for (i=0; i<shell->sphere->npoints; i++)
    if (  (pointflags[i] & (BURIED | UNWANTED | BNDRY | PATCH )) == 0)
      break;

  assert(i<shell->sphere->npoints);	/* should have been found */
  qp = &shell->sphere->Points[i];	/* qp now an acceptable point */
  
  /* From here, walk towards a a BURIED | UNWANTED point */
  /* walk in some random direction, but keep it approximately same: */
  dir= -1;
  i=w=0;				
  maxlen=shell->sphere->maxringlen;
  do {
    pp=qp;				/* old pt */
    nngbs = pp->nngbs;
    dir=MOD( w + nngbs/2, nngbs);	/* at other side*/
    w = pp->whoamis[ dir ];
    if (i++ > maxlen) {			/* running around in circle; perturb */
      i=0;
      w = rand() % nngbs;		/* infrequent, so can afford call */
    }
    qp=pp->Ngbs[ dir ];
  } while_not ( (pointflags[qp->nr] & (BURIED | UNWANTED) ));
  assert(dir>=0);
  assert(!(pointflags[pp->nr] & (BURIED | UNWANTED) ));

  *dirp=dir;
  return pp;
} /* find_start_pt */

static int next_bndry_pt(point_p pp, SHORT *pointflags, int indir) {
  int i, nngbs, outdir, nr;
  SHORT flags;

  nngbs=pp->nngbs;
  for (i=1; i<nngbs; i++) {
    outdir= MOD(indir+i, nngbs);
    nr=pp->ngbs[outdir];
    flags = pointflags[nr];
    if_not(flags & (BURIED | UNWANTED)) /* done */
      return outdir;			/* that's all */
  } /* for i */
  return -1;				/* not found: ISOLATED point */
} /* next_bndry_pt */

static bndry_p isolated_pt_bndry(shell_p shell, bndry_p bndry,  patch_p patch,
				 point_p pp, const geom_t *geom) {
  shell->pointflags[pp->nr] |= (ISOLATED | PATCH | BNDRY); /* all in one */
  patch->points[0]=bndry->path[0]=pp->nr;
  bndry->ncons=1;
  patch->npoints=1;
  save_ptngbs(shell, bndry, pp, 0, pp->nngbs-1, geom);
  bndry->length=1;
  return save_bndry(shell, bndry);
} /* isolated_pt_bndry */

static int add_bndry_pt(bndry_p bndry, point_p pp, 
			SHORT *flagp, patch_p patch) {
  /* sets flags, and updates counts. Returns wether shared */
  int nr, s;

  nr=pp->nr;
  bndry->path[ bndry->length ]=nr;

  s=0;
  if (  (*flagp & (BNDRY | PATCH) ) == (BNDRY | PATCH) ) { 
    /* seen point before, current invocation */
      *flagp |= TWOWAY;			
  }
  else { 
    s = ((*flagp & (BNDRY | PATCH) )== BNDRY);
    /* s is normally 0, but 1 if the point is SHARED. This is because
       in the current invocation, all boundary points are marked both
       PATCH and BNDRY; if multiple boundaries are marked, however, the
       PATCH bit is unset.
     */
    *flagp |= ( BNDRY | PATCH | (s*SHARED) ) ;
    bndry->flags |= (s*SHARED);         /* phew! */
    bndry->ncons++;			/* npoints here !!!!! */
    patch->points[ patch->npoints++ ]=nr;
  }
  return s;
} /* add_bndry_pt */

static bndry_p find_boundary(shell_p shell, int start, patch_p patch,
			     const geom_t * geom, int * directions,
			     int *sharedp) {
  /* finds a boundary , and returns allocated version of that. */
  /* *sharedp is incremented by the nr of SHARED points found */
  int indir, outdir, s;
  point_p pp;
  void * start_edge, *e;
  SHORT *pointflags;
  bndry_t bndryt;
  SHORT bndrypath[ MAX_PATHLENGTH ];	/* temporary storage */
  shell_p *ptngbs[ MAX_PATHLENGTH ];
  shell_p bndryngbs[ MAX_NNGBS ];

  bzero((char*)&bndryt, sizeof(bndry_t));
  bndryt.path= &bndrypath[0];
  bndryt.pt.ngbs= &ptngbs[0];
  bndryt.ngbs.shell= &bndryngbs[0];

  pointflags=shell->pointflags;

  pp=find_start_pt(shell, start, &directions[0]);
  outdir=next_bndry_pt(pp, pointflags, directions[0] );
  if (outdir < 0)		    /* vvv provides storage space */
    return isolated_pt_bndry(shell, &bndryt, patch, pp, geom);

  e=NULL; 
  indir=pp->whoamis[outdir];		/* indir for next point */
  start_edge = &(pp->Ngbs[ outdir ]);
  pp=pp->Ngbs[ outdir ];		/* real starting pt */
					/* find_start_pt found start edge ! */
  do {
    s=add_bndry_pt(&bndryt, pp, pointflags+pp->nr, patch);
    *sharedp += s;
    shell->flags |= SHARED*s;
    outdir=next_bndry_pt(pp, pointflags, indir); /* new direction */
    assert(indir != outdir);
    if ( outdir < 0 ) {
      pointflags[ pp->nr ] |= UTURN;
      outdir=indir;
    }
    directions[ bndryt.length ]=outdir;
    save_ptngbs(shell, &bndryt, pp, indir, outdir, geom);
    bndryt.length++;
    indir=pp->whoamis[outdir];		/* indir for next point */
    e= &(pp->Ngbs[outdir]);
    pp = pp->Ngbs[outdir];		/* new pt */

  } while_not ( e == start_edge );
  return save_bndry(shell, &bndryt);
} /* find_boundary */

static point_p find_interior_pt(shell_p shell, bndry_p bndry, 
				const int *directions,  
				int *indirp, int *outdirp) {
  int i, j, outdir;
  point_p pp, Points;
  SHORT *pointflags;

  Points=shell->sphere->Points;
  pointflags=shell->pointflags;

  for (i=0; i<bndry->length; i++) {	/* find interior point */
    pp= &Points[ bndry->path[i] ];
    for (j=1; j<4; j++) {		/* wild guess */
      outdir=MOD(directions[i] + j, pp->nngbs);
      if ( ( pointflags[ pp->ngbs[ outdir ] ] &
	    (BURIED | UNWANTED | PATCH))==0 ) {
	*outdirp=outdir;
	*indirp=pp->whoamis[ outdir ];
	return pp->Ngbs[ outdir ];
      }	/* if pointflags ... */
    } /* for j */
  } /* for i */
  return NULL;
} /* find_interior_pt */

static int check_history(int m, int *indir_history, int *outdir_history, 
		      SHORT *pointflags, point_p *point, 
		      int *indirp, int*outdirp) {
  int i, nngbs, in, in2, out, nr, flags, outdir, indir;
  point_p pp;
  /* here we end up if we can't find an unmarked point: backtrack! */
  pp= *point; 

  while( --m ) {
    /* from here, role of h.in and h.out is reversed ! */
    in=indir_history[m];
    out=outdir_history[m];
    pp=pp->Ngbs[ in ] ;			/* previous point */
    in2=indir_history[m-1];		/* next point */
    nngbs=pp->nngbs;
      
    for (i=1; i<nngbs; i++) {		/* revolve around point in the */
      outdir=MOD(out+i,nngbs);		/* possibly unmarked part ('left'), */
      if (outdir==in2)			/* back on the traced path without */
	break;				/* having found unmarked point  */
      nr=pp->ngbs[outdir];
      flags=pointflags[nr];
      if ( ! (flags & (PATCH | BURIED | UNWANTED))) { 
	/* point unmarked: accept */
	indir=pp->whoamis[outdir];
	pp=pp->Ngbs[outdir];
        *indirp = indir;
        *outdirp = outdir;
	*point=pp;
	return m;
      }
    } /* for i */
  } /* while(--m) */
  return 0;
} /* check_history */

static int mark_interior(shell_p shell, int wanted, 
			 patch_p patch, int *directions) {
  int i, outdir, indir, nngbs, m, nr;
  SHORT flags, oldflags, *patchpoints, *pointflags;
  int outdir_history[ MAX_TMPLT_PNTS ], indir_history[MAX_TMPLT_PNTS ];
  point_p pp;

  if (patch->npoints == wanted || patch->bndries->ncons <= 4 ) /* npoints! */
    return 0;
  pp=find_interior_pt(shell, patch->bndries, directions, &indir, &outdir);
  if (pp == NULL)
    return 0;

  /* init: */
  patchpoints=patch->points;
  pointflags=shell->pointflags;
  assert(patch->npoints > 0);
  assert(patch->npoints < wanted);

  m=0;
  nr=pp->nr;
  flags = pointflags[nr];

 NEW_POINT:				/* think of this as a big while(1) */
  assert((flags & (BURIED | UNWANTED | PATCH ))==0); 
					/* can be BNDRY ! (if MBNDRIES) */
  patchpoints[ patch->npoints++]=nr;	/* includes BNDRY points !!!  */

  outdir_history[m]=outdir;
  indir_history[m]=indir;
  m++;

  pointflags[nr] |= PATCH;

  nngbs=pp->nngbs;
  oldflags = flags;
  for (i=1; i<nngbs; i++) { 
    outdir = MOD(indir + i, nngbs);
    nr= pp->ngbs[outdir];
    flags = pointflags[nr];
    if ( flags & (BURIED | UNWANTED) )  {
      if_not (oldflags&BNDRY) {		/* new, unknown boundry */
	directions[0] = outdir;		/* find_boundary fills rest of it */
	return 1;
      } /* else */
      continue;
    } /* if BURIED | UNWANTED */
    if ( ! (flags & (PATCH) )) { 
      indir = pp->whoamis[outdir];
      pp=pp->Ngbs[outdir];
      goto NEW_POINT;
    }
  } /* for i < nngbs */
  
  assert(patch->npoints <= wanted);		/* cannot find more ! */
  if (patch->npoints == wanted)
    return 0;				/* ready */
  
  m=check_history(m, indir_history, outdir_history, pointflags, 
		  &pp, &indir, &outdir); /* their new values */
  if(m) {				/* if m==0, no missed branch was */    
    nr=pp->nr;
    flags=pointflags[nr];
    goto NEW_POINT;			/* found; otherwise, it's new m */
  }

  /* may have missed interior point; check rest of the boundary:  */
  pp=find_interior_pt(shell, patch->bndries, directions, &indir, &outdir);
  if (pp) {
    m=0;
    nr=pp->nr;
    flags=pointflags[nr];
    goto NEW_POINT;
  }
  /* didn't find anything during backtrack: finished */
  return 0;				/* all OK */
} /* mark_interior */

void free_neibordata(shell_p shell) {
  /* frees data on point and boundary neibors */
  int i;
  bndry_p bn;
  patch_p pa;

  for (pa=shell->patches; pa; pa=pa->next)
    for(bn=pa->bndries; bn; bn=bn->next) {
      if (bn->pt.ngbs) {
	for (i=0; i<bn->length; i++) 
	  if(bn->pt.ngbs[i])
	    FREE(bn->pt.ngbs[i]);
	FREE(bn->pt.ngbs);		/* becomes normals pointer */
	bn->pt.ngbs=NULL;
      }
      if (bn->ngbs.shell) {
	FREE(bn->ngbs.shell);
	bn->ngbs.shell=NULL;		/* becomes coordinate pointer */
      }
    } /* for pa, bn  */
  return;
} /* free_neibordata */

void free_all_neibordata(box_p box) {
  int ii;
  shell_p shell;

  for (ii=0; ii<box->natoms; ii++) { 
    shell=box->shells + ii;
    if (INVALID_SHELL(shell))
      continue;
    free_neibordata(shell);
  }
  return;
} /* free_all_neibordata */

int free_bndry(bndry_p bn) { 
  int i;

  assert(bn);
  if ( bn->pt.ngbs ) {			/* may have been free()ed by repair */
    if_not (bn->flags & COORDINATES)	/* NOT IF it's coordinates ! */
      for (i=0; i<bn->length; i++) {
	if( bn->pt.ngbs[i] ) { 
	  FREE(bn->pt.ngbs[i]);		/* TWOWAY point */
	  bn->pt.ngbs[i]=NULL;		/* so we don't do it next time */
	} else
	  assert( ((shell_p)bn->patch->atom->pointer)
		 ->pointflags[bn->path[i]] & TWOWAY );
      } /* for i< length */
   FREE(bn->pt.ngbs);			/* later reallocated as ngbs.bndry */
    bn->pt.ngbs=NULL;
  } /* if pt.ngbs */

  if(bn->ngbs.bndry) { 
    FREE(bn->ngbs.bndry);		/* allocated first as ngbs.shell, */
    bn->ngbs.bndry=NULL;
  }
  if (bn->cons)
    FREE(bn->cons);
  assert(bn->path);
  FREE(bn->path);			/* in con.c */
  FREE(bn);
  return 0;
} /* free_bndry */

static void free_bndries(bndry_p bndries) {
  /* frees bndries. Always after a 'stitch': pt.ngbs are normals !   */
  bndry_p bn, next;

  for(bn=bndries; bn; bn=next) {
    next=bn->next;
    free_bndry(bn);
  } /* for bndries */
  return;
} /* free_bndries */

int free_patch(patch_p patch) { 
  /* frees a complete patch; returns nr of points it had */
  int n;

  n=0;
  free_bndries(patch->bndries);
  if (patch->points)			/* need not be case (?) */
    FREE(patch->points);
  if(patch->ngbs)			/* need not be case */
    FREE(patch->ngbs);			/* doesn't belong here, strictly */
  if (patch->points)
    FREE(patch->points);
  FREE(patch);
  return n;
} /* free_patch */

int free_patches(patch_p patches) { 
  int n;
  patch_p pa, next;

  n=0;
  if(patches==NULL)
    return 0;
  for(pa=patches; pa; pa=next) {
    next=pa->next;
    n += free_patch(pa);
  } /* for patches */
  return n;
} /* free_patches */

static void free_conlists(shell_p shell) { 
  int i;
  list_p ***conlists;

  conlists=shell->next.conlists;
  if (conlists==NULL)
    return;
  for (i=0; i<shell->sphere->npoints; i++) 
    if(conlists[i])			/* if not BURIED! */
      FREE(conlists[i]);

  FREE(conlists);
  return;
} /* free_conlists */

void free_topology(box_p box) {
  int i;
   
  for (i=0; i<box->natoms; i++) {
    free_conlists(box->shells + i);
    free_patches(box->shells[i].patches);
  }
  return;
} /* free_topology */

static double patch_area(shell_p shell, const patch_t * patch) {
  int i, npoints;
  SHORT * points; 
  point_p Points;
  double area;

  npoints=patch->npoints;
  points=patch->points;
  Points=shell->sphere->Points;
  area=0.0;
  for (i=0; i<npoints;  i++) 
    area += Points[ points[i] ].area;
  return area;
} /* patch_area */

static patch_p save_patch(shell_p shell, const patch_t * patch,
			  int save_internal) {
  /*
   * does    : saves (allocates and copies) newly found patch, and returns it
   *
   */
  patch_p pa;
  bndry_p bn;

  pa = memsav(patch, sizeof(patch_t));
  if(save_internal)
    pa->points=memsav(pa->points, pa->npoints*sizeof(SHORT));
  else
    pa->points=NULL;

  for (bn=patch->bndries; bn; bn=bn->next)
    bn->patch=pa;			/* patch they belong to */

  pa->area = patch_area(shell, patch);

  pa->atom=shell->atom;			/* atom patch belongs to */
  return pa;
} /* save_patch */

static int set_conlists(shell_p shell) {
  /* allocates and sets lists of pointer to roots of conlists  */
  int i, l, npoints, p, lengths[MAX_TMPLT_PNTS];
  patch_p pa;
  bndry_p bn;
  list_p *** conlists;			/* YO */

  npoints=shell->sphere->npoints;

  bzero((char*)lengths, npoints * sizeof(int));
  
  if(shell->next.conlists==NULL)
    shell->next.conlists=
      (list_p***)CALLOC(npoints, sizeof(list_p**));
  else					/* doing this shell 2nd time */
    ; /* assert(shell->flags & UNWANTED); */

  conlists=shell->next.conlists;

  for(pa=shell->patches; pa; pa=pa->next) {
    for(bn=pa->bndries; bn; bn=bn->next) {
      assert(bn->cons);
      for (i=0; i<bn->length; i++) {
	p=bn->path[i];
	if (conlists[p]==NULL)
	  conlists[p]=(list_p**)MALLOC(4*sizeof(list_p*) );
	l=lengths[p]++;			/* current length */
	assert( l < 3);
	conlists[p][l] = bn->cons+i;
	conlists[p][l+1]=NULL;
      }	/* for i */
    } /* for bn */
  } /* for pa */
  return 0;
} /* set_conlists */

static int shared2unwanted(shell_p shell)  {
  /* makes all SHARED points UNWANTED */
  int i, n;
  SHORT *pointflags;			/*  */

  pointflags=shell->pointflags;
  shell->flags |= UNWANTED;

  n=0;
  for (i=0; i< shell->sphere->npoints; i++) 
    if (pointflags[i] & SHARED) {
      warn("SHARED; deleting point %d %s\n", i,
		 atom_spec(fileglobal.box, shell - fileglobal.box->shells));
      n++;
      assert( ! (pointflags[i] & UNWANTED));
      shell->nburied++;
      pointflags[i] &= ~SHARED;
      pointflags[i] |= UNWANTED;
    }
  return n;
} /* shared2unwanted */

static int do_find_patches(shell_p shell, int recover, int *nbndriesp, 
			   int firsttime) {
  /*
   * does    : find all shells patches;
   *
   * gets    : a shell, and a flags that tells wether to save the found points
   *
   * affects : shell->pointflags[*], *nbndriesp will be nr of BOUNDARIES found
   *
   * returns : status
   *
   * warns   : about multiple boundries for a point
   *
   * comment : 
   */
  int i, wanted, start, npatches, nbndries, shared;
  ushort shflags;
  SHORT *pointflags;
  patch_t patcht,  *patch;
  bndry_p bn;
  int directions[ MAX_TMPLT_PNTS ];
  SHORT patchpoints[ MAX_TMPLT_PNTS ];
  geom_t geomt, *geom;   

  /* init */
  free_patches(shell->patches);	/* may be NULL */
  shell->patches=NULL;		/* ! */
  fileglobal.firsttime=firsttime;

  geom= &geomt;
  pointflags = shell->pointflags;
  npatches=nbndries=0;
  shflags=0;
  ngb_geometry(shell, &geomt);		/* local data concerning neigbours */
  assert(! shell->patches); 

  shared=0;
  wanted = shell->sphere->npoints - shell->nburied; /* to be accounted for */
  for ( ; wanted != 0; wanted -= patcht.npoints) { /*one patch per traversal */
    start= -1;				/* meaning: unknown */
    assert(wanted > 0 );
    bzero((char*)&patcht, sizeof(patch_t));
/*    patcht.area=0.0;			 need not be \0\0\0\0 */
    patcht.points= &patchpoints[0];	/* temporary storage */

    while (1) {				/* each traversal finds one bndry */
      patcht.npoints=0;
      bn=find_boundary(shell, start, &patcht, &geomt, directions, &shared);
      assert(bn);			/* ??? */
      nbndries++;
      assert(nbndries < 10);
      bn->next=patcht.bndries;			/* prepend */
      patcht.bndries=bn;
      if ( mark_interior(shell, wanted, &patcht, &directions[0]) )
	/* multiple bndry, i.e., a patch with a 'hole' in it (i.e., 
         * one that runs as a band around the whole atom; fairly common)
         */;
      else 
	break;				/* else: patch has multiple boundry */
      shflags |= MBNDRIES;
      for (i=0; i<patcht.npoints; i++)	/* unmark this patch's points */
	pointflags[ patchpoints[i] ] &= ~PATCH ;
      start = patchpoints[ patcht.npoints-1 ]; /* new starting point */
      /* directions[0] should have been set in mark_interior */
    } /* while 1 */
    npatches++;				/* done collecting all patch bndries */
    patch = save_patch(shell, &patcht, recover);
    for(bn=patch->bndries; bn; bn=bn->next)
      bn->ncons=0;			/* not npoints anymore */
    patch->next=shell->patches; 
    shell->patches=patch;
  } /*  for (wanted != 0) */

  shflags |= (npatches>1)*MPATCHES;
/*  assert( check_shared(shell) == 0); */
  set_conlists(shell);
/*   shflags  |= (mark_shared(shell) >1) *SHARED; */

  *nbndriesp= nbndries;			/* just for statistics */
  shell->flags |= shflags;

  return shared;
} /* do_find_patches */


int find_patches(shell_p shell, int recover, int *nbndriesp, 
		    int firsttime) {
/* for each atom (==shell), finds all the atomic patches. If SHARED
   points are found, these points are first marked UNWANTED, then the
   process is repeated. In rare cases, this has to be done several times. */
  int shared, deleted;
  
  shared=do_find_patches(shell, recover, nbndriesp, firsttime);
  if (shared) { 
    deleted = shared2unwanted(shell);
    clear_top_flags(shell);
    shared=do_find_patches(shell, recover, nbndriesp, 1);
    if (shared) { 
      warn("RARE: can't get rid of SHARED points, now iterating until gone\n");
      do {
        deleted = shared2unwanted(shell);
        warn("got get rid of %d additional SHARED points\n", deleted);
        clear_top_flags(shell);
        shared=do_find_patches(shell, recover, nbndriesp, 1);
        if (shared)
          warn("but found %d additional SHARED points\n", shared);
      } while (shared > 0);
    }
    return deleted;
  } else 
    return 0;
} /* find_patches */

void clear_top_flags(shell_p shell) { 
  /* resets BNDRY and PATCH flags */
  int i;
  SHORT * pointflags;
  
  pointflags=shell->pointflags;

  for (i=0; i<shell->sphere->npoints; i++) 
    pointflags[i] &= ( BURIED | UNWANTED); /* clears all flags but those */
  return;
} /* clear_top_flags */

static int do_widows(shell_p shell) {
  /* checks widows that are due to absence (i.e. being buried) of any */
  /* neigbours; returns nr of deleted points */
  int i, n, u, ntot, p, atomnr;
  patch_p patch, pnext, *wp;
  bndry_p bndry, bnext, *wb;

  ntot=n=0;
  atomnr = shell - fileglobal.box->shells;
  pnext=NULL; bnext=NULL;
  wp= & shell->patches;
  for(patch=shell->patches; patch; patch=pnext) { 
    pnext=patch->next;
    wb = &patch->bndries;
    for(bndry=patch->bndries; bndry; bndry=bnext) {
      bnext=bndry->next;
      if (bndry->flags & ISWIDOW) {	/* delete that */
	n++;
	assert(bndry->ncons==0);
	assert(bndry->path);		/* must be ... */
	if (bndry->length > 3)
	  warn("bndry of length %d becomes WIDOW\n", bndry->length);
	for (i=0; i<bndry->length; i++) {
	  p=bndry->path[i];
	  assert( !(shell->pointflags[p]  & SHARED));
	  u= ( (shell->pointflags[ p ] & UNWANTED) ==0 ); /* newly unwanted */
	  if (u) { 
	    shell->flags |= UNWANTED;
	    shell->pointflags[ p ] |= UNWANTED;
	    shell->nburied ++;
	    ntot ++;
	    warn("no connections; deleted point %d %s\n", p, 
		 atom_spec(fileglobal.box, atomnr));
	  } /* if not already unwanted */
	} /* for i < length */
	free_bndry(bndry);

	/* excise bndry from patch's list of bndries: */
	*wb = bnext;
	patch->flags |= HASWIDOWS;
	bndry=NULL;			/*  */
      }	/* if a widow */ 
      else
	wb = &bndry->next;
    } /* for bndry */
    if (patch->bndries) {
      wp = & patch->next;
      continue;
    }
    /* else: patch has become a widow itself; excise */
    shell->flags |= HASWIDOWS;
    free_patch(patch);
    *wp = pnext;
  } /* for patch */

  if (shell->patches==NULL)
    shell->flags |= ISWIDOW;
  else 
    set_conlists(shell);
  assert(n);
/*  assert(ntot); NO ! */
  return ntot;
} /* do_widows */

#define INVALID ( BURIED | ISOLATED | ISWIDOW )
int remove_widows(shell_p shell, int *deletedp) { /* macroize ? */
  /* 1 if shell cannot be used any more, 0 otherwise */
  if (shell->flags & INVALID )
    return 1;
  if ( ! (shell->flags & REMOVE_WIDOWS) ) /* apparently no work to be done */
    return 0;
  shell->flags &=  ~REMOVE_WIDOWS;
  *deletedp += do_widows(shell);
  return  shell->flags & INVALID;
} /* remove_widows */

static int reset_next_ptrs(box_p box) {
  int i;
  /* make next.next pointer NULL, because they'll become next.conlists */
  for (i=0; i<box->natoms; i++) 
    box->shells[i].next.next=NULL;
  return 0;
} /* reset_next_ptrs */


int atom_topology(box_p box, int mode) {
  /* establishes patches on atoms; returns total nr of point deleted */
  int nbndries, ii, recover, deleted;
  shell_p sh;

  deleted=0;
  fileglobal.box=box;
  recover=mode-1;

  reset_next_ptrs(box);

  /* establish topology: */
  for (ii=0; ii<box->natoms; ii++) {
    sh=box->shells+ii;
    if ( remove_widows(sh, &deleted) )
      continue;
    deleted +=
      find_patches(sh, recover, &nbndries, 1);

    remove_widows(sh, &deleted);
  } /* for ii */
  return deleted;
} /* atom_topology */

