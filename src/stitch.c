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

/* routines for stitching the seams between patches */

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

static struct  { 
  box_p box;
  shell_p *queue;
  int qsize;
  int nwidows;
  int nunsuited;
} local;

/* flags for br status:  */
#define SAME_BIT	0 
#define SAME 		BIT(SAME_BIT)
#define ALREADY 	BIT(1)

#define PHI0_BIT 	2
#define PHI1_BIT 	3
#define PHI2_BIT 	4

#define PHI0		BIT(PHI0_BIT)
#define PHI1	 	BIT(PHI1_BIT)
#define PHI2	 	BIT(PHI2_BIT)
#define PHI		(PHI0 | PHI1 | PHI2)

#define SUB0_BIT 	5
#define SUB1_BIT 	6
#define SUB0 		BIT(SUB0_BIT)
#define SUB1 		BIT(SUB1_BIT)
#define SUB	  	(SUB0 | SUB1)

#define TOOLONG_BIT 	7
#define TOOLONG 	BIT(TOOLONG_BIT)

#define DUB_BIT		8
#define DUB		BIT(DUB_BIT)

#define VERYDUB_BIT	(DUB_BIT + 1)
#define VERYDUB		BIT(VERYDUB_BIT)

#define WRONG (SAME | ALREADY | PHI | SUB | TOOLONG | VERYDUB)

/* following has to be defined: */

static void queue_shell(shell_p shell) { 
  /* inserts shell at end of queue, if it has to be reexamined */
  if (shell->flags & STITCHED) {
    shell->flags &= ~STITCHED;		/* mark as done */
    local.queue [ local.qsize++ ] = shell; /* add to queue */
  }  /* else: is already in queue */
  return;
} /* queue_shell */
/* the matching dequeue_shell() can be found just before stitch() */


const float * bndry_coordinates(bndry_p bndry) {
  /*
   * does    : if not already allocated, allocates bndry->ngbs.coordinates
   *   as well as bndry->pt.normals, and sets them; returns coordinates;
   *
   * gets    : 
   *
   * affects : 
   *
   * returns : 
   *
   * warns   : 
   *
   * comment : 
   */
  int i, k, length, p;  const atom_t *atom;
  SHORT * path;
  float * coordinates, R, *normals, x;
  const float * atomxyz, *spherepoints;

  if(bndry->ngbs.coordinates)
    return bndry->ngbs.coordinates;

  atom=bndry->patch->atom;
  atomxyz=atom->xyz;
  R=atom->RADIUS;

  spherepoints=((shell_p)atom->pointer)->sphere->coordinates;

  length=bndry->length;
  path=bndry->path;
#define SIZE sizeof(float[3])
  bndry->ngbs.coordinates=coordinates=(float *)CALLOC(length, SIZE);
  bndry->pt.normals=normals=(float *)CALLOC(length, SIZE);
#undef SIZE
  for (i=0; i<length; i++) {
    p=path[i];				/* id of point along path */
    for (k=0; k<3; k++) { 
      x=spherepoints[ 3*p +k ];
      normals[3*i + k] = x;
      coordinates[ 3*i + k ] = atomxyz[k] + R*x; 
    }
  }
  bndry->flags |= COORDINATES;
  return bndry->ngbs.coordinates;
} /* bndry_coordinates */

static void free_coordinates(bndry_p bndry) {
  if(bndry->ngbs.coordinates) {		/* don't do if already NULL */
    assert( bndry->flags & COORDINATES ); /* otherwise no sense */
/*    bndry->flags &= ~COORDINATES; */
    assert(bndry->pt.normals);		/* they come in pairs */
    FREE(bndry->ngbs.coordinates);
    bndry->ngbs.coordinates=NULL;
    FREE(bndry->pt.normals);
    bndry->pt.normals=NULL;
    return;
  } /* else: */
  assert(! bndry->pt.normals);	/* should have been deleted */
  return;
} /* free_coordinates */

static int delete_crossing_cons(int nseams, seam_p seams, int maxnseams) {
  /* deletes all crossing edges from seam. queues all shells affected */
  int i, ncross, end;
  con_p con, *crossing;

  crossing=(con_p*)ALLOCA(maxnseams/2 * sizeof(con_p));

  ncross=0;
  for (i=0; i<nseams; i++) {		/* find crossings: */
    con=seams[i].con;
    if(con->flags[0].seen==1 && con->flags[1].seen==1)
      crossing[ ncross++ ]=con;		/* keep for deletion after this */
    else {
      seams[i].con->flags[0].seen=0;	/* reset flags on the others */
      seams[i].con->flags[1].seen=0; 
      end= seams[i].end;
      seams[i].con->flags[ end ].done = 0; /*  */
    }
  } /* for i < nseams */

  for (i=0; i<ncross; i++) {		/* and delete them */
    con=crossing[i];
    if(con->flags[0].seen==1)		/* every con appears twice ! in list */
      con->flags[0].seen=0;		/* skip it the first time  */
    else {
      queue_shell( con->bndries[0]->patch->atom->pointer );
      queue_shell( con->bndries[1]->patch->atom->pointer );
      remove_con(con);
    }
  } /* for i < ncross */
  AFREEA(crossing);
  return ncross/2;
} /* delete_crossing_cons */

static /* double */ void 
  set_dir( double *dir, const float *y, const float *x) {
  /* sets dir to be unit vector from y to x; dir[3] is set to the distance */
  /* returns 0 */
  int k;
  double d, d2;

  if (y==x) {
/*    (dir[3]= HUGE_VAL); */
    assert(0);
  }

  d2=0.0;
  for (k=0; k<3; k++) {
    d = dir[k]=y[k]-x[k];
    d2 += SQR(d);
  }
  d= sqrt(d2);
  dir[3]=d;
  for (k=0; k<3; k++)
    dir[k] /= d;
/*  return d; */
} /* set_dir */

static br_p closest_bridge(int nbridges, br_p bridges) {
  /* finds the next closest bridge in this ring of bridges */
  int i;
  uint flags;
  br_p br, minbr;
  double d, min;

  assert(nbridges >= 4);
  min=HUGE_VAL;				/* first look for minimum distance */
  minbr=NULL;

  for (i=0,br=bridges; i<nbridges; i++,br=br->next) {
    flags=br->flags;
    if( flags & WRONG )
      continue;
    d=br->newdir[3];
    d *= (1 + ((flags & DUB)!=0));
    if (d < min) { 
      min=d;
      minbr=br;
    }
  }
  return minbr;			/* may be NULL */
} /* closest_bridge */

static con_p new_con(br_p br) {
  int i, end;
  br_p nextnext;
  double sign;
  con_p con;

  con=alloc_con(NULL);

  nextnext=br->next->next;
  end=(br->bndry->patch->atom > nextnext->bndry->patch->atom );
  con->bndries[ end ]=br->bndry;
  con->bndries[ !end ]=nextnext->bndry;
  con->flags[ end ].offset=br->offset;
  con->flags[ !end ].offset=nextnext->offset;

  sign= 1.0 - 2.0*end;

  for (i=0; i<3; i++) 
    con->direction[i]= sign* br->newdir[i];
  con->direction[3]= br->newdir[3] ; 
  con->flags[0].done=1;
  con->flags[1].done=1;
  return con;
} /* new_con */

static uint normals_flags(br_p br) {
  /* checks if candidates are not at funny angles with points they */
  /* connect, and returns flags to model this */
  double *dir, d1, d2;
  const float *n1, *n2;

  dir=br->newdir;
  n1=br->normal;
  n2=br->next->next->normal;

  d1=DOT3P(n1, dir);
  br->dotp[0]=d1;
  d2= -DOT3P(n2, dir);
  br->dotp[1]=d2;

  return (SUB0*(d1 < br->mindotp)) | (SUB1 * (d2 < br->next->next->mindotp));
} /* normals_flags */

static double torsion(br_p cand) {
  /* returns dotproduct of normals at and points of current cand */
  const float *n1, *n2;

  n1=cand->normal;
  n2=cand->next->next->normal;
  return DOT3P(n1, n2);
} /* torsion */


static int same_atom(br_p br) {
  return 
    (br->bndry == br->next->next->bndry) /* same bndry */
      ||				/* or same atom */
	(br->bndry->patch->atom == br->next->next->bndry->patch->atom);
} /* same atom */

static int test_angle(double* x, double*y, const float *normal) {
  /* returns flags based on wether connection lies within current vectors */
  double c[3], d;
  static double MINDOTP = 0.0, MINMINDOTP= -0.2;

/* looks like:
 *  n  y  
 *   +--->
 *   ^  ^ 
 * x | / c
 *   |/ 
 *   .
 */
    
  CROSS3P(c, y, x);
  d=DOT3P(normal, c);
  if(d < MINDOTP )
    if (d < MINMINDOTP )
      return -1;			/* really too small */
    else
      return 1;				/* dubious */
  else 
    return 0;
} /* test_angle */

static int point_already(br_p br, br_p obr) { 
  /* checks if end point already present among connections */
  int end;
  SHORT /*path,*/ *opath, opoint;
  list_p l;
  bndry_p bndry, obndry;
  con_p con;

  bndry=br->bndry;
/*  path=bndry->path; */
  opath=obr->bndry->path;
  opoint=obr->bndry->path[ obr->offset];
  for (l=bndry->cons[br->offset]; l; l=l->next) {
    con=l->ptr;
    end=(con->bndries[0]==bndry);	/* other end */
    obndry=con->bndries[end];
    if (obndry->patch->atom->pointer != obr->shell)
      continue;
    if ( obndry->path[ con->flags[end].offset ] == opoint )
      return 1;
  }
  return 0;
} /* point_already */

static int con_already(br_p br) {
  /* checks if con is not already represented in one way or other */
  br_p obr;

  obr=br->next->next;
  return point_already(br, obr) || point_already(obr, br);
} /* con_already */

static con_p con_tooclose(br_p br, double*newdir, int sign) {
  /* walks lists of cons on other bndry passing through current point */
  int s, n;
  list_p l, **conlists, *conlist;
  shell_p shell;
  bndry_p bn;
  con_p con;
  double dot;
  static double t=0.8;

  if ( ! (*br->pointflags & (SHARED | TWOWAY | UTURN | ISOLATED)))
    return NULL;

  conlists=br->conlists;
  assert(conlists);

  /* walk list of all conlist available for this point: */
  n=0;
  while( (conlist=conlists[n++]) ) {
    l= *conlist;
    if (l==NULL)			/* still empty */
      continue;
    con=l->ptr;
    bn=con->bndries			/* find out this bndry's id */
      [ (con->bndries[0]->patch->atom->pointer != shell) ];
    for (; l; l=l->next) {
      /* avoid cons immediately adjacent: */
      if( &l->next == br->where)	/* preceeding con */
	continue;
      if (*br->where == NULL)		/* end of list */
	continue;
      if (*br->where == l)		/* next con  */
	continue;
      con=l->ptr;
      s = (2*(con->bndries[0]==bn)-1)*sign;
      dot = ((double)s)*DOT3P(newdir, con->direction);
      if (dot > t)
	return con;
    } /* for l */
  } /* while conlist */
  assert(n<=4);
  return NULL;
} /* con_tooclose */

static int check_ply(br_p br) { 
  con_p con;

  if( (con=con_tooclose( br, br->newdir, 1)) )
    return 1;
  if( (con=con_tooclose( br->next->next, br->newdir, -1)) ) { 
    return 2;
  }
  return 0;
} /* check_ply */

static void cand_geom(br_p br /*, int nbr */) {
  /* determines admissability of, and geom date of a new candidate bridge */ 
  int i, dub;
  uint f, flags;
  double d;
 static double MINTORSION= -0.8;

#define BITSHIFT(word, bit, amount) (((word) & (BIT(bit)) )<<(amount))

#define DUB_ORTWICE(word, cond) ((DUB*dub) | BITSHIFT((word), DUB_BIT, (cond)))

  if (same_atom(br)) { 
    br->flags = SAME; 
    return; 
  } 
  if ( con_already(br) ) { 
    br->flags = ALREADY;
    return; 
  }

  flags=0;
  i=test_angle(br->dir, br->next->dir, br->next->normal);
  dub= (i > 0);
  flags |= DUB_ORTWICE(flags, dub);
  flags |= (PHI0 * (i<0) );

  /* d= */ set_dir( br->newdir, br->next->next->coord, br->coord);

  f=normals_flags(br);			/* no dubious here */
  flags |= f;

  if ((i=check_ply(br)))
    flags |= VERYDUB;

  d=torsion(br);
  dub= (d < MINTORSION);
  flags |= DUB_ORTWICE(flags, dub);
  br->flags = flags;
  return;
} /* cand_geom */

static void set_seams(int nseams, seam_p seams) {
  /*
   * does    : makes real seams out of the data in and ends, writing it
   *   to seams. Mainly sets from and till pointers
   *
   * gets    : array of cons and ends that constitute the seam
   *
   * affects : seams[*], cons[*].done
   *
   * returns : -
   *
   * warns   : 
   *
   * comment : this should be optimized a bit more ...
   */
  int i, n, end, oend;
  con_p con;

  for (i=0; i<nseams; i++) {
    n=(i+1)%nseams;			/* index of next connection */
    end= seams[i].end;
    con= seams[i].con;
    con->flags[ end ].seen=0;
/*    con->flags[0].seen=0;  */
/*    con->flags[1].seen=0;  */

    seams[i].from= con->flags[ end ].offset;
    oend=end;
    end = !seams[n].end;		/* begin of the next seam */
    assert( seams[n].con->bndries[ end]== con->bndries[oend] ); /* of course */
    seams[i].till = seams[n].con->flags[ end ].offset;
  } /* for i < nseams */

  return;
} /* set_seams */

static list_p * insert_place(const seam_t * seam, int offset) {
  list_p l;
  con_p con;

  if ( offset == seam->from ) {
    if(seam->length==1) {		/* special case */
      con=seam->con;
      for (l=seam->bndry->cons[offset]; l; l=l->next) 
	if(l->ptr==con)
	  return &l->next;
      assert(0);			/* not found: can't be */
    } else {
      for (l=seam->bndry->cons[offset]; l->next; l=l->next) 
	;
      return &l->next;
    }
  } /* else */				/* simple */
  return seam->bndry->cons + offset;
} /* insert_place */

static double find_dotp(bndry_p bn) {
  /* returns a mininmal dotproduct for this bndry */
  sphere_p sph;
  
  sph= ((shell_p)bn->patch->atom->pointer)->sphere;
  return -sph->sinbeta*0.8;
} /* find_dotp */

static int delete_point(br_p br) { 
  /* marks point as UNWANTED; returns 1 if it wasn't already, 0 otherwise */
  /* queues shell */
  int p, newly, atomnr;
  SHORT flags; 
  shell_p shell;

  flags= *br->pointflags;
  assert( flags & BNDRY);		/* only deleting bndry points */

  shell=br->shell;
  *br->pointflags |= UNWANTED;
  shell->flags |= UNWANTED;
  newly=( (flags & UNWANTED) ==0);	/* often happens more than once */
  shell->nburied += newly;		/*  !!! */

  shell->flags |= (UNWANTED | RETOP);
  queue_shell(shell);

  if(newly) { 
    p=br->bndry->path[ br->offset ];
    atomnr=shell - local.box->shells;

    warn("not suited: deleted point %d %s\n", p
	 , atom_spec(local.box, atomnr) );
  }
  return newly;
} /* delete_point */

static int init_bridges(int nseams, seam_p seams, br_p bridges, 
			int *deletedp) {
  int i, j, nbr, len, totlen, blen, offset, deleted, p;
  br_p br;
  shell_p shell;
  seam_p seam;
  bndry_p bndry;
  const float *coordinates, *normals;
  SHORT *path;
  double dotp, atomradius;

  set_seams(nseams, seams);		/* sets from and till pointers */
  *deletedp=deleted=totlen=nbr=0;
  br=bridges;
  for (i=0; i<nseams; i++) {		/* don't forget the from and till */
    seam= &seams[i];
    len=seam->length;
    totlen += len;			/* keep count of total length */
    bndry=seam->bndry;
    shell=bndry->patch->atom->pointer;
    atomradius= bndry->patch->atom->RADIUS;

    dotp=find_dotp(bndry);
    blen=bndry->length;
    coordinates=bndry_coordinates(bndry);
    normals=bndry->pt.normals;		/* also set in bndry_coordinates */
    path=bndry->path;
    for (j=0; j<len; j++) {
      br->next=br+1;
      br->prev=br-1;
      br->bndry=bndry;
      br->mindotp= dotp;
      br->shell=shell;
      br->atomradius=atomradius;
/*      br->seam=seam; */
      offset=(seam->from+j)%blen;
      br->where=insert_place(seam, offset);
      br->offset=offset;
      br->coord=coordinates + 3*offset;
      br->normal=normals + 3*offset;
      p=path[offset];
      br->conlist= &bndry->cons[offset];
      br->conlists=shell->next.conlists[ p ];
      br->pointflags = &shell->pointflags[ p ];
      if ( (*(br->pointflags)) & BURIED) /* seen before ! (no other way?) */
	 deleted += delete_point(br);
      *(br->pointflags) |= BURIED;	/* mark */
      br++;
      nbr++;
    } /* for j < len */
  } /* for i <nseams */

  if (deleted) {
    for (i=0; i<nbr; i++) 
      *(bridges[i].pointflags) &= ~BURIED; /* unmark */
    *deletedp = deleted;
    return nbr;
  }

  assert(br - bridges == totlen && nbr==totlen);
  br--;					/* one bridge too far ... */
  bridges->prev=br;			/* circularize */
  br->next=bridges;

  for (i=0; i<nbr; i++) {
    *(br->pointflags) &= ~BURIED;	/* unmark */
    br=bridges+i;
    set_dir(br->dir, br->next->coord, br->coord); /* directions along seam */
  }
  for (i=0; i<nbr; i++)
    cand_geom(bridges+i /*, nbr*/);	/* directions and weights of */
					/* candidate bridges */
  assert(nbr==totlen);			/* (may change) */

  return totlen;
} /* init_bridges */

static void unmark_sides(con_p con, int side) { 
  /* resets the 'done' bit of first con on the 'side' side of this one */
  int offset, len, i, j, end;
  bndry_p bndry;
  
  bndry=con->bndries[ side ];
  len=bndry->length;
  offset=con->flags[ side ].offset;

  for (i=1; i<len; i++) {		/* or <= len ? */
    j=(offset+i)%len;
    if ( bndry->cons[j] ) {
      con=bndry->cons[j]->ptr;
      end= (con->bndries[0]==bndry);	/* (hopefully got them right) */
      con->flags[end].done=0;
      return;
    } /* if */
  } /* for i  */
  return;
} /* unmark_sides */

static int delete_point_cons(bndry_p bndry, int offset) { 
  /* deletes all cons from and to this point; returns nr of points deleted */
  /* queues shell on both sides (although one should be enough ?)  */
  int n;
  list_p l, next;
  con_p con;

  n=0;
  for (l=bndry->cons[offset]; l; l=next) { /* delete all con's */
    n++;
    next=l->next;
    con=l->ptr;
    unmark_sides(con, 0);
    unmark_sides(con, 1);
    queue_shell(con->bndries[0]->patch->atom->pointer);
    queue_shell(con->bndries[1]->patch->atom->pointer);
    remove_con_halfknown(l, bndry);
  }
  bndry->cons[ offset ]=NULL;
  return n;
} /* delete_point_cons */

static int add_conlist(shell_p shell, bndry_p oldbndry, int offset) {
  /* puts conlist at offset in oldbndry into newbndry; returns nr of them */
  int i, end, n, point;
  patch_p patch;
  bndry_p bndry;
  list_p l, conlist;
  con_p con;

  n=0;
  point=oldbndry->path[offset];
  conlist=NULL;
  assert( !(shell->pointflags[point] & (SHARED | TWOWAY | BURIED | UNWANTED)));
  /* search place to insert:  */
  for(patch=shell->patches; patch; patch=patch->next) /* NEW patches */
    for(bndry=patch->bndries; bndry; bndry=bndry->next) { /* NEW bndries */
      assert( !  (bndry->flags & ISWIDOW));
      for (i=0; i<bndry->length; i++)
	if (bndry->path[i]==point)	/* gotcha */
	  goto END;			/* threefold break */
    } /* for patches, bndries */
  assert(0);				/* not found ... */
 END:
  assert(bndry->cons[i] == NULL );	/* must be first encounter */
  conlist=bndry->cons[i] = oldbndry->cons[offset];
  oldbndry->cons[offset]=NULL;

  /* assign the new bndries + offsets: */
  for (l=conlist; l; l=l->next) {
    n++;
    con=l->ptr;
    end= (con->bndries[0] != oldbndry);	/* this end */
    con->bndries[ end ]=bndry;		/* differs from old bndry ! */
    con->flags[ end ].offset=i;		/* offset usually different ! */
  }
  bndry->ncons += n;
  return n;
} /* add_conlist */

static int transfer_cons(shell_p shell, patch_p old_patches, 
			 SHORT * old_pointflags) {
  /* re-establishes connections stuff as far as possible; deletes rest */
  int i, point, m, n;
  patch_p oldpatch;
  bndry_p oldbndry;

  n=m=0;
  for(oldpatch=old_patches; oldpatch; oldpatch=oldpatch->next) {
    for(oldbndry=oldpatch->bndries; oldbndry; oldbndry=oldbndry->next) {
      assert( !(oldbndry->flags & ISWIDOW));
      for (i=0; i<oldbndry->length; i++) {
	point=oldbndry->path[ i ];
	if (oldbndry->cons[i] == NULL)
	  continue;
	if ( (old_pointflags[ point ] | shell->pointflags[ point ]) 
	    & (UNWANTED | SHARED | TWOWAY ) ) {	/* throw away */
	  m += delete_point_cons(oldbndry, i); /* all cons FROM and TO point */
	} else				/* keep */
	  n += add_conlist(shell, oldbndry, i);
      }	/* for i */
    } /* for bndry */
  } /* for patch */
/*  warn("accomodated %d, rejected %d connections\n", n, m); */
  return n;
} /* transfer_cons */

static int mark_widows(shell_p shell) { 
  /* judges wether bndries are WIDOW on the basis of their nr of new */
  /* connections; returns nr of widow bndries */
  int n;
  patch_p pa;
  bndry_p bn;

  n=0;
  for(pa=shell->patches; pa; pa=pa->next)
    for(bn=pa->bndries; bn; bn=bn->next)
      if (bn->ncons==0) { 
	n++;
	bn->flags |= ISWIDOW;
	shell->flags |= REMOVE_WIDOWS;
      }
  assert(n==0 || shell->flags & REMOVE_WIDOWS);
  return n;
} /* mark_widows */

static int redo_topcon(shell_p shell) { 
  /* does everything needed to restore this shells patches etc. as part of */
  /* the triangulation; returns nr of widow points */
  int i, nwidows;
  char *old_pointflags;
  int n, dummy;
  patch_p old_patches;

  nwidows=0;
  assert(shell->flags & RETOP);
  old_patches=shell->patches;
  shell->patches=NULL;			/* otherwise find_patches frees it */
  /* preserve pointflags, since SHARED + TWOWAY are cleared: */
  MEMSAV(old_pointflags, (char*)shell->pointflags,
	 shell->sphere->npoints*sizeof(short));	/* (alloc on stack) */

  clear_top_flags(shell);		/* throw away everything but BURIED */

  i=find_patches(shell, 0, &dummy, 0); /* new topology */
  /* what to do with deleted points ??? */
  free_neibordata(shell);

  assert( ! (shell->flags & (ISWIDOW | REMOVE_WIDOWS) ));
  n=transfer_cons(shell, old_patches, (SHORT*)old_pointflags);
  shell->flags &= ~(ISWIDOW | REMOVE_WIDOWS); /* pertains to old flags ! */

  mark_widows(shell);
  remove_widows(shell, &nwidows);

  free_patches(old_patches);
  AFREEA(old_pointflags);
  return nwidows;
} /* redo_topcon */

static int delete_usual_suspects(int nbridges, br_p bridges) {
  int i, n;
  SHORT f, allflags;
  br_p br;

  allflags=0;
  for(i=0, br=bridges; i<nbridges; br=br->next,i++)
    allflags |= *br->pointflags;

#define SETFLAGS(FLAGS) if (allflags & (FLAGS) ) f=(FLAGS);


  if_not (allflags & ( SHARED | UTURN | ISOLATED | TWOWAY) )
    return 0;

  SETFLAGS(SHARED) else SETFLAGS(UTURN) else SETFLAGS(ISOLATED) 
    else SETFLAGS(TWOWAY);
  assert(f==SHARED || f==UTURN || f==ISOLATED || f==TWOWAY);

  n=0;
  for (i=0, br=bridges; i<nbridges; br=br->next,i++)
    if ( (*br->pointflags) & f)
      n+=delete_point(br);

  assert(n);
  return n;				/* may be 0: only VERYDUB's found */
} /* delete_usual_suspects */

static int delete_points(int nbridges, br_p bridges) { 
  /* decides which points are to be deleted, sets flags */
  /* accordingly. Returns nr of points deleted */
  int n;
  uint flags;
  br_p br;
  int i;

  n=delete_usual_suspects(nbridges, bridges);
  if (n)
    return n;

  for (i=0,br=bridges; i<nbridges; i++,br=br->next) {
    flags=br->flags;
    if (flags & PHI0)
      n += delete_point(br->next);
    
    assert( ! (flags & (PHI1 | PHI2)));

    if (flags & SUB0)
      n += delete_point(br->next->next);
    if (flags & SUB1)
      n += delete_point(br);
    
  } /* for i < nbridges */
  if (n)
    return n;
  /* else: not found; only VERYDUB available; clear them */
  warn("clearing VERYDUB flags\n");
  for (i=0,br=bridges; i<nbridges; i++,br=br->next)
    br->flags &= ~VERYDUB;
  return 0;
} /* delete_points */

static int repair(int nbridges, br_p bridges) { 
  /* re-establishes topology of atoms where deletions occurred; returns nr */
  /* of points deleted (due to becoming widows) */
  /* queues shells that had to be redone */
  int i, nwidows;
  br_p br;
  shell_p shell;

  nwidows=0;				/* remove_widows increments ! */
  for (i=0,br=bridges; i<nbridges; i++, br=br->next)  {
    shell=br->shell;
    if ( remove_widows(shell, &nwidows) )
      continue;
    if (shell->flags & RETOP) { 
      nwidows += redo_topcon(shell);	/* new topology and connections */
      shell->flags &= ~RETOP;
      queue_shell(shell);
    } /* if RETOP */		        /* deletion of connections ?? */
  } /* for i < nbridges */
  return nwidows;
} /* repair */

static int triangulate_seam(int nseams, seam_p seams, int * widows) {
  /* triangulates interior of seam; returns nr of points deleted */
  int k, n, nbridges, deleted, maxnbridges, nwidows;
  br_t *bridges, *br, *nextnext, *newbr;
  con_p con;
  list_p *where;

  maxnbridges=nseams*MAX_PATHLENGTH;
  bridges = (br_p)ALLOCA(maxnbridges * sizeof(br_t));

  nwidows=0;
  n=nbridges=init_bridges(nseams, seams, bridges, &deleted);
					/* now a circular list */
  assert(nbridges < maxnbridges);

  if (deleted)				/* already some deleted */
    nwidows=repair(nbridges, bridges);
  else {				/* normal case */
    for(br=bridges; 1; n-- ) {		/* exit is via a break  */
      newbr=closest_bridge( n, br);
      if(newbr==NULL) {			/* could not find accetable br */
	deleted=delete_points(n, br);
	if (deleted) {
	  nwidows=repair(n, br);	/* re-does topology and connections */
	  break;			/* out of for loop */
	} /* else: VERYDUB flags are reset; try again */
	newbr=closest_bridge(n, br);
	assert(newbr);
      }
      br=newbr;
      con=new_con(br);
      br->bndry->ncons++;
      br->next->next->bndry->ncons++;

      /* insert con into bndries: */
      where = br->where; 
      *where=Lcons(con, *where);	/* br->where unchanged */

      where=br->next->next->where;
      *where=Lcons(con, *where);

      if(n==4)				/* ready ! */
	break;

      nextnext=br->next->next;
      nextnext->where= &((*where)->next); /* this one does change */

      /* rewire (effectively shortcutting br->next): */
      br->next=nextnext;
      nextnext->prev=br;
      for (k=0; k<4; k++)
	br->dir[k] = br->newdir[k];

      /* following weigths and dir distances/directions change: */
      cand_geom(br);			/* this order !  */
      cand_geom(br->prev);
    } /* for n >= 4 */
  } /* if deleted==0 */
  
  AFREEA(bridges);
  *widows=nwidows;
  return deleted;
} /* triangulate_seam */

static int do_seam(con_p con, int end, 
		   int maxnseams, seam_p seams, int *widows) {  
  /* triangulates the interior of the seam started by con in
   * direction end (ie. con->bndries[end] is the one we are not on). returns
   * 0 if all OK, or -1 if con's have been taken out due to
   * crossing. Otherwise, nr of shells affected during repair action is
   * returned 
   */

  int nseams, n;

  n=find_seam(con, end, &nseams, maxnseams, seams);
  switch(n) {
  case 0: 
    break;				/* found a triangle: done */
  case 1:				/* length > 2  (normal case) */
    n=triangulate_seam(nseams, seams, widows); /* triangulate inner */
    break;
  case 4:				/* con crossing */
    n=delete_crossing_cons(nseams, seams, maxnseams);
    n= -n;				/* negative */
    break;
  default:
    assert(0);				/* should not occur */
  } /* switch */

  return n;
} /* do_seam */

static int stitch_shell(shell_p shell, int maxnseams, seam_p seams, 
			int * widows) { 
  int i, end, n;
  patch_p patch;
  bndry_p bndry;
  list_p l;
  con_p con;

  assert(! ( shell->flags & STITCHED )); /* completely done with */

  for(patch=shell->patches; patch; patch=patch->next)
    for(bndry=patch->bndries; bndry; bndry=bndry->next) {
      assert (! (bndry->flags & ISWIDOW));
      for (i=0; i<bndry->length; i++) {
	for (l=bndry->cons[i]; l; l=l->next) {
	  con=l->ptr;
	  assert(con->bndries[0] && con->bndries[1]);
	  end=(con->bndries[0]==bndry); /* id of the bndry we're NOT on */
	  if(con->flags[end].done==0) { /* found new untraversed seam */
	    n=do_seam(con, end, maxnseams, seams, widows);
	    if(n != 0 || *widows)	/* < 0:crossing, > 0: points deleted */
	      return n;
	  } /* if done==0 */
	} /* for l=cons[i] */
      } /* for i < length */
      free_coordinates(bndry);
    } /* for patch, bndry */
  for(patch=shell->patches; patch; patch=patch->next) /* free bndries on  */
    for(bndry=patch->bndries; bndry; bndry=bndry->next)	/* original shell */
      free_coordinates(bndry);
  return 0;
} /* stitch_shell */

static void dequeue_shell(shell_p shell) { 
  /* removes shell from end of queue. only to be used from within stitch() */

  assert(local.queue [ local.qsize-1 ] == shell);

  shell->flags |= STITCHED;		/* mark as done */
  local.queue [ --local.qsize ] = NULL; /* remove from queue */
  return;
} /* dequeue_shell */

int stitch(box_p box) {
  int ii, err, n, ntot, ncons, maxnseams, nwidows, w;
  shell_p shell, *queue;
  seam_t *seams;
  
  local.box=box;

  local.queue = (shell_p*)ALLOCA(2*box->natoms*sizeof(shell_p));
  queue=local.queue;

  /* fill list of shells to be done:  */
  nwidows=n=0;
  for (ii=0; ii<box->natoms; ii++) { 
    shell=box->shells + ii;
    if ( remove_widows(shell, &nwidows) ) 
      continue;
    queue[n++]=shell;
  }
  local.qsize=n;			/* nr of non-buried atoms */
  maxnseams= MAX(32, n/8);		/* wild guess */
  seams=(seam_p)ALLOCA( maxnseams*sizeof(seam_t));

  ntot=ncons=0;
  do {					/* drain list till ready */
    assert( local.qsize < 2*box->natoms);
    n=local.qsize;			/* keep for comparing */
    shell=queue[ n-1 ] ;		/* always the last of the list */
    if (shell->flags & STITCHED ) {	/* already done */
      dequeue_shell(shell);
      continue;
    }
    
    if ( remove_widows(shell, &nwidows) ) { /* has become unavailable */
      dequeue_shell(shell);
      continue;
    }

    w=0;
    err=stitch_shell(shell, maxnseams, seams, &w); /* call */

    if (err==0 && w==0 ) { 
      dequeue_shell(shell);
      continue;				/* at do */
    }
    if(err < 0)			/* connections deleted */
      ncons -= err;			/* since negative */
    else
      ntot += err;
      
    nwidows += w;
    if (err>0)
      warn("unable to close gap; deleted %d points around %s\n",
	   err, atom_spec(box, shell - box->shells));
    n=local.qsize - n;
    if (n)
      warn("  (an extra %d atoms got involved)\n", n);
  } while_not (local.qsize == 0);

  warn("deleted %d points (%d unsuitable, %d widow), %d connections\n", 
       ntot+nwidows, ntot, nwidows, ncons);

  AFREEA(seams);
  AFREEA(local.queue);
  return 0;
} /* stitch */

