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
 * Philip Lijnzaad, plijnzaad@gmail.com			      *
 * ------------------------------------------------------------------ */

/* routines to perform numerical surface calculation */

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

static void gridify(box_p box) { 
		   
  /*
   * does    : allocates a 3D block of somepointer, and assigns each
   *   atom to a (possibly empty) linked list of atoms belonging to a
   *   cell of the 3D grid that constitues the block
   *
   * gets    : atoms, their number, and the dimensions of the protein
   *
   * affects : box->ncells, box->cells, shell->next pointers. cells are
   *   laid out like a straight 3D array  
   *
   * warns   : box and cell size, and average and maximum number of atoms
   *   per box. 
   *
   * comment : 
   */
  int i, j, ncells, dims[3], n, offset, size, natoms;
  float maxrad, spacing, *limits, f;
  shell_p  * cells;
  shell_p shell;
  atom_p atom, atoms;
  
  limits=box->limits;
  /* establish dimensions etc. */
  maxrad=limits[ MAXRADIUS ]; 
  spacing = 2.0 * maxrad;

  warn("grid cell size is %.3f Angstrom\n", spacing);

  ncells=1;
  for (i=0; i<3; i++) { 
    f=(limits[2*i +1] - limits[2*i])/spacing;
    n=(int)floor(f) + 1;		/* if eg. alway z==0, still nz==1 */
    dims[i]=n;
    ncells *=n;
  }

  warn("dimensions of box holding protein: %d x %d x %d\n", 
       dims[0], dims[1], dims[2]);

  size=ncells*sizeof(shell_p);
  cells=box->cells=(shell_p*)AREALLOC(box->cells, size);
  bzero((char*)cells, size);

  /* put in box: */
  natoms=box->natoms;
  atoms=box->atoms;

  for (i=0, atom=atoms; i<natoms; i++,atom++) {
    offset=0;
    for (j=0; j<3; j++)  {
      n=(int)floor( (atom->xyz[j] - limits[2*j])/spacing);
      offset = offset*dims[j] + n;
    }
/*  is same as
 *    icoord[j]=ffloor((Atoms[i]->xyz[j] - limits[2*j])/spacing);
 *    offset=dims[2]*dims[1]*icoord[0] + dims[2]*icoord[1] + icoord[2]; 
 */
    assert(offset < ncells);
    shell=atom->pointer;
    shell->next.next=cells[offset];
    cells[offset]=shell;
  } /* for i < natoms */

  /* transfer stuff to box */

  box->cells=cells;
  box->ncells=ncells;
  box->spacing=spacing;
  box->maxrad=maxrad;

  memcpy(box->dims, dims, sizeof(int[3]));

  return;
} /* gridify */


static void purge_ngbs(box_p box) {
  /*
   * does    : purges neighbour lists of neighbours that have become BURIED
   *  
   * gets    : a box
   *
   * returns : 
   *
   * warns   : about isolated atoms and widow atoms
   *
   * comment : 
   */

  int i, nburied, nngbs, j, k, natoms;
  shell_p sh, *ngbs;
  atom_p atoms;

  natoms=box->natoms;
  atoms=box->atoms;

  nburied=0;				/* keep totals */
  for (i=0; i<natoms; i++)  {
    sh=atoms[i].pointer;
      
    if (!sh) {				/* buried: nothing to be done */
      nburied++;
      continue;
    }

    /* count unburied neibors: */
    nngbs=sh->nngbs;
    ngbs=sh->ngbs;
    k=0;
    for (j=0; j<nngbs; j++)		/* get rid of buried neibors */
      if_not (ngbs[j]->flags & BURIED)
	ngbs[k++]=ngbs[j];
    assert(k<=sh->nngbs);
    
    if(k) {				/* still some neigbours left:  */
      sh->nngbs=k;			/* new nr of neibors */
      continue;
    }

    /* at this point we have something funny: */

    if (sh->nburied==0 || sh->flags & ISOLATED) { 
      sh->flags |= ISOLATED;
      warn("atom isolated %s\n", atom_spec(box, i) );
    } else {
      sh->flags |= ISWIDOW;		/* no non-BURIED neighbours */
      warn("atom has %d widow points; considered buried\n",
	   sh->sphere->npoints  - sh->nburied, 
	   atom_spec(box, i) );
    }
  } /* for i */
  return;
} /* purge_ngbs */

static void isolated_atom(shell_p sh, box_p box) {
  int i;

  i=sh - box->shells;

  box->atoms[i].pointer=sh;

  sh->flags &= ~BURIED;
  sh->flags |= ISOLATED;
  sh->nburied=0;
  sh->area= M_4PI;			/* area as if radius were 1.0 */

  warn("atom %d lies completely isolated\n", i);
  return;
} /* isolated_atom */

static void calc_acc(box_p box) { 
  /*
   * does    : loops over all shells in all cells of the box, and buries
   *   all occluded points. See also cal_acc2, below; could be faster than
   *   this routine
   *
   * gets    : a box with atoms
   *
   * affects : box[*]->pointflags[*]
   *
   * returns : number of atoms buried
   *
   * warns   : 
   *
   * comment : this is the place where shell_p->pointflags are allocated
   *   and initialized
   * 
   */

  int i, j, ix, iy, iz, nx, ny, nz, jx, jy, jz, ncells, offseta, offsetb
    , nngbs, nexposed, natomsexposed, k, size, npoints, natoms;
  shell_p *cells, a, b;
  sphere_p sphere;
  point_p Points;
  atom_p atom, atoms;
  float *xa, *xb, *coordinates, *xp;
  double d2, ra, rb, d, area;
  shell_p ngblist[ MAX_NNGBS ];
  float axes[ 4* MAX_NNGBS ], *lastaxis, *axis;
  SHORT pointflags[ MAX_TMPLT_PNTS ], *pf;

  natoms=box->natoms;
  atoms=box->atoms;

  ncells=box->ncells;
  cells=box->cells;

  nx=box->dims[0];
  ny=box->dims[1];
  nz=box->dims[2];

#define OFFSET(X,Y,Z) nz*(ny*X + Y)+Z

  natomsexposed=0;
  /* loop over all atoms: */
  offseta=0;
  for(ix=0; ix<nx; ix++) {
    for(iy=0; iy<ny; iy++) {
      for(iz=0; iz<nz; iz++, offseta++) {
	assert(offseta == OFFSET(ix, iy, iz));

	for(a = cells[offseta] ; a!=NULL; a=a->next.next) {
	  atom=a->atom;
	  xa=atom->xyz;
	  ra=xa[3];
	  
	  /* first assumption: pretend atom is buried */
	  atom=a->atom;
	  atom->pointer=NULL;
	  a->flags |= BURIED;
	  a->nburied=a->sphere->npoints;

	  /* find neighbours: */
	  axis=axes;
	  nngbs=0;
	  /* walk all cells adjacent to that of a */
	  for ( jx= (ix ? ix-1:0); (jx<=ix+1 && jx<nx); jx++) {
	    for ( jy= (iy ? iy-1:0); jy<=iy+1 && jy<ny; jy++) {
	      for ( jz= (iz ? iz-1:0); jz<=iz+1 && jz<nz; jz++) {
		offsetb=OFFSET(jx, jy, jz);
		assert(offsetb<ncells);
		/* walk all its contents: */
		for(b = cells[offsetb] ; b!=NULL; b=b->next.next) {
		  if (a==b)		/* same atom: skip */
		    continue;
		  xb=b->atom->xyz;
		  rb=xb[3];
		  d2=0.0;
		  for (k=0; k<3; k++) {
		    d=xb[k] - xa[k];
		    axis[k]=d;		/* store already now ! */
		    d2 += SQR(d);
		  }		    
		  if ( d2 < SQR(ra + rb) ) { /* a neighbour */
		    if ( d2 < SQR(ra-rb) ) { /* one buried other completely */
		      if (ra >= rb) {	/* a buries b: b can't be neighbour */
			continue;	/* next b */
		      } else {		/* b buries a: skip to next a */
/*			warn("atom %d buries atom %d completely\n",
 *			     b->atom - atoms, a->atom - atoms);
 */			
			if(a->pointflags)
			  FREE(a->pointflags);
			if(a->ngbs)
			  FREE(a->ngbs);
			goto NEXT_A;
		      }
		    }
		    axis[3]=((d2 - SQR(rb))/ra + ra)*0.5;
		    axis+=4;
		    ngblist[nngbs++]=b;
		    assert(nngbs <= MAX_NNGBS);
		  } /* if neibor */
		} /* for b */
	      }	/* for jz */
	    } /* for jy */
	  } /* for jx */

	  if(nngbs==0) {		/* isolated atom */
	    natomsexposed++;
	    isolated_atom(a, box);
	    goto NEXT_A;
	  }

	  /* loop over all points of atom a: */
	  area=0.0;
	  sphere=a->sphere;
	  coordinates=sphere->coordinates;
	  npoints=sphere->npoints;
	  size=npoints*sizeof(SHORT);
	  memset(pointflags, 1, size);
	  Points=sphere->Points;

	  nexposed=0;
	  lastaxis=axes;
	  for (i=npoints-1, xp=coordinates+3*i;
	       i>=0; 
	       i--, xp -=3) {
	    if ( DOT3P(xp, lastaxis ) < lastaxis[3] ) {  /* point not buried */
	      for (j=nngbs-1, axis=axes+4*j ;
		   j>=0; 
		   j--, axis -=4) {	/* walk whole list */
		if (DOT3P(xp, axis) > axis[3]) { /* point is buried */
		  lastaxis=axis;
		  goto NEXT_I;
		}
	      }	/* for j */
	      /* point is apparently exposed: */
	      area += Points[i].area ;
	      pointflags[i]=0;		/* unset the BURIED flag */
	      nexposed++;
	    } /* if dot < lastdot */
	  NEXT_I:			/* or after loop over atoms: */
	    ;
	  } /* for i */

	  if(nexposed) {		/* SURPRISE: atom not buried ! */
	    natomsexposed++;
	    if (nexposed==npoints) {
	      isolated_atom(a, box);
	      goto NEXT_A;
	    }	    

	    a->flags &= ~BURIED;	/* reset BURIED flag */
	    a->area=area;		/* area as if radius were 1.0 ! */
	    a->nburied= npoints - nexposed ;
	    atom->pointer=a;

	    /* ONLY NOW ALLOCATE POINTFLAGS AND NGBS:   */
	    /* if (keep_pointflags) ...  */
	    a->pointflags=(SHORT*)AREALLOC(a->pointflags, size );
					/* since nr of points may have */
					/* changed in mean time ! */
	    pf=a->pointflags+(npoints-1); /* end of mem, for counting down */
	    for (i=npoints-1; i>=0; i--) /* copy found flags to allocated */
	      *pf-- = (pointflags[i]&BURIED); /* mem */

	    size=nngbs*sizeof(shell_p);
	    a->nngbs=nngbs;
	    
	    a->ngbs=(shell_p *)AREALLOC(a->ngbs, size);
	    memcpy(a->ngbs, ngblist, size);
	  } else {			/* atom is buried */
	    if(a->pointflags) { 
	      FREE(a->pointflags);
	      a->pointflags=NULL;	/* avoid future free()ings */
	    }
	    if(a->ngbs) {
	      FREE(a->ngbs);
	      a->ngbs=NULL;		/* avoid future free()ings */
	    }
	  } /* if nexposed */
	NEXT_A:
	  ;
	} /* for(a) */
      }	/* for(iz) */
    } /* for(iy) */
  } /* for(ix) */
#undef OFFSET				/* avoid trouble later on */

  box->nburied= box->natoms - natomsexposed ;
  return;
} /* calc_acc */

int accessability(box_p box) { 
  /*
   * does    : buries atoms, set BURIED flags on buried points, and
   *   assigns accessible areas to shells (not to atom->AREAS !)
   *
   * gets    : atoms
   *
   * affects : ((shell_p)atom[i].pointer)->points[i] is set BURIED if point
   *   is buried
   *
   * returns : nr of atoms buried
   *
   * warns   : not
   *
   * comment : 
   */

  gridify(box);
  calc_acc(box); 
  purge_ngbs(box);

  return box->nburied;
} /* accessability */

void free_accessability(box_p box) {
  /*
   * does    : frees shell[*]->ngbs and shell[*]->pointflags
   *
   * gets    : the usual suspects
   *
   * affects : shell->ngbs & shell->pointflags are set to NULL
   *
   */

  int i, natoms;
  shell_p sh;

  natoms=box->natoms;
  
  /* in fact, following has to be done by free_accessability() or so: */
  for (i=0, sh=box->shells; i<natoms; i++, sh++) { 
    if (sh->pointflags)	{		/* if not buried, or so ... */
      FREE(sh->pointflags);
      sh->pointflags=NULL;
    }
    if (sh->ngbs) {			/* if not buried, or so ... */ 
      FREE(sh->ngbs);
      sh->ngbs=NULL;
    }
  }
  return;
} /* free_accessability */
