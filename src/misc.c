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

/* miscellaneous routines for quilt */

#include <math.h>
#include "extra-math.h" 

#include "utils.h"
#include "alloc.h"
#include "list.h"

#include "quilt.h"

/*** program global variables ***/
/* options are in main.args ! */
/* things i couldn't really make local ... : */

FILE * logfile;
int overall_npatches			/* set in rims_to_patches */
,  maxnpoints;				/* set in read_test_file */

atom_p root_of_atoms;			/* base of atoms; sometimes used */
shell_p root_of_shells;			/* base of shell_p allocation  */


SHORT *twohistories[2];			/* for use in mark_interior */

int ** fast_mod;			/* matrix for fast modulo-taking */

/* local prototype */
static double * set_limits(double density);

int ** fast_modulo_matrix(int from, int to) { 
  /* returns an array that can be used to speed up 'w%ngbs' operations. This */
  /* routine just returns a matrix that would look like 
   * array[3] = { 0, 1, 2,           0, 1, 2, ... }  
   * array[4] = { 0, 1, 2, 3,        0, 1, 2, 3 ... }  
   * and would be used as mod = array[ngbs][w] (reverse of the % notation !!)
   * rows are repeated FASTMOD_NTIMES times to the right
   */


  int ** arrp, *ip;
  int i, denom, k, l, n;
  
  die("fast_modulo_matrix is bugged ...");

  n = to - from +1;
  arrp = (int**)MALLOC(n*sizeof(int*));
  
  for (i=0; i <= n; i++) {
    denom=from+i;
    ip = (int*)MALLOC(denom * FASTMOD_NTIMES * sizeof(int));
    arrp[i] = ip;
    for (k=0; k<FASTMOD_NTIMES; k++) 
      for (l=0; l<denom; l++) 
	*ip++ = l; 
  }
    
  return arrp - from;
} /* fast_modulo_matrix */

#ifdef NOT_DEFINED
int fast_mod_func(int i, int j) {	/* for debugging */
  int m;

  assert(i>=0);
  assert(j==5 || j==6);
  assert( i < FASTMOD_NTIMES*j );

  m=fast_mod[j][i];
  assert(m == i%j);
  return m;
}
#endif

sphere_p choose_uni_sphere(int pts_per_atom) { 
  /* choose universal sphere type for all atoms, and set maxnpoints */
  int i;
  sphere_p t;

  for (i=0; i < N_TMPLT_SPHERES; i++) { 
    t=Template_spheres[ i ];
    if (t->npoints == pts_per_atom) { 
      maxnpoints = t->npoints;
      return t;
    }
  }
  BUG;					/* since should have been caught */
					/* during command line processing */
  return NULL;				/* make gcc -Wall shut up */
} /* choose_uni_sphere */

sphere_p choose_sphere(atom_p atomp, double density) {
  /* choose sphere type per atom. This has to be not too slow, so use */
  /* table (set up by set_limits). Also keep track of maximum sphere size */
  /* used so far. */

  int i=0;
  double r=atomp->RADIUS;
  sphere_p t=NULL;
  static double * limits;

  if (!limits) 
    limits= set_limits(density);

  if (r <= limits[0])			/* too small */
    t=Template_spheres[0];
  else if (r > limits[ N_TMPLT_SPHERES -2 ]) /* too big */
    t = Template_spheres[ N_TMPLT_SPHERES -1 ];
  else {				/* find in table */
    while( r > limits[i] )
      i++;
    t = Template_spheres[ i ] ; 
  }
  assert(t);

  if (t->npoints > maxnpoints)
    maxnpoints = t->npoints;

  return t; 
} /* choose_sphere */

static double * set_limits(double density) {
  /* sets up table to do easy lookup of which template sphere to use, */
  /* given a radius. Also writes table to stderr */
  int i,n,l;
  static int halfway[ N_TMPLT_SPHERES -1 ] = TEMPLATE_NPOINTS_HALFWAY;
  double * rlimit; 
  double r, prevr, d1, d2; 


  rlimit = (double*)MALLOC( sizeof(double)*(N_TMPLT_SPHERES -1) );

  for (i=0; i < N_TMPLT_SPHERES -1 ; i++)  { 
    l = halfway[i]; 
    rlimit[i] = sqrt( (double)l / (M_4PI * density) );
  }

  /* print out table with radius limits per sphere type and corresponding */
  /* density interval: */
  n = Template_spheres[0]->npoints;
  r = rlimit[0];
  d1 = n/(r*r*M_4PI);

  warn("\
setting densities ...\n\
radius        nr  pts/(A^2)\n\
r <= %5.3g:  %3d  [%5.3g - 99.99]\n", r, n, d1);

  for (i=1; i < N_TMPLT_SPHERES -1 ; i++) {
    n = Template_spheres[i]->npoints;
    prevr=rlimit[i-1];
    r = rlimit[i];
    d1 = (double)n/(r*r*M_4PI);
    d2 = (double)n/(SQR(prevr)*M_4PI);
    warn("r <= %5.3g:  %3d  [%5.3g - %5.3g]\n", r, n, d1, d2);
  }

  /* Don't forget last one; not represented in table (r should be ok still) */
  n = Template_spheres[ N_TMPLT_SPHERES -1 ]->npoints;
  d1 = (double)n/(r*r*M_4PI);
  warn("r >  %5.3g:  %3d  [ 0.00 - %5.3g]\n", r, n, d1); 
  return rlimit;
} /* set_limits */

const char * atom_spec(box_p box, int i) { 
  static char spec[80];
  atom_p atom;

  atom=box->atoms+i;

  if (atom->residue)			/* give full atom name */
    sprintf(spec, "[ atom %d, ( %c %c%hd%c@%s ) ]",
	    i,
	    atom->residue->type, 
	    atom->residue->chain->name, 
	    atom->residue->nr, atom->residue->inscode,
	    atom->name);
  else
    sprintf(spec, "[ atom %d ]", i);

/*  sprintf(spec, ")"); */
  
  return &spec[0];
} /* atom_spec */



