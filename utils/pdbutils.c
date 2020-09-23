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

/* module to do common additional things with protein structures, such as */
/* surf. accessability calc. */
/* When using the accessible surface calculation, please refer to the */
/* following paper:
 *   Eisenhaber, F., Lijnzaad, P., Argos, P., Sander, C. & Scharf, M. (1995)
 *   J. Comp. Chem. _16_, 273-284. "The double cubic lattice method:
 *   Efficient approaches to numerical integration of surface area and
 *   volume, and to dot surface contouring of molecular assemblies"
 */


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vecmat.h" 
#include "extra-math.h" 
#include <ctype.h>
#include <string.h> 

#include "vecmat.h"			/* macros for vector things */
#include "utils.h"			/* some basic definitions/prototypes */
#include "list.h"
#include "alloc.h"
#include "hash.h"

#include "pdb.h"			/* the data structure etc. */
#include "pdbutils.h"

/* file global variables & data  (none) */


#define MINDIST 0.0001                         /* no atoms should be closer than this */

#define ISOLATED(x)  warn("%s lies isolated\n", full_atom_name(x)),\
				(x)->AREA=M_4PI * SQR((x)->RADIUS)

static int do_access(int nverts, float * vertices,
		     atom_p *cells, int * boxdims, 
                     atom_p start_of_atoms /* for reporting only */ ) { 
  /* calculate surface accessability and put areas into atom->area; return */
  /* nr of buried atoms */ 
  int i, j, ix, iy, iz, nx, ny, nz, jx, jy, jz, offseta, offsetb
    , nngbs, nexposed, nburied, k, iarea, nei[3];
  atom_p a, b, ngblist[ MAX_NNGBS ];
  float *xa, *xb;			/* atom coordinates */
  float *xp;				/* vertex coordinates */
  float axes[ 4* MAX_NNGBS ], *lastaxis, *axis; /* @ */
  double areapp, d2, ra, rb, d;
  static int neibors [ 27*3 ] = { 
    0,  0,  0, 
    0,  0, -1, 
    0,  0,  1, 
    0, -1,  0, 
    0,  1,  0, 
   -1,  0,  0,      
    1,  0,  0,      
    	  
    0, -1, -1,     
    0, -1,  1,      
    0,  1, -1,     
    0,  1,  1,      
   -1,  0, -1,     
   -1,  0,  1,      
    1,  0, -1,     
    1,  0,  1,      
   -1, -1,  0,      
   -1,  1,  0,      
    1, -1,  0,      
    1,  1,  0,      
    	  
   -1, -1, -1,
   -1, -1,  1,      
   -1,  1, -1,     
   -1,  1,  1,      
    1, -1, -1,     
    1, -1,  1,      
    1,  1, -1,     
    1,  1,  1,      
  };
  
  nx=boxdims[0];
  ny=boxdims[1];
  nz=boxdims[2];
  
#define OFFSET(X,Y,Z) nz*(ny*(X) + (Y))+(Z)
  
  areapp = M_4PI/nverts;
  
  nexposed=nburied=0; 
  nei[0]= 0;				/* directions to look in for finding */
  nei[1]= -1;				/* neighbour cells. */
  nei[2]= 1;


/* loop over all atoms: */
  
  offseta=0;
  for(ix=nx-2; ix>=1; ix--) {		/* unit offset, to allow a row of */
    for(iy=ny-2; iy>=1; iy--) {		/* empty cells around the box */
      for(iz=nz-2; iz>=1;  iz--) {
	offseta= OFFSET(ix, iy, iz);
	for(a = cells[offseta] ; a!=NULL; a=a->pointer) {
	  xa=a->xyz;
	  ra=xa[3];
	  
	  /* find all neighbours of a: */
	  axis=axes;
	  nngbs=0;
	  /* walk all cells adjacent to that of a */
#define OFFSET_J(X, Y, Z) (OFFSET( ix + nei[X], iy + nei[Y], iz + nei[Z]))
	  for (i=26; i>=0; i-- ) {
	    jz = iz + neibors[ 3*i+2 ];
	    jy = iy + neibors[ 3*i+1 ];
	    jx = ix + neibors[ 3*i ];
	    offsetb=OFFSET(jx, jy, jz);

/* 	  for ( jx = 2; jx >= 0; jx-- ) {
 * 	    for ( jy = 2; jy >= 0; jy-- ) {
 * 	      for ( jz = 2; jz >= 0; jz-- ) {
 *	  offsetb=OFFSET_J(jx, jy, jz);
 */

/* 		printf("%d %d %d - %d %d %d : %d\n",
 * 		       jx, jy, jz, ix+nei[jx], iy+nei[jy],iz+nei[jz], offsetb);
 */		/* the jx, jy, and jz may and up in the empty margins! */
		/*		assert(offsetb<ncells); */
		/* walk all its contents: */

	    for(b = cells[offsetb] ; b!=NULL; b=b->pointer) {
	      if (a==b)			/* same atom: skip */
		continue;
	      xb=b->xyz;
	      rb=xb[3];
	      d2=0.0;
	      for (k=2; k>=0; k--) {
		d=xb[k] - xa[k];
		axis[k]=d;		/* store already */
		d2 += SQR(d);
	      }		    
              /* sanity check here: no distance should be smaller than
               * MINDIST */
              if (d2 < MINDIST) 
                DIE("(squared)distance between atoms %d and %d is %7.3g,",
                    a-start_of_atoms, b-start_of_atoms, d2);

	      if ( d2 < SQR(ra + rb) ) { /* a neighbour */
		if ( d2 < SQR(ra-rb) ) { /* one buried other completely */
		  if (ra >= rb)		/* a buries b: b can't be neighbour */
		    continue;		/* next b */
		  else {		/* b buries a: skip to next a */
		    nburied++;
		    goto NEXT_A;
		  }
		}
		axis[3]=((d2 - SQR(rb))/ra + ra)*0.5; /* scaled dotprod. */
		axis+=4;
		ngblist[nngbs++]=b;
	      } /* if neibor */
	    } /* for b */
	  } /* for i */
/* 	    } /# for jy #/
 * 	  } /# for jx #/
 */
	  if (nngbs > MAX_NNGBS)	/* hopefully didn't crash earlier */
	    die("\
number of neighbouring atoms too large in do_acces(); [at %s]\n\
reduce atom and/or solvent radii, or change MAX_NNGBS and recompile\n");

	  if(nngbs==0) {		/* isolated atom */
	    ISOLATED(a);
	    nexposed++;
	    goto NEXT_A;
	  }

	  iarea=0;			/* nr of accessible dots */

	  /* Now here comes Michael Scharf's tremendous trick: walk first */
	  /* over points; if one atoms buries a point, the neigbour point */
	  /* often is also buried by the same atom. Tightest part of the */
	  /* algorithm. Also, dotproduct and counting down slightly faster. */ 

	  lastaxis=axes;
	  for (i=nverts-1, xp=vertices+3*i; /* iterate over points */
	       i>=0; 
	       i--, xp -=3) {
	    if ( DOT3P(xp, lastaxis ) < lastaxis[3] ) {	/* point not buried */
	      for (j=nngbs-1, axis=axes+4*j ; /* iterate over neibor atoms */
		   j>=0; 
		   j--, axis -=4) {
		if (DOT3P(xp, axis) > axis[3]) { /* point is buried */
		  lastaxis=axis;
		  goto NEXT_I;
		}
	      }	/* for j */
	      /* point is apparently exposed: */
	      iarea++; 

/* 	      { 
 * 		float *x=a->xyz, r=x[3];
 * 		fprintf(stdout, ".dot %7.3f %7.3f %7.3f\n", 
 * 		       x[0]+r*xp[0], x[1]+r*xp[1], x[2]+r*xp[2]);
 * 	      }
 */
	    } /* if dot < lastdot */
	  NEXT_I:			/* or after loop over atoms: */
	    ;
	  } /* for i */

	  if(iarea) {
	    nexposed++;
	    if (iarea == nverts) {
	      ISOLATED(a);
	      goto NEXT_A;
	    }	    
	    a->AREA= iarea* areapp * SQR(a->RADIUS); /* that's all folks */
	  } else { 
	    nburied++; 
	    a->AREA= 0.0;
	  }
	NEXT_A:
	  ;
	} /* for(a) */
      }	/* for(iz) */
    } /* for(iy) */
  } /* for(ix) */
#undef OFFSET				/* avoid trouble later on */
  /*  return nexposed; */
  return nburied;
} /* calc_acc */

#define VERT_SORT_DIM 4			/* unit sphere grid dimension */
#define VERT_SORT_SPACING (2.0/VERT_SORT_DIM)

static int vertex_cmp(const void * a, const void * b) {
  int ax, ay, az, bx, by, bz, offseta, offsetb;

#define INTDIM(X,I) ((int)floor((((float*)(X))[(I)] + 1.0)/VERT_SORT_SPACING))
#define OFFSET(DIM, X,Y,Z) DIM*(DIM*X + Y)+Z

  ax=INTDIM(a, 0);
  ay=INTDIM(a, 1);
  az=INTDIM(a, 2);

  bx=INTDIM(b, 0);
  by=INTDIM(b, 1);
  bz=INTDIM(b, 2);

  offseta=OFFSET(VERT_SORT_DIM, ax,ay,az);
  offsetb=OFFSET(VERT_SORT_DIM, bx,by,bz);

  return offseta - offsetb;
} /* int vertex_cmp */
#undef OFFSET

static int fibo_vertices(int a, int b, float*vertices) { 
  /* writes b points in to 1dim array vertices */
  int i,j;
  double s,f,r,phi,z;
  
  s=1.0/b;				/* factor scaling i to 0.0 - 1.0 */
  f=s*2*M_PI;				/* factor scaling j to 0.0 - 2PI */
  j=0;
  for (i=0; i<b; i++) {
    /* 'wrap' i around a circle of circumference b, with step size a; */
    /* since a and b do not have a common factor (property of fibonacci */
    /* numbers, by induction), this will in exactly b steps bring you back */
    /* to b (or rather b%b==0), without selecting any intermediate value */
    /* twice. When this series is translated to points on the */
    /* circumference of a circle, every consecutive point has the same */
    /* distance to its predecessor, approaching the Golden Ratio times the */
    /* circles' circumference for large b. */

    j=(i*a)%b;
    phi=f*j;				/* non-regular jumping around circle */
    z=1.0 - 2*i*s;			/* z-value */
    r=sqrt((1-SQR(z)));			/* radius projected to xy-plane */
    vertices[3*i]=   r*cos(phi);
    vertices[3*i+1]= r*sin(phi);
    vertices[3*i+2]= z;
  }
  return 0;
} /* fibo_vertices */

static float * gen_vertices(int *npoints) {
  int a, b, c, n;
  float *vertices;

  if ( 0 ) { 
    int i=0;
    FILE *file=fopen("points", "r");
    char line[132];
    
    if ( fscanf(file, "%d\n", &n) != 1)
      die("no nr of points\n");
    vertices=(float*)MALLOC(n*sizeof(float[3]));
    while( fgets(line, 132, file) ) { 
      sscanf(line, "%f %f %f", vertices+3*i,vertices+3*i+1,vertices+3*i+2);
      i++;
    }
    assert(i == n+1);
    fclose(file);
  } else {
    n= *npoints;
    a=0; b=1;				/* Lucas numbers: a=0, b=2,3,4,etc.  */
    do {
      c=a+b;
      a=b;
      b=c;
    } while (b<n);
    n=b;
    
    vertices=(float*)MALLOC(n*sizeof(float[3]));
    
    fibo_vertices(a, b, vertices);
    qsort(vertices, n, sizeof(float[3]), vertex_cmp);
  }

  *npoints=n;
  return vertices;
} /* gen_vertices */

static int box_dimensions(structure_p str, uint flags,
			  int *boxdims, double *limits) { 
  /*
   * does    : establish dimensions of protein and hence of box (grid) to
   *   be used (depending on flags; not documented)
   *
   * gets    : protein STR
   *
   * affects : BOXDIMS[0,1,2] is set to x,y,z integer dimensions of box;
   *   LIMITS is set to xmin,xmax,ymin,ymax,zmin,zmax,Rmin,Rmax of protein
   *
   * returns : no of cells to be used for grid
   *
   * warns   : 
   *
   * comment : 
   */
  int i, j, n, k, ncells;
  atom_p atom;
  double f, spacing, maxrad;		/* xmin, xmax, etc., radmin, radmax. */

  n=str->natoms; atom=str->atoms; k=0;

  if (flags & CHECK_ATOMS) {		/* check which atoms are needed */
    maxrad=0.0;
    for (i=0; i<3; i++) {
      limits[ 2*i ] = HUGE_VAL;
      limits[2*i+1] = -HUGE_VAL;
    }
    for (i=0; i<n; i++, atom++) {
      if_not (atom->flags & CALC_AREA)
	continue;
      k++;

      for (j=0; j<3; j++) {		/* find extrema of coordinates */
	f=atom->xyz[j];
	if (f < limits[2*j])		/* minima */
	  limits[2*j]= f;
	if ( f > limits[2*j+1])		/* maxima */
	limits[2*j+1]=f;
      }
      if ( atom->RADIUS > maxrad ) maxrad=atom->RADIUS;
    } /* for i< n */
  } else {				/* use all atoms */
    for (j=0; j<6; j++) 		/* extrema of coordinates are known */
      limits[j]=str->bounding_box[j];
    maxrad=str->maxradius;
    assert(maxrad > 0.001);		/* if not: set_radii() forgotten ? */
  }
  
  limits[7] = maxrad;
  spacing = 2.0 * maxrad;

  ncells=1;				/* establish nr of cells needed */
  for (i=0; i<3; i++) { 
    f=(limits[2*i +1] - limits[2*i])/spacing;
    n=(int)floor(f) + 1 + 2 ;
    boxdims[i]=n;
    ncells *=n;
  }
  return ncells;
} /* box_dimensions */

static int box_atoms(int natoms, atom_p atoms, uint flags,
		     int *boxdims, double *limits,
		     atom_p * cells, void **pointers) { 
  /* puts atoms into cells, and saves atom->pointers to pointers */
  int i, j, offset, n;
  atom_p atom=atoms;
  double spacing = 2.0*limits[7];	/* twice maximal radius */
  
  for (i=0; i<natoms; i++,atom++) {
    if( (flags & CHECK_ATOMS) && !(atom->flags & CALC_AREA))
      continue;

    offset=0;
    for (j=0; j<3; j++)  {
      n=(int)floor( (atom->xyz[j] - limits[2*j])/spacing) + 1; /* (margin) */
      offset = offset*boxdims[j] + n;
    }
    /* assert(offset < ncells); */

    pointers[i]=atom->pointer;		/* save old pointer */
    atom->pointer=cells[offset];	/* add to linked list */
    cells[offset]=atom;
  } /* for i < natoms */

  return 0;
} /* box_atoms */

static void debox_atoms(structure_p str, uint flags, void ** pointers) {
  /* restores atom->pointers */
  int i, natoms=str->natoms;
  atom_p atom=str->atoms;
  
  for (i=0; i<natoms; i++, atom++) { 
    if ((flags & CHECK_ATOMS) & !(atom->flags & CALC_AREA) )
      continue;
    atom->pointer=pointers[i];
  }
  return;
} /* debox_atoms */

int surface_areas(int npoints, structure_p str, uint flags) { 
  /* please refer to the following paper when using this routine:
   *   Eisenhaber, F., Lijnzaad, P., Argos, P., Sander, C. & Scharf, M. (1995)
   *   J. Comp. Chem. _16_, 273-284. "The double cubic lattice method:
   *   Eficient approaches to numerical integration of surface area and
   *   volume, and to dot surface contouring of molecular assemblies"
   */
  /* calculates atomic areas based on atom->RADIUS and puts them in */
  /* atom->AREA, using a fibonacci point distribution of more than npoints */
  /* vertices. If ! (flags&CHECK_ATOMS), all atoms are used; if */
  /* (flags&CHECK_ATOMS), only the atoms for which atom->flags & CALC_AREA */
  /* are considered */

  int boxdims[3], size, ncells, natoms;
  atom_p *cells;
  void **pointers;
  double limits[8];
  static int nverts;
  static float *vertices;		/* @ */

  if(str->maxradius <= 0.001)
    die("Bug at %s:%d: maxradius <= 0.001 in surface_areas(); set_radii() forgotten?\n",
	__FILE__, __LINE__);

  natoms=str->natoms;
  if (! vertices) {			/* get point distribution */
    nverts=npoints;			/* approximate nr of points */
    vertices=gen_vertices(&nverts);
    warn("using %d points per sphere for area calculation\n", nverts);
  }

  ncells=box_dimensions(str, flags, boxdims, limits); /* determine boxing params */
  size=ncells*sizeof(atom_p);
  cells=(atom_p*)ALLOCA(size);		/* allocate on stack */
  bzero((char*)cells, size);		/* important */

  pointers=(void**)ALLOCA(natoms*sizeof(void*)); /* to save pointers */

  box_atoms(natoms, str->atoms, flags, boxdims,	/* do the boxing */
	    limits, cells, pointers);
/*   { int i;
 *   for (i=0; i<ncells; i++) 
 *     if (cells[i])
 *       warn("cell %d length %d\n", i, Llength((list_p)cells[i]));
 *   }
 */  
  do_access(nverts, vertices, cells, boxdims, str->atoms);
  debox_atoms(str, flags, pointers);		/* restore pointers */

  AFREEA(pointers);
  return 0;
} /* surface_areas */

int phi_psi(structure_p str) {
  /* calculates usual phi-psi angles, and puts results in residue->phi,psi */
  int i, j, nres;
  residue_p res;
  chain_p ch;
  double phi, psi;
  float *n, *ca, *c;

  for (i=0; i<str->nchains; i++)  { 
    ch=str->chains+i;
    res=ch->residues;
    nres=ch->nresidues;
    for (j=0,res=ch->residues; j<nres; j++, res++) {
      if(res->flags.bb_missing) {
	res->phi = res->psi = -360.0;
	continue;
      }
      n= &res[0].atoms[0].xyz[0];
      ca= &res[0].atoms[1].xyz[0];  
      c= &res[0].atoms[2].xyz[0];
      if (j!=0)				/* not for first of chain */
	phi=fdihedral(res[-1].atoms[2].xyz, n, ca, c); /* c(-1), ... */
      else
	phi= M_2PI;

      if (j != nres-1 )			/* not for last of chain */
	psi= fdihedral(n, ca, c, res[1].atoms[0].xyz); /* ..., n(+1) */
      else
	psi= M_2PI;

      res->phi = DEGPERRAD*phi;
      res->psi = DEGPERRAD*psi;

    } /* for j */
  } /* for i */
  return 0;
} /* phi_psi */
