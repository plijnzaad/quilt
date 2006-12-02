/* ------------------------------------------------------------------ *
 *								      *
 * Permission to copy all or part of this work is granted, provided   *
 *  that this NO WARRANTY and this copyright notice are retained      *
 *  verbatim and are displayed conspicuously.  If anyone needs other  *
 *  permissions that aren't covered by the above, please contact the  *
 *  author.							      *
 *								      *
 * NO WARRANTY: THIS WORK IS PROVIDED ON AN "AS IS" BASIS.  THE       *
 *  AUTHOR PROVIDES NO WARRANTY WHATSOEVER, EITHER EXPRESS OR         *
 *  IMPLIED, REGARDING THE WORK, INCLUDING WARRANTIES WITH THE WORK,  *
 *  INCLUDING WARRANTIES WITH RESPECT TO ITS MERCHANTABILITY OR       *
 *  FITNESS FOR ANY PARTICULAR PURPOSE.				      *
 *								      *
 * Philip Lijnzaad, lijnzaad@embl-heidelberg.de			      *
 *								      *
 * ------------------------------------------------------------------ */

/* the main routines for finding patches on a molecular surface */

#include <stdio.h>
/* #define _BSD 1				 for proper RAND_MAX !!!!  */
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "extra-math.h"

#include "getargs.h"
#include "utils.h"
#include "alloc.h" 
/* #include "files.h" */

#include "pdb.h"
#include "quilt.h"

#define ATOM_SURFACC_THRESHOLD 10.0	/* redistrib atom if more
					   accessible than this (Angstroms) */
#define ATOM_SURFACC_THRESHOLD2 5.0	/* for use in randomize2 */

#define LOW_A  42.0			/* yes, LOW_A is higher than HIGH_A */
#define HIGH_A 38.0
#define LOW_B  0.7
#define HIGH_B 0.8

#define MID_A 37.167
#define MID_B 0.753

/* following is not implemented very nicely, i admit */
#define PROBE 1.4

#define Nradius (1.7 + PROBE)
#define Oradius (1.4 + PROBE)
#define Cradius (1.9 + PROBE)
#define Sradius (1.8 + PROBE)

#define Ncode 0
#define Ocode 1
#define Ccode 2
#define CScode 2			/* cheating: take them together */
#define Scode 3


static int box_acc(box_p box) {
  /* quick hack: reading pdb file with last column accessability; get totals */
  int i;
  char c;
  atom_p at;
  double area, totarea, Carea, Sarea, Narea, Oarea;

  totarea = Carea = Sarea = Oarea = Narea=0.0;
  for (i=0; i<box->natoms; i++) { 
    at = &box->atoms[i];
    area=at->AREA;
    c=at->name[1];
    switch(c) { 
    case 'C': Carea += area; break;
    case 'S': Sarea += area; break;
    case 'N': Narea += area; break;
    case 'O': Oarea += area; break;
    default: break;			/* or category X ? */
    }
    totarea  += area;
  } /* for i < natoms */
  box->totarea=totarea;
  box->Narea=Narea;
  box->Oarea=Oarea;
  box->Carea=Carea;
  box->Sarea=Sarea;
  return 0;
  /*  return check_area(box); */
} /* box_acc */

static int redistrib(char type, double radius, double *target_fractionp,
		     double * new_areap, int index, atom_p mutable[] ) { 
  /*
   * does    : randomly redistributes TOTAREA over atoms in MUTABLE,
   *   starting at INDEX and turning them into type TYPE with radius
   *   RADIUS, until *TARGET_FRACTIONP is attained
   *
   * gets    : 
   *
   * affects : much of mutable[*]->{type, radius, area}, *target_fractionp
   *
   * returns : number of atoms actually mutated, to be used as index
   *
   * warns   : when there's not enough to be mutated
   *
   * comment : end of mutable-array is signalled by (void*)-1
   */
  int i;
  atom_p at;
  double new_area, sofar, area1, areaR, target_fraction, frac, sofar2
    , new_area2;

  new_area = *new_areap;
  target_fraction = *target_fractionp; /* the one to be reached */

  i=index;				/* where to start */
  sofar =0.0;
  while( 1  ) {
    at=mutable[i];
    if (at == (void*)-1) 		/* sentinel: end of array */
      return warn("***** mutable atoms  (%d) exhausted\n", i);

    area1=at->AREA;			/* with radius == 1.0 */
    areaR = area1 * SQR(radius);	/* with radius R */

    /* now find out what the new fraction would be; if too large, exit loop */
    sofar2 = sofar + areaR;
    new_area2 = new_area -  area1 * SQR(Cradius)  + areaR;
    frac=sofar2/new_area2;
    if ( frac > target_fraction)	/* exceeding the target: leave */
      break;

    sofar = sofar2;
    new_area = new_area2;
    at->name[1]=type;			/* 'N' or 'O' */
    i++;
  } /* while */

  *new_areap = new_area;
  *target_fractionp = sofar/new_area;

  return i;
} /* redistrib */

/* #include <sys/types.h> */
/* #include <time.h> */
/* #include <unistd.h> */
extern int time(void *);
extern int getpid(void);


double radii[4] = { Nradius, Oradius, Cradius, Sradius };

int adjust_hfobicity(double percentage, box_p box, 
		     double *Nfracp, double *Ofracp) {
  double ratio, hfyl, frac, Nfrac, Ofrac;
  
  frac= 1.0 - percentage;
  hfyl= (box->Oarea + box->Narea)/box->totarea;
  ratio= frac/hfyl;
  Nfrac = *Nfracp * ratio;
  Ofrac = *Ofracp * ratio;
  warn("adjusting Nfrac  %.4f -> $.4f, Ofrac %.4f %.4f\n", 
       *Nfracp, Nfrac, *Ofracp, Ofrac);
  *Nfracp = Nfrac;
  *Ofracp = Ofrac;
  return 0;
} /* adjust_hfobicity */

#define ERR_CORR_SLOPE (1.16)
#define ERR_CORR_INTERCEPT (-0.0673 )
int compensate(box_p box, double * Nfracp, double * Ofracp) { 
  /* try and correct systematic error in hfobicity thing */
  double hfob, new_hfob;
  
  hfob = (box->Carea + box->Sarea) / box->totarea;
  new_hfob = ERR_CORR_SLOPE * hfob + ERR_CORR_INTERCEPT;
  adjust_hfobicity(new_hfob, box, Nfracp, Ofracp);
  return 0;
}

static int seed_it(int seed) { 
  if (seed)
    srandom( seed );			/* gives same random sequence */
  else
    srandom( (int)time(NULL) * (getpid() % 19) ); /* seed random generator */
  return 0;
}

int ran_atoms(box_p box, int percentage) {  
  /*
   * does    : randomly changes atoms from C->{N,O} and back. If
   *   PERCENTAGE >= 0, the new hydrophobicity will be PERCENTAGE/100.0;
   *   If percentage  < 0, the proteins hydrophobicity is kept constant
   *
   * gets    : a protein, the atoms of which have the actual area in at->AREA.
   *   The box->*area's etc. should be set
   *
   * affects : at->AREA and at->name
   *
   * returns : 0
   *
   * warns   : about numers of atoms changed, etc.
   *
   * comment : not very nicely implemented
   */
  /* exspects areas in pdb file in last field !!!! */
  int i, nC, nS, nO, nN, n, index;
  atom_p at, *mutable;
  char type;
  double area, Narea, Oarea, Nignore, Oignore, new_area
    , CSignore, CSarea, Nfrac, Ofrac, x;

  box_acc(box);
  Narea=box->Narea;			/* area to be redistributed */
  Oarea=box->Oarea;
  CSarea=box->Carea +box->Sarea;
  Nfrac=Narea/box->totarea;
  Ofrac=Oarea/box->totarea;

  new_area=Nignore=Oignore=CSignore=0.0;

  mutable=(atom_p*)ALLOCA(box->natoms*sizeof(atom_p)); /* temp array */

  nC=nS=nO=nN=0;
  n=0;					/* nr of mutable atoms */
  for (at=box->atoms,i=0; i<box->natoms; at++,i++) {
    area= at->AREA;			/* read from file */
    type=at->name[1];
    if (area > ATOM_SURFACC_THRESHOLD ) { /* square Angstroms */
      switch(type) { 
      case 'N': at->AREA = area/SQR(Nradius); nN++; break; /* make 1 A radii */
      case 'O': at->AREA = area/SQR(Oradius); nO++; break;
      case 'C': at->AREA = area/SQR(Cradius); nC++; break;
      case 'S': at->AREA = area/SQR(Sradius); nS++; break;
      default: break;
      }
      at->name[1] = 'c';		/* means: can be changed */
      mutable[ n++] = at;
      new_area += SQR(Cradius)*at->AREA; /* area as if all atoms are 'C' */
    } else { 
      new_area += area;
    }
  } /* for i */

  mutable[n] = (void*) -1;		/* sentinel */

  seed_it(opt_general);

  randomize(n, sizeof(atom_p), mutable); /* simple shuffle of all pointers */

  warn("%d N, %d O, %d C, %d S, total %d *surface* atoms\n", 
       nN, nO, nC, nS, nN + nO + nC + nS);

  if (percentage > 100)                 /* this is an undocumented feature... */
    compensate(box, &Nfrac, &Ofrac);
  else if (percentage >= 0)
    adjust_hfobicity((double)percentage/100.0, box, &Nfrac, &Ofrac);

  index=0;
  x=Nfrac;
  index=redistrib('N', Nradius, &Nfrac, &new_area, index, mutable);
  warn("redistributed %d N atoms, area fraction %.12f -> %.12f\n",
       index, x, Nfrac);

  n=index;
  x=Ofrac;
  index=redistrib('O', Oradius, &Ofrac, &new_area, index, mutable);
  warn("redistributed %d O atoms, area fraction %.12f -> %.12f\n",
       index-n, x, Ofrac);

  for (at=box->atoms,i=0; i<box->natoms; at++,i++)
    at->name[1]=toupper(at->name[1]);

  set_radii( box->structure, AtomTypes, 0.0 );

  AFREEA(mutable);
  return 0;
} /* ran_atoms */

