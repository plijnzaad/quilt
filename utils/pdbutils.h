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


#include "pdb.h"

#define CHECK_ATOMS BIT(0)		/* consult atom->flags & CALC_AREA */
					/* if atom should be considered */

#define MAX_NNGBS 256			/* max nr. of neigbours to any atom */

int surface_areas(int nvertices, structure_p, uint flags);
  /* calculates atomic areas based on atom->RADIUS and puts them in */
  /* atom->AREA, using a fibonacci point distribution of more than npoints */
  /* vertices. If ! (flags&CHECK_ATOMS), all atoms are used; if */
  /* (flags&CHECK_ATOMS), only the atoms for which (atom->flags&CALC_AREA) */
  /* != 0 are considered */

int phi_psi(structure_p str);
  /* calculates usual phi-psi angles, and puts results (in degrees) in */
  /* residue->phi,psi */
