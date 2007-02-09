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
#include "main.args"			/* command line arguments, options */


#define PATH ".:..:/junk/lijnzaad"	/* where to look for files */
#define BRKFLAGS (KEEP_COFACTORS  /* |  KEEP_IONS */) 
#define READFLAGS ( SET_RADII | BRKFLAGS | QUICKNDIRTY /* READ_DSSP */ )

/* global variables: all defined in misc.c */

#define AREA xyz[4]

static void pre_process_atoms(box_p box,
			      double extension,
			      double threshold,
			      double polar_extension,
			      int pts_per_atom,
			      double pt_density ) {
  /* allocate atomtr etc. for atoms, set radius extensions, and choose */
  /* template sphere to be used. Returns maximum radius used */

  /*
   * does    : allocates a box; allocates shell_p for atoms, sets
   *   atom[*]->pointer accordingly; saves original radii; extends radii
   *   according to arguments; chooses a sphere type to be used, and
   *   allocates shell->pointflags; determines min, max of x,y,z
   *
   * gets    : 
   *
   * affects : 
   *
   * returns : box
   *
   * warns   : about min, max radii
   *
   * comment : 
   */
  
  int i, j, size, natoms, extend_polar;
  shell_p shell;
  sphere_p uni_sphere = NULL;
  atom_p atom;
  float f;
  static float limits[8];			/* xyzr min,max */
  float *radii; 

  natoms=box->natoms;

  if (pts_per_atom) 			/*  use one template sphere of this */
    uni_sphere = choose_uni_sphere(pts_per_atom); /*  size for all atoms */

  assert(! box->shells);
  box->shells = (shell_p)CALLOC(box->natoms, sizeof(shell_t));

  if (! box->original_radii) {		/* first call: save them */
    set_radii(box->structure, AtomTypes, 0.0);
    radii = (float *)CALLOC(natoms, sizeof(float));
    for (i=0; i<natoms; i++)
      radii[i]=box->atoms[i].RADIUS;	/* save the original values */
    box->original_radii=radii;
  } else
    radii=box->original_radii;

#include <limits.h>
  for (i=0; i<4; i++) {			/* set extrema */
    limits[2*i]= HUGE_VAL;		/* so everything's less than this */
    limits[2*i+1]= -HUGE_VAL;		/* so everything's greater than this */
  }
  
  extend_polar = (polar_extension > 0.0 );
  
  for (i=0; i<natoms; i++) {
    atom=box->atoms+i;
    shell=box->shells+i;

    atom->RADIUS =  radii[i] + (float)extension;

    if( strchr("NO", atom->name[1]) ) {	/* a polar atom */
      shell->flags |= POLAR;		/* for convenience */
      if (extend_polar && atom->AREA > threshold)	/* square A's */
	atom->RADIUS += (float)polar_extension;	
    }

    for (j=0; j<4; j++) {		/* update extrema */
      f=atom->xyz[j];
      if (f < limits[2*j])		/* minima */
	limits[2*j]= f;
      if ( f > limits[2*j+1])		/* maxima */
	limits[2*j+1]= f;
    }

    shell->sphere = pts_per_atom ?	/* choose appropriate sphere */
      uni_sphere : choose_sphere(atom, pt_density);

    size=shell->sphere->npoints*sizeof(SHORT);
/*     if(!shell->pointflags)		/# 2nd time: don't allocate #/
 *       shell->pointflags=(SHORT*)MALLOC(size);
 * 
 *     bzero(shell->pointflags, size);
 */    
/*     shell->flags &= POLAR;		/# reset #/
 *     shell->nburied=0;
 */    
    shell->atom=atom;
    atom->pointer=shell;
  } /* for i */

  warn("radii after pre-processing between %5.3g and % 5.3g\n",
       limits[ MINRADIUS ], limits[ MAXRADIUS ]);
  memcpy(box->limits, limits, sizeof(float[8]));

} /* pre_process_atoms */

static void post_process_atoms(box_p box,
			       double reduction, double polar_reduction) {
  /* reduce atom radii */
  int i, natoms;
  atom_p atom, atoms;
  
  atoms=box->atoms;
  natoms=box->natoms;
  for (i=0, atom=atoms; i < natoms; i++, atom++) {
    atom->RADIUS -= (float)reduction;
    if(strchr("NO", atom->name[1]))
	    atom->RADIUS -= (float)polar_reduction;
  }
} /* post_process_atoms */

static void reset_atompointers(box_p box) {
  int i;

  for (i=0; i<box->natoms; i++) 
    box->atoms[i].pointer=NULL;
  return;
} /* reset_atompointers */


void free_triangulation(box_p box) { 
  /*
   * does    : frees all data associated with one triangulation
   *
   */

  free_connections(box);
  free_topology(box);
  free_accessability(box);
  FREE(box->shells);

  if(! box->old_shells)
    return;

  box->shells=box->old_shells;
  free_connections(box);
  free_topology(box);
  free_accessability(box);
  FREE(box->shells);
  
  return;
} /* free_triangulation */

static void free_rest(box_p box) {
  /*
   * does    : frees miscellaneous other stuff
   */

  FREE(box->cells);
  FREE(box->original_radii);


  if (box->structure) {
    free_pdb(box->structure);
    free_pdb_related();
  }
  else
    FREE(box->atoms);

  FREE(box);

  free_templates();

  return;
} /* free_rest */

static void do_patches(FILE *file, box_p box, char * nature) {
  /*
   * does    : calc.s all patches of chemical NATURE (eg. "CS", "NO", or
   *   "all", and writes to FILE
   *
   * gets    : BOX with patches
   *
   * affects : a lot. 
   *
   * returns : 
   *
   * warns   : 
   *
   * comment : mainly dispatches to get_surfaces(), which see. 
   */
  surf_p surf;

  if (!strcmp(opt_patches, "all")) {
    surf=get_surfaces(box, "");	/* surface of any kind: checks */
    print_surfaces(file, surf, "");
    if(dmw_on)
      free_surface(surf);		/*   surface components */

    surf=get_surfaces(box, "CS");		/* apolar surface */
    print_surfaces(file, surf, "CS");
    if(dmw_on)
      free_surface(surf);

    surf=get_surfaces(box, "NO");		/* polar surface */
    print_surfaces(file, surf, "NO");
    if(dmw_on)
      free_surface(surf);
  } else {
    surf=get_surfaces(box, opt_patches);
    print_surfaces(file, surf, opt_patches);
    if(dmw_on)
      free_surface(surf);		/*   surface components */
  }  
} /* do_patches */

static int transfer_acc(box_p box) {
  /* puts sh->area into atom->AREA; still UNIT SPHERE !  */
  int i;
  char c;
  atom_p at;
  double area, uarea, totarea, Carea, Sarea, Narea, Oarea;

  totarea = Carea = Sarea = Oarea = Narea=0.0;
  for (i=0; i<box->natoms; i++) { 
    uarea=box->shells[i].area;		/* unit sphere area */
    /* area *= SQR(box->atoms[i].RADIUS); */
    at = &box->atoms[i];
    at->AREA=uarea;
    area = SQR(at->RADIUS) * uarea;
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
} /* transfer_acc */

#ifndef catch_fpe
#define catch_fpe(X)			/* mask it */
#endif

#define BLURP "\n\nquilt version " VERSION "\n\nCopyright Philip Lijnzaad\n\n\n"

int main(int argc, char**argv) {
  int patchmode, n;
  box_p box;

  fprintf(stderr, BLURP);

  DOdmw(1);				/* debugging purposes */
  catch_fpe(FPE_DEFAULT);		/* floating point exception handling */

  usage_message=USAGE;
  argc=getargs(argc, argv, argtab);	/* parse arguments */
  if (check_arguments(argc) )
    die("%s\n", usage_message);
  patchmode= (opt_recover || (opt_patches[0]!='\0'));

  box=get_atoms(opt_infile, opt_protein, opt_bin, PDB_DEFAULT_FLAGS );

  if (opt_randomize || opt_ranval >= 0) {
      warn("\
NOTE: the randomizing options require atomic suface accessibilities in the 5th\n\
PDB field, but this is not checked. See the -a file.area option");
      ran_atoms(box, opt_ranval);      /* redistribute N,O over surface */
  }

  if(opt_infile)
    fclose(opt_infile);

  assert(box->atoms);			/* is set even if Str!=NULL */

  if(opt_donothing)			/* for timing purposes */
     return 0;

  get_template_spheres(templates_file);	/* get (or generate) the templates */


  pre_process_atoms(box, 
		    opt_extension, 
		    0.0 /*=threshold*/, 
		    0.0 /*=opt_polar_extension*/,
		    opt_pts_per_atom, 
		    opt_pt_density); 
  
  n=accessability(box);
  warn("%d atoms got completely buried (first time)\n", n);
  transfer_acc(box);		/* put areas from sh to atom->xyz[4] */
  warn("\
total, hydrophobic, hydrophylic area, hfob ratio: %.3f %.3f %.3f %.3f\n", 
       box->totarea, box->Carea + box->Sarea, box->Narea + box->Oarea, 
       (box->Carea + box->Sarea) / box->totarea );

  if(opt_areasfile) {
    print_areas(opt_areasfile, box);
    fclose(opt_areasfile);
  }
  if (opt_dotsfile) {
    print_dots(opt_dotsfile, box);
    fclose(opt_dotsfile);
  }

  n=atom_topology(box, opt_recover + patchmode );
  if(n)warn("deleted %d widow points\n", n);
  n=connections(box, patchmode );
  warn("found %d connections among boundries (first time)\n\n", n);


  if (opt_polar_extension != 0.0) {	/* have to do thing 2nd time */
    /* extend radii of some atoms, based on their accessability found so far */
    box->old_shells=box->shells;	/* for comparing later */
    box->shells=NULL;
    reset_atompointers(box);
    pre_process_atoms(box,		/* opt_radius, OBSOLETE */
		      opt_extension,	/* probe radius */
		      opt_threshold*M_4PI/100.0, /* percentage of full area */
		      opt_polar_extension,
		      opt_pts_per_atom, 
		      opt_pt_density); 
    n=accessability(box);			/* 2nd time */
    warn("%d atoms got completely buried (2nd time)\n", n);
    atom_topology(box, patchmode);	/* 2nd time don't need internal pts */
    n=connections(box, patchmode );
    warn("found %d connections among boundries (2nd time)\n", n);
  }

  /* patches stuff: */
  if(opt_recover)
    recover(stdout, box, opt_patdots);
  else if(opt_patches[0])
    do_patches(stdout, box, opt_patches);
  
  /* miscellaneous output */
  post_process_atoms(box, opt_reduction, opt_polar_reduction); 

  if (opt_outfile) { 
    print_triangulines(opt_outfile, box);
    fclose(opt_outfile);
  }
  if (opt_facesfile) {
    warn("written %d faces\n",
	 print_faces(opt_facesfile, box));
    fclose(opt_facesfile);
  }

  if (dmw_on) {				/* only then does it make sense */
    free_triangulation(box);
    free_rest(box);			/* misc. stuff */
  }

  return 0;
} /* main */
