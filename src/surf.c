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

/* routines for doing traversals of the whole triangulated surface */

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
#include "quilt.h"
#include "hash.h"

#define MAXNPATCHES 1024		/* just a guess (fails for 1rus */
					/* -th 10, but that's quite extreme) */
#define P_COLLECTED BIT(0)		/* patch found by collect_patches */

#define P_RECOVERED BIT(1)		/* patch not available for recovery */
#define P_PRELIM    BIT(3)		/* BIT(2) is HASWIDOWS */
#define P_DONT_RECOVER (P_RECOVERED | P_PRELIM)

static struct {
  box_p box;
  patch_p * patches_parray;
  int patches_index, maxnpatches;
  const char * patches_nature;
  int shared;
} local;

static int uniqify_atoms(surf_p surf) { 
  /* fills surf->atoms with uniq atom pointers that constitute patch */
  int i, natoms;
  atom_p *atoms;
  htable_p atom_table;
  hnode_p node;

  atom_table=Hcreate(NULL, surf->npatches, -1 );

  for (i=0; i<surf->npatches; i++)
    Hset(atom_table, surf->patches[i]->atom); /* just create table entry;
						 don't do anything with it */

  surf->natoms=natoms=atom_table->nkeys;
  atoms=surf->atoms=(atom_p*)MALLOC(natoms*sizeof(atom_p));

  for (i=0,node=Hsorted_keys(atom_table, NULL); node; node=node->next,i++)
    surf->atoms[i]=node->key;
  
  Hfree(atom_table, NULL);
    
  return natoms;
} /* uniqify_atoms */

static void sum_areas(surf_p surf) {
  /*
   * does    : puts totals of npoints and areas into newly allocated surf
   *
   * gets    : a surf
   *
   * affects : surf->npoints, surf->area
   *
   * returns : nothing
   *
   * warns   : 
   *
   * comment : can we trust the radii here ??
   */

  patch_p patch, *patches;
  int i, totnpoints=0, npatches;
  double totarea, area, radius;

  totarea=0.0;
  patches= surf->patches;
  npatches = surf->npatches;

  for (i=0; i<npatches; i++) {		/* collect areas by going to the */
    patch = patches[i];			/* patch on one atom */
    radius = patch->atom->RADIUS;	/* that atom's radius */
    area = SQR(radius)*patch->area;
    totarea += area; 
    totnpoints += patch->npoints;
  }

  surf->area = totarea;
  surf->npoints = totnpoints;

  return;
} /* sum_areas */

static surf_p make_surf(int npatches, patch_p * patches_parray) {
  /*
   * does    : allocates a surf for the newly found surface, and saves
   *   results of collect_patches. Sets areas
   *
   * gets    : results of collect_patches
   *
   * affects : nothing
   *
   * returns : a surf
   *
   * warns   : 
   *
   * comment : 
   */

  /* allocates a surf, and saves results of collect_patches in it. Areas */
  /* are set in sum_areas. surf is returned */
  surf_p surf;

  assert(npatches < local.maxnpatches);

  surf = (surf_p)CALLOC(1, sizeof(surf_t));
  surf->npatches = npatches;
  surf->patches = (patch_p*)MALLOC(npatches*sizeof(patch_p));
  memcpy(surf->patches, patches_parray, npatches*sizeof(patch_p));

  sum_areas(surf);
  uniqify_atoms(surf);
  return surf;
} /* make_surf */

static surf_p isolated_atom_patch(atom_p atp) {
  /* forms a surface out of an isolated atom, and returns it */

  surf_p surf = (surf_p)CALLOC(1, sizeof(surf_t));

  surf->npatches=1;
  surf->area = M_4PI*SQR(atp->RADIUS);
  return surf;
} /* isolated_atom_patch */

#define RIGHT_NATURE(ATOM) (local.patches_nature[0] == '\0' \
			    || strchr( local.patches_nature,\
				      (ATOM)->name[1]) != NULL)

static void collect_patches(patch_p patch) {
  /* try and extend this patch recursively, until we can't go any further */
  /* check if the patch we're moving to is of the right nature. To keep it */
  /* fast & slim, just keep track of the patches in one big global array + */
  /* index  */
  int i, nngbs;
  patch_p p;
  
  assert_not(patch->flags & P_COLLECTED);
  patch->flags |= P_COLLECTED;
  local.patches_parray[ local.patches_index ++ ] = patch;
  /* just store pointer to patch in huge array, and update later */

  nngbs=patch->nngbs;
  for(i=0; i<nngbs; i++ ) { 
    p=patch->ngbs[i];
    if ( !(p->flags & P_COLLECTED) 		/* not yet done */
	&& RIGHT_NATURE(p->atom)) 
      collect_patches(p);
  } /* for i */
  return;
} /* collect_patches */

static int print_surface(FILE *file, int nr, surf_p surf) { 
  /* prints one complete patch */
  int i, k, n;
  list_p l;
  shell_p sh, shells;
  atom_p atom, atoms;
  residue_p res;
  char ch, ss;
  static char secstr[3]={ '-', 'a','b', };
  double area;
  char *buf;

  atoms = local.box->atoms;
  shells= local.box->old_shells;	/* for areas */
  fprintf(file, "\
# %d %d atoms %d patches %d points area %.3f %d bndries ABC %.3f %.3f %.3f\n", 
nr, surf->natoms, surf->npatches, surf->npoints, surf->area,
	  surf->nbndries, 
	  surf->ellips_axes[0], surf->ellips_axes[1], surf->ellips_axes[2]); 
  i=0;
  for (l=surf->bndries; l; l=l->next, i++) {
    fprintf(file, "# bndry %d, %ld patches\n",
	    i, ((ulong*)l->ptr)[0]);
  }  
  k=0;
  buf=(char*)ALLOCA(24*surf->natoms);
  for (i=0; i<surf->natoms; i++)
    k += sprintf(buf+k, "%d ", surf->atoms[i] - atoms);
  fill(buf, "% ", " ", 72);
  fprintf(file, "%s", buf);

  k=0;
  for (i=0; i<surf->natoms; i++) {
/*    patch=surf->patches[i]; # may not be uniq atom */ 
    atom=surf->atoms[i];
    n=atom - atoms;			/* index of this atom */
    sh= &shells[n];			/* the sh from before expansion */
/*    area= SQR(atom->RADIUS)*patch->area; */
    area= SQR(atom->RADIUS)*sh->area;	/* the sum */
    res=atom->residue;
    ss =  secstr [ res->flags.SS ] ;
    ch=res->chain->name;
    k+= sprintf(buf+k, "%c %c%hd%c@%s=%.0f%c;", 
	    res->type, ch, res->nr, res->inscode, atom->name,
	    100.0*area, ss);
  } /* for i */
  fill(buf, " ", ";", 72);	/* only break every atoms, ie. on ';' */
  fprintf(file, "%s", buf);
  AFREEA(buf);
  return 0;
} /* print_surface */

int print_surfaces(FILE * file, surf_p Surfaces, const char *nature) {
  /*
   * does    : prints the whole thing to stdout
   *
   */
  int k, totnpatches, totnpoints, nsurfaces;
  surf_p surf;
  double totarea;
/*  char line[256]; */

  nsurfaces=totnpatches=totnpoints=0;
  totarea=0.0;

  /* first get totals */
  for(k=0,surf=Surfaces; surf; surf=surf->next, k++) {
    totnpatches += surf->npatches;
    totarea += surf->area;
    totnpoints += surf->npoints;
    nsurfaces++;
  }

  fprintf(file, "\n\
surface '%-3s': %3d surfaces %4d atom patches %5d points area %10.3f\n"
	 , nature[0] ? nature : "any", 
	 nsurfaces, totnpatches, totnpoints, totarea);

  for(k=0,surf=Surfaces; surf; surf=surf->next, k++)
    print_surface(file, k, surf);
  return nsurfaces;
} /* print_surfaces */

#define SURF_PRECISION 1000.0		/* ignore last three decimals */
					/* when comparing */

#define AREA_OF(X) ((int) ((*(surf_p*)(X))->area * SURF_PRECISION ) )
#define ATOM_OF(X)  ((*(surf_p*)(X))->patches[0]->atom)
static int surf_cmp(const void *a , const void *b) {
  /* sort surfs by size, but ignore differences smaller than */
  /* 1.0/SURF_PRECISION; after this, sort on atom nr of first atom  */ 
  int i;
  i=AREA_OF(b) - AREA_OF(a);
  if(i)
    return i;
  return (int)((ulong)ATOM_OF(a) - (ulong)ATOM_OF(b)); /* contiguous ! */
} /* surf_cmp */

#undef  ATOM_OF
#define ATOM_OF(X) ((*(patch_p*)(X))->atom)
static int patch_ptr_cmp(const void *a, const void *b) {
  if (ATOM_OF(a) < ATOM_OF(b))
    return -1;
  if (ATOM_OF(a) > ATOM_OF(b))
    return 1;
  return 0;
} /* patch_ptr_cmp */

static surf_p sort_surfaces(surf_p surfaces) {
  /*
   * does    : sorts patches in order of descending size, and the atomic
   *   patches contributing to it to ascending atom nr. Also gives
   *   the surfaces their rank nr.
   *
   * gets    : a list of surfaces (patches)
   *
   * affects : order of surfaces->patches, and order of surfaces; a
   *
   * returns : resulting list
   *
   * warns   : not
   *
   * comment : 
   */
  int i;
  surf_p surf, surfs;

  for(surf=surfaces; surf; surf=surf->next) /* order atoms within patch */
    qsort(surf->patches, surf->npatches, sizeof(patch_p), patch_ptr_cmp);

  /* sort by size and atomnr: */
  surfs=(surf_p)Lsort((list_p)surfaces, surf_cmp); 

  for(i=0,surf=surfs; surf; surf=surf->next)	/* add rank numbers */
    surf->rank=i++;

  return surfs;
} /* sort_surfaces */

static void reset_patch_flags(surf_p surfaces) {
  /* resets the 'done' flags of SURFACES  */
  int j;
  surf_p surf;

  for (surf=surfaces; surf; surf=surf->next)
    for (j=0; j<surf->npatches; j++) { 
      if (surf->patches)
        surf->patches[j]->flags &= ~P_COLLECTED;	/* leave other flags
                                                           intact */
      else { 
        warn("*** suspect: surface without patches, must be isolated atom, right? Removing surface altogether.");
        Lexcise((list_p*)&surfaces, (list_p)surf);
      }
    }
  return;
} /* reset_patch_flags */

static int surf_ellips(surf_p surf) {
  /*
   * does    : calcs lengths of ellipsoid axes containing all surf atoms.
   *
   * gets    : surf list 
   *
   * affects : surf->ellips_axes[0-2] of first surf. 
   *
   * returns : 1 if ellipsoid too flat or so
   *
   * warns   : if ellipsoid is degenerate
   *
   * comment : uses just the atom centers, not the dots (may be in future)
   */
  int i,n, err;
  float **xyz;
  double axes[3];

  n=surf->natoms;
  xyz=(float**)ALLOCA(n*sizeof(float*));
  for (i=0; i<n; i++)
    xyz[i]=surf->atoms[i]->xyz;

  err=ellipsaxes_3D(n, xyz, 1.0, axes);

  if(!err)
    for (i=2; i>=0; i--)
      surf->ellips_axes[i] = axes[i]; /*  + 3.3; */

  AFREEA(xyz);
  return err || axes[0] < 0.01 ;
} /* surf_ellips */


surf_p get_surfaces(box_p box, const char * nature) {
  /*
   * does    : traverse the surface of specified nature, and assembles one
   *   or more contingous pieces out of it. The atom faces are those of
   *   after any polar expansion
   *
   * gets    : atoms, etc. and nature (eg "NO" for oxygen and nitrogen).
   *   If nature == "", nature is not checked, giving back all surface
   *   components
   *
   * affects : nothing
   *
   * returns : a linked list of connected surfaces of specified nature. A
   *   surf_t contains an array of pointers to the atomic patches
   *   constituting it. The components of the surface are sorted in order 
   *   desceding area
   *
   * warns   : 
   *
   * returns: list of surfaces
   * comment : ignores WIDOW atoms (@change this)
   */
  int ii;
  atom_p atom;
  shell_p shell;
  patch_p patch, *localarray;
  surf_p Surfaces=NULL, surf;

  local.box=box;
  local.patches_nature=nature;
  local.maxnpatches=box->natoms + 50; /* wild guess */
  /* assume there will usually be (much)less than natoms patches, unless */
  /* it's a small structure, so add the 50 */

  localarray= (patch_p*)ALLOCA(local.maxnpatches * sizeof(patch_p));
  bzero((char*) localarray, local.maxnpatches*sizeof(patch_p));
  local.patches_parray = localarray;  /* global pointer into it */

  for (ii=0; ii<box->natoms; ii++) {
    shell=box->shells+ii;
    /*    radius=atp->RADIUS; */
    if (shell->flags & BURIED)
      continue;
    atom=box->atoms + ii;

    if( nature[0]			/* if at all a nature is specified */
       && ! strchr(nature, atom->name[1])) /* wrong type of surface */
      continue;

    if (shell->flags & ISOLATED ) {	/* degenerate case */
      surf = isolated_atom_patch(atom);
      Lprepend(surf, Surfaces);
      continue;
    }

    /* go look for a patch we didn't reach yet: */
    for(patch = shell->patches; patch; patch=patch->next) {
      if ( (patch->flags & P_COLLECTED))
	continue;

      /* found new, unseen patch; use as entry into surface of this type */
      local.patches_index = 0;	/* reset ! */
      collect_patches(patch);		/* traverse surface recursively */
      assert( local.patches_index <= local.maxnpatches); 
					/* otherwise array corrupt */
      surf=make_surf(local.patches_index, 
		     local.patches_parray); /* save */
/*     surftop(surf); not ready */
      Lprepend(surf, Surfaces);		/* put into list */
    } /* for patch */
  } /* for ii */
  reset_patch_flags(Surfaces);
  AFREEA(localarray);
  return sort_surfaces(Surfaces); 
} /* get_surfaces */

static patch_p match_patch(shell_p shell, shell_p orgshell, patch_p patch) {
  /*
   * does : find out the original patch in the original triangulation that
   *   corresponds to patch
   *
   * gets    : the SHELL and PATCH after the expansion, the ORGSHELL before
   *  the expansion, 
   *  
   * affects : 
   *
   * returns : the patch in the original triangulation 
   *
   * warns   : 
   *
   * comment : works by looking at point numbers, and therefore 
   *   the nr's of points per sphere have to be same in both triangulations
   */
  int i, j;
  bndry_p bndry;
  patch_p orgpatch;
  SHORT point;


  if ( (shell->flags | orgshell->flags) & MPATCHES)
    ;					/* deal with that in a moment */
  else {				/* just one face */
    assert_not( orgshell->patches->next);
    if (orgshell->patches->flags & P_DONT_RECOVER) { /* bad luck: used before */
      local.shared++;
      return NULL;
    }
    return orgshell->patches;		/* the only one */
  }

  /* else : one or both of the shells contain more than one face */
  for(bndry=patch->bndries; bndry; bndry=bndry->next) {
    for (i=0; i<bndry->length; i++) {	/* only need to check bndry points ! */
      point=bndry->path[i];
      for(orgpatch=orgshell->patches; orgpatch; orgpatch=orgpatch->next) {
	if (orgpatch->flags & P_DONT_RECOVER)
	  continue;
	for (j=orgpatch->npoints-1; j>=0; j--) /* faster: starts in middle ! */
	  if ( point==orgpatch->points[j] ) /* found it */
	    return orgpatch;
      }	/* for orgpatch */
    } /* for i< length */
  } /* for bndry */
  /* not found: */
  local.shared++;
  return NULL;
} /* match_patch */

static void mark_recovered(patch_p patch) {
  /* mark patch as not available for further recovery from other patches */
  assert_not(patch->flags &  P_DONT_RECOVER );
  patch->flags |= P_RECOVERED;
}
static void mark_prelim(patch_p patch) {
  /* mark patch as not available for further recovery from other patches */
  assert_not(patch->flags & P_DONT_RECOVER );
  patch->flags |= P_PRELIM;
}

static int translate_patches(box_p box, int norgpatches, patch_p * patches, 
			     patch_p * newpatches) {  
  /*
   * does    : given an array of atom patches, finds the corresponding
   *   atom patches from before the hydrophylic extension
   *
   * gets    : patches
   *
   * affects : contents of newpatches, and sets patch->flags & P_RECOVERED
   *
   * returns : nr of new patches
   *
   * warns   : dies if finds anything with multiple patches on one atom
   *
   */
  int i, atomnr, n;
  shell_p shell, shells, orgshell, org_shells;
  patch_p patch, orgpatch;
  
  shells=box->shells;
  org_shells=box->old_shells;

  n=0;
  for (i=0; i<norgpatches; i++) {
    patch=patches[i];			/* atomic patch belonging to surf */
    atomnr= patch->atom - box->atoms;
    shell=shells+atomnr;
    orgshell=org_shells+atomnr;
    if (orgshell->flags & ISWIDOW)
      continue;
    orgpatch=match_patch(shell, orgshell, patch); /* find corresponding orgpatch */
    if (!orgpatch) { 
/*       warn("corresponding atom face not found: %s\n",
 * 	   atom_spec(local.box, atomnr)); /# happens quite often, dont warn #/
 */
      /* assert(patch->npoints <= 2); NO: can easily happen at border */
      continue;				/* what else ? */
    }
    mark_prelim(orgpatch);
    newpatches[n++]=orgpatch;
  } /* for i < noldpatches */
  return n;				/* nr of new patches */
} /* translate_patches */

static int patch_there(patch_p patch, int npatches, patch_p*newpatches) {
  /*
   * does    : check patch is already among newpatches
   *
   * gets    : the patch to be tested, index NPATCHES into array NEWPATCHES
   *
   * affects : local.shared is incremented if trying to recover a patch
   *   that is already 'in use' (recovered) by another prelim-patch
   *
   * returns : 1 if it's already there, 0 if not
   * comment : 
   */
  int k;

/*  assert_not(patch->flags & P_COLLECTED); */
  if (patch->flags & P_DONT_RECOVER) {	/* already used once before ! */
    local.shared++;
    return 1;
  }
  for (k=npatches-1; k>=0; k--)		/* check those accepted */
    if (patch==newpatches[k])		/* already present */
      return 1;

  return 0;
} /* patch_there */

static surf_p extend_surface(surf_p surf, box_p box) {
  /*
   * does    : tries to extend a surface so as to recover area that got
   *   buried by extending the hydrophylic radii
   *
   * gets    : a surface to work on; the original surfaces is obtained
   *   from box
   *
   * affects : 
   *
   * returns : a newly allocated surface
   *
   * warns   : 
   *
   * comment : works by looking at the original neighbours of a particular
   *   atom patch
   */

  int i,j, npatches, atomnr, nnew;
  patch_p orgpatch, ngb_patch, orgpatches[ MAXNPATCHES ];
  shell_p shell, org_shells;
  atom_p atoms;

  npatches=surf->npatches;
  if ( npatches >= MAXNPATCHES )	/* unlikely, but warn */
    DIE("\
number of atoms in patch too large (%d); change MAXNPATCHES and recompile,\n\
or run with lower threshold (-th),", npatches);

  npatches=translate_patches(box, npatches, surf->patches, &orgpatches[0]);
  if (npatches == 0)			/* found no match ... */
    return NULL;
  /* all surf->patches[0] are now translated to what used to be.*/
  /* orgpatches are the original patches, i.e. those from before expansion */
  atoms=box->atoms;
  org_shells=box->old_shells;
  
  nnew=npatches;
  /* look for new patches, being only the first ngbs of the original patches */
  for (i=0; i<npatches; i++) {
    orgpatch=orgpatches[i];		/* atomic patch belonging to surf */
    for (j=0; j<orgpatch->nngbs; j++)  { /* check all atom's ngbs */
      ngb_patch=orgpatch->ngbs[j];	/* new patch to try */
      atomnr = ngb_patch->atom - atoms;	/* index of atom, and of shell */
      shell= org_shells + atomnr;
      assert(shell->nngbs && !(shell->flags&BURIED));
      if( (shell->flags & POLAR )	/* ignore  */
	 ||
	 patch_there(ngb_patch, nnew, orgpatches)) /* already seen */
	continue;
      orgpatches[nnew++]=ngb_patch;	/* add it */
      mark_recovered(ngb_patch);
    }	/* for j */
  } /* for i */
  return make_surf(nnew, orgpatches);
} /* extend_surface */

#ifdef undefined
static void purge_multiples(surf_p surfaces) {
  /*
   * does    : gets rid of atomic patches that are shared between surfaces
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
  int i, len, nmultiples;
  surf_p surf, *victims;
  patch_p multiples[ MAXNPATCHES ];

  len=Llength(surfaces);
  victims=(surf_p*)ALLOCA(len*sizeof(surf_p));

  nmultiples=0;
  for(surf=surfaces; surf; surf=surf->next) {
    
  }
  AFREEA(victims);
  return;
} /* purge_multiples */
#endif /* undefined */

static int print_patch_dots(surf_p surfaces, int patdots) {
  /* print dots of the first patdots patches to patN.dots */
  int i, j, n;
  surf_p surf;
  atom_p atom;
  FILE *file;
  char filename[100];

  n=0;
  for (i=0,surf=surfaces; i<patdots; i++,surf=surf->next) {
      sprintf(filename, "/junk/lijnzaad/pat%d.dots", i);
      warn("writing to %s\n", filename);
      file=fopen(filename, "w");
      if(!file)
	perror(filename),exit(1);
      fprintf(file, ".# patch %d\n", i);
      for (j=0; j<surf->natoms; j++) {
	atom=surf->atoms[j];
	fprintf(file, ".cmov %.3f %.3f %.3f\nx%d\n", 
		atom->xyz[0], atom->xyz[1], atom->xyz[2], j);
	n+=print_atom_dots(file, atom, ".dot", "");
      }
      fclose(file);
  } /* for i < patdots */
  return n;
} /* print_patch_dots */

int recover(FILE * file, box_p box, int patdots) {
  /*
   * does    : find patches for both triangulations, and tries to recover
   *   the bits of patches that got buried in the 2nd thing
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
  int i;
  surf_p surf, newsurf, Surfaces, newsurfaces;

  Surfaces=get_surfaces(box, "CS");	/* a whole, sorted, list of surfaces */
  /* this is done on the topology with polar atoms extended */
  
  newsurfaces=NULL;

  for(i=0,surf=Surfaces; surf; surf=surf->next, i++) {
    newsurf=extend_surface(surf, box);	/* recover what got buried */
    if(!newsurf) {			/* no matching patch */
/*      warn("could not find matching atom face ...\n"); is already warned */
      continue;
    }
    if(newsurf->natoms >= 5)	/* 4 still tricky; ellipsoid easily singular */
      if ( surf_ellips(newsurf) )
	warn("can't use ellipsoid (%d %d %.3lf)\n", 
	     newsurf->npatches, newsurf->npoints, newsurf->area);
    newsurf->old=surf;			/* establish correspondence */
    Lprepend(newsurf, newsurfaces);
  }

  /* (after this, newsurfaces is sorted more or less in ascending order) */
  if(local.shared)
    warn("\
%d atomic faces border several patches; assigned to largest\n\
\n",
	 local.shared);
  
  newsurfaces=sort_surfaces(newsurfaces);

  /* print new ranks: */
  fprintf(file, "# new patch ranks, after recovering\n");
  for(i=0,surf=newsurfaces; surf; surf=surf->next, i++) { 
    assert(i==surf->rank);
    fprintf(file, "#%d: %d (%.3f) WAS #%d: %d (%.3f)\n",
	   surf->rank, surf->npatches, surf->area, 
	   surf->old->rank, surf->old->npatches, surf->old->area);
  }

  print_surfaces(file, newsurfaces, "after recover");
  if (patdots)				/* dump dots to files */
    print_patch_dots(Surfaces, patdots);

  if(dmw_on) { 
    free_surface(Surfaces);
    free_surface(newsurfaces);
  }
  return 0;
} /* recover */

/* following doesn't really belong here but in io.c / debug.c ? */
static void write_point(FILE *file, int n,
			atom_p atom, int point) { 
  int i;
  shell_p shell;
  double xyz[3], r;
  static char string [80];

  shell=atom->pointer;
  r=atom->RADIUS;
  if (n==0)
    string[0]='\0';
  sprintf(string, "%s %d-%d", string, (int)(atom - local.box->atoms),
	  point);

  for (i=0; i<3; i++) 
    xyz[i]=atom->xyz[i] + r*shell->sphere->coordinates[3*point +i];
  switch (n) {
  case 0:
    fprintf(file, ".m %8.3f %8.3f %8.3f\n", xyz[0], xyz[1], xyz[2]);
    break;
  case 1:
    fprintf(file, ".d %8.3f %8.3f %8.3f\n", xyz[0], xyz[1], xyz[2]);
    break;
  case 2:
    fprintf(file, ".d %8.3f %8.3f %8.3f # %s\n", xyz[0], xyz[1], xyz[2],
	    string);
    break;
  default:
    assert(0);				/* shouldn't reach this point */
  } /* switch */
  return;
} /* write_point */

static int atomfaces(FILE *file) { 
  /* writes faces on atoms to file */
  int i, j, k, ntot;
  SHORT *pointflags;
  box_p box;
  shell_p shell;
  atom_p atom;
  sphere_p sphere;
  face_p face;

  fprintf(file, "### faces on atoms\n");
  box=local.box;
  ntot = 0;
  for (i=0; i<box->natoms; i++)  {
    atom=box->atoms+i;
    shell=atom->pointer;
    if (!shell)
      continue;
    if(INVALID_SHELL(shell))
      continue;
    fprintf(file, "# atom %d\n", i);
    sphere=shell->sphere;
    pointflags=shell->pointflags;
    for (j=0; j<sphere->nfaces; j++)  { 
      face = sphere->Faces+j;
      for (k=0; k<3; k++)
	if ( pointflags [ face->points[k] ] & (BURIED | UNWANTED) ) 
	  goto NEXT_FACE;
      ntot++;
      for (k=0; k<3; k++)
	write_point(file, k, atom, face->points[k]);
    NEXT_FACE:
      ;
    }
  } /* for i */
  return ntot;
} /* atomfaces */

static void set_seams(int nseams, seam_p seams) {
  int i, n, end, oend;
  con_p con;

  for (i=0; i<nseams; i++) {
    n=(i+1)%nseams;			/* index of next connection */
    end= seams[i].end;
    con= seams[i].con;
    seams[i].from= con->flags[ end ].offset;
    oend=end;
    end = !seams[n].end;		/* begin of the next seam */
    assert( seams[n].con->bndries[ end]== con->bndries[oend] ); /* of course */
    seams[i].till = seams[n].con->flags[ end ].offset;
  } /* for i < nseams */
  return;
} /* set_seams */

static int print_face(FILE *file, con_p con, int end) { 
  /* prints face 'to the left' of con's end */
  int i,j,n, nseams;
#define MAXNSEAMS 4			/* or 3 ? */
  seam_t *seam, seams[ MAXNSEAMS ];
  bndry_p bn;

  n=find_seam(con, end, &nseams, MAXNSEAMS, seams);
  if (n !=0) {
    assert(n==4);
    assert(nseams == 2);		/* but could be more ... find out */
    return 1;
  }
  set_seams(nseams, seams);
  n=0;
  for (i=0; i<nseams; i++) {
    seam=seams+i;
    bn=seam->bndry;
    for (j=0; j<seam->length; j++) {
      write_point(file, n,
		  bn->patch->atom, bn->path[(seam->from+j)%bn->length]);
      n++;
    }
    seam->con->flags [ seam->end ].done =0; /* mark (unmarked is 1 ! */
    seam->con->flags [ seam->end ].seen =0; /* reset */
  } /* for i  */
  assert(n==3);
  return 0;
} /* print_face */

static int interfaces(FILE *file) { 
  /* writes faces between atoms to file */
  int i, ii,n,end;
  atom_p atoms;
  shell_p shell;
  patch_p patch;
  bndry_p bndry;
  list_p l;
  con_p con;
  box_p box;


#define POINT_SPEC(I) con->bndries[(I)]->path[con->flags[(I)].offset],\
  atom_spec(box, con->bndries[(I)]->patch->atom - atoms)

  fprintf(file, "### faces between atoms\n");
  box=local.box;
  atoms=box->atoms;
  for (ii=0; ii<box->natoms; ii++) {
    shell=box->shells+ii;
    if ( INVALID_SHELL(shell))
      continue;
    for(patch=shell->patches; patch; patch=patch->next)
      for(bndry=patch->bndries; bndry; bndry=bndry->next)
	for (i=0; i<bndry->length; i++) 
	  for (l=bndry->cons[i]; l; l=l->next) { 
	    con=l->ptr;
	    end = (con->bndries[0]==bndry);
	    if (con->flags[end].done==0) /* marking is done with 0 here ! */
	      continue;
	    if( print_face(file, con, end) )
	      warn("no proper face around points %d %s and %d %s\n",
		   POINT_SPEC(0), POINT_SPEC(1));
	    else
	      n++;
	  } /* for l */
  } /* for ii */
  return n;
} /* interfaces */

int print_faces(FILE *file, box_p box) { 
  local.box=box;
  return atomfaces(file) + interfaces(file);
} /* print_faces */


void free_surface(surf_p surface) {
  surf_p s;
  void * next;

  for(s=surface; s; s=next) {
    next=s->next;
    FREE(s->patches);			/* array of pointers to patches */
    FREE(s->atoms);			/* array of pointers to atoms */
    FREE(s);
    Lfree(s->bndries, _free);
  }
  return;
} /* free_surface */
