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

#ifndef _QUILT_H_
#define _QUILT_H_ 1

#include "compat.h"

/* change to needs: */

/* where template spheres are to be read from, or written to: */
/* TEMPLATES_DIR and TEMPLATES_FILE can be defined from the compilation */
/* commandline; if they're not defined, use following */
#define DFLT_TEMPLATES_DIR ""		/* (current directory, else home) */
/* if not defined from command line, use following: */
#define DFLT_TEMPLATES_FILE ".quilt-spheres.bin" /* filename */

#include "utils.h"			/* for uchar etc. typedefs */
#include "list.h"			/* for list typedefs */
#include "pdb.h"			/* for pdb typedefs */

#ifdef __GNUC__				/* gnu cc, knows directive 'inline' */
#  ifdef DEBUG
#  define inline
#  endif
#else  					/* others don't */
#  define inline
#endif

/* constants: */
/* floats (tolerances etc.) */
#define DEF_PROBE_RADIUS 1.40		/* default radius of probe */

#define COS_SMALL_ANGLE 0.99		/* cosine of 8.1 degrees. Used in */
					/* zip_pair to rule out too sharp */
					/* corners */ 
/* tolerances: */
#define TINY_ANGLE_TOL 0.005		/* fractional tolerance for */
					/* judgeing an angle too tiny */
#define INTERSECT_TOL 0.001		/* atom overlap should be at least */
					/* this for intersect_pair() to be */
					/* called */

#define SWALLOW_TOL 0.001		/* if the difference with being */
					/* completely swallowed is less than */
					/* this, consider it to be so */

#ifdef SPH3
#  define MAX_NGBS_PER_PNT 7		/* maximal coordination number */
#else
#  define MAX_NGBS_PER_PNT 6		/* maximal coordination number */
#endif /* SPH3 */

#define NDAT_ITEMS 3			/* 3 coordinates per point */

/* types of spheres: */
#define SPH_ICOS    0			/* icosahedron */
#define SPH_PDOD    1			/* pentakisdodecahedron */
#define SPH_ZIEL    2			/* for sph3. (not used actually) */
#define SPH_UNKNOWN 3			/* when read from file */

#define N_TMPLT_SPHERES 7		/* number of template spheres */
#define MAX_TMPLT_PNTS 256 /* 252 */	/* highest number of points on */
					/* template spheres (== 10*d*d+2; */
					/* 30*d*d+2 if dodecaedron-based) */
#define MAX_PATHLENGTH MAX_TMPLT_PNTS/2
#define MAX_NCELLS 343			/* max nr of cells per sphere */
#define MINRADIUS 3*2			/* index in limits array */
#define MAXRADIUS 3*2 + 1		/* index in limits array */

#define MAX_NNGBS 256			/* max nr of neibors to any atom */

#define MAX_NCONS 16			/* max nr. of cons out of a point */

/* flags  ( #define BIT(n) (1<<(n)) is already defined) */
/* general purpose flag */

/* flags for points */
#define BURIED 		BIT(0)		/* points/edges buried */
					/* DON'T CHANGE THIS ONE !!! */
					/* set in acc.c:calc_acc() */
#define BNDRY		BIT(1)		/* part of a bndry */
					/* set in top.c:find_bndry(), */
					/* top.c:hide_flags()  */
#define PATCH		BIT(2)		/* marked as part of patch */
					/* set in top.c:find_bndry(), */
					/* top.c:mark_interior()  */
#define ISOLATED	BIT(3)		/* single pt bndry & patch */
					/* for shell: completely exposed ! */
					/* set in top.c:find_bndry() */
#define TWOWAY		BIT(4)		/* bndry passes twice (or more) */
					/* through through this pt */ 
#define UTURN		BIT(5)
#define SHARED		BIT(6)		/* point shared by 2 bndries */
#define UNWANTED 	BIT(7)		/* to be neglected */

/* flags for bndries, patches, and/or shells */
#define COORDINATES	BIT(8)

/* #define BURIED	BIT(0)  should be same everywhere */
#define DONE		BIT(0)		/* patch has been traversed once */

#define ISWIDOW		BIT(1)		/* no neighbours (shell only) */
					/* set in acc.c:calc_acc() */
#define HASWIDOWS	BIT(2)		/* atom contains a widow point */
/* #define ISOLATED	BIT(3)		 completely exposed atom */
					/* set in acc.c:purge_ngbs(), */
					/* acc.c:isolated_atom()  */
#define MBNDRIES	BIT(4)		/* shell contains > 1 bndrys */
					/* top.c:find_patches */
#define MPATCHES	BIT(5)		/* shell contains > 1 patches */
#define COMPLEX		(MBNDRIES | MPATCHES) /* both */
#define POLAR		BIT(6)		/* atom is polar */
#define REMOVE_WIDOWS	BIT(8)		/* to flag need for remove_widows() */

/*  */
#define RETOP 		BIT(9)		/* shell has to triangulated again */
#define STITCHED	BIT(10)		/* shell does not have to be "  " */

#define INVALID_SHELL(S)  ( ((S)->flags) & (BURIED | ISOLATED | ISWIDOW) )
#define VALID_SHELL(S)  (! INVALID_SHELL((S)))

/* macros: */
#define FASTMOD_NTIMES 3		/* nr of times to repeat row of 0..n */

#define FASTMOD(W, NGBS) fast_mod[ (NGBS) ][ (W) ] 

/* #define MOD(A,B) FASTMOD((A),(B)) */

#define MOD(A,B) ((A)%(B))		/* fast_mod thing is bugged */

/* int fast_mod_func(int, int); */
/* #define MOD(A,B) fast_mod_func((A),(B)) */

/* typedefs: */
#include <limits.h>

typedef unsigned short SHORT;
#define POINT_MAX USHRT_MAX

/* to obtain +- 25 % reduction in memory usage, use following instead: */

/* typedef unsigned char SHORT;
 * #define POINT_MAX UCHAR_MAX
 */

typedef unsigned short EDGE_ID;
#define EDGE_MAX USHRT_MAX

#define AVAIL (POINT_MAX)		/* can't use that as point id */
#define UNAVAIL (POINT_MAX-1)		/* neither this */
#define MAXNROFPOINTS MIN((POINT_MAX-2), ((EDGE_MAX+6)/3))
					/* since 3N - 6 edges !! */
typedef struct point_ {
  float *xyz;
  float *ngbs_xyz[MAX_NGBS_PER_PNT];	/* coords of ngb points */
  SHORT nr;				/* its number; mainly for debugging */
  SHORT flags;
  SHORT nngbs;				/* number of neibors or edges */
  SHORT ngbs	 [MAX_NGBS_PER_PNT];	/* ids of ngbs */
  EDGE_ID edges	 [MAX_NGBS_PER_PNT];	/* ids of edges */
/*  was: SHORT edges [MAX_NGBS_PER_PNT]; */
  SHORT whoamis  [MAX_NGBS_PER_PNT];	/* index of this point in its */
  float area;
  struct point_ 
    *Ngbs	 [MAX_NGBS_PER_PNT];	/* pointers to neighbour points */
} point_t, *point_p;

typedef struct edge_ {			/* little used, but available */
  SHORT ends[2];			/* endpoints; always ends[0]<ends[1] */
  SHORT flags; 
  union sw_ {
    SHORT s; 
    uchar w[2]; 
  } sw; 
#define face_with sw.s			/* when doing find_faces */
#define whoami sw.w			/* othewise */

  float direction[3];
#define midpoint direction		/* needed during tr._sphere() */
  float length2;
  float arclen;
} edge_t, *edge_p;

typedef struct face_ {			/* not used, but available */
  ushort edges[3];
  ushort points[3];
/* was:  SHORT edges[3]; */
  SHORT flags;
  float area, midpoint[3];
} face_t, *face_p;

typedef struct sphere_ {
  int npoints, nedges, nfaces, order, maxringlen;
  int dim;				/* dimension of grid */
/*  int ncells; */
  int type;				/* SPH_ICOS or SPH_PDOD */

/*   point_p * cells;			/# array of length dim^3, holding #/
 * 					/# linked lists of poinst in cell #/
 *   uchar * pts_per_cell;			/# nr of atoms per cells #/
 */

  float * coordinates;			/* coordinates of all points, as  */
  point_p Points;			/* series of float[3]'s */
  edge_p Edges;
  face_p Faces;

  double edgelengths[3];		/* min, max, ideal edge length */
  double arclengths[3];			/* min, max arclen, ideal length */
  double beta, halfalpha;		/* cos(alpha)/cos(alpa/2); alpha/2 */
  double sinbeta;
  /* halfalpha = 0.5*maxarclen, beta is cos(minarclen)/cos(halfalpha) */
#define ideal_arclen arclengths[2]

} sphere_t, *sphere_p;

typedef struct bndry_ {			/* boundary of a patch */
  struct bndry_ *next;			/* next bndry belonging to patch */
  struct patch_ *patch;			/* patch it belongs to */
  SHORT *path;				/* point numbers */
  list_p * cons;			/*  */
  union { 
    struct shell_ ***ngbs;		/* during connect_2_bndries */
    float * normals;			/* during stitch */
  } pt;					/* shells (NULL terminated) */
  union { 
    struct shell_ ** shell;		/* during connect_bndries() */
    float *coordinates;			/* during stitch() */
    struct bndry_ ** bndry;		/* during connect_elements() */
  } ngbs;

  SHORT length;				/* length of path */
  SHORT ncons;				/* nr of connections */
  SHORT nngbs;				/* nr of bndry ngbs */ 
  SHORT flags;
} bndry_t, *bndry_p;

typedef struct patch_ {			/* only on one atom */
  struct patch_ * next;
  bndry_p bndries;
  atom_p atom;				/* to which it belongs */
  struct patch_ ** ngbs;		/* array of ngbs  */
  SHORT nngbs;				/* nr of ngbs */
  SHORT npoints;			/* all points of this patch */
  SHORT internals;			/* just the internals */
  SHORT flags;
  SHORT * points;			/* only internal points ! */
  float area;				/* area as if R where 1.0 ! */
} patch_t, *patch_p;

typedef struct shell_ {                 /* partial triangulation of atom;
                                           holds all the additional
                                           information per atom */ 
  union {
    struct shell_ * next;		/* next atom in grid cell */
    list_p *** conlists;
  } next;
  sphere_p sphere;			/* what sort of sphere */
  atom_p atom;				/* the atom it belongs to */
  struct shell_ ** ngbs;		/* neigbouring shells */
  patch_p patches;
  SHORT *pointflags;			/* point->flags per point */
  SHORT nngbs;				/* nr of neibors */
  SHORT flags;				/* eg. BURIED */
  SHORT nburied;			/* nr of atoms buried */
  float area;				/* area if it were a unit sphere */
} shell_t, *shell_p;

typedef struct surf_ {			/* usually covers more than one atom */
  struct surf_ * next;
  int npatches;				/* nr of patches forming surface */
  int npoints;				/* nr of points sampling it */
  patch_p *patches;			/* array of ptrs to patches */
  atom_p *atoms;			/* array of ptrs to atoms (unique) */
  double area;				/* calc'd in sum_area; real radii now*/
  int natoms;				/* nr. of unique atoms */
  int rank;				/* rank in sorting */
  struct surf_ * old;			/* pointer to equivalent surf before */
  float ellips_axes[3];
  int nbndries;
  list_p bndries;
} surf_t, *surf_p;			/* rescovering buried area */

typedef struct box_ {
  atom_p atoms;
  int natoms;
  structure_p structure;
  float * original_radii;
  shell_p * cells;
  shell_p shells;
  shell_p old_shells;			/* for 2nd time around ... */
  int ncells;
  int dims[3];				/* dimensions: nr of x,y and z boxes */
  float maxrad;				/* maximum radius */
  float spacing;
  float limits[8];			/* min-x, max-x; min-y, max-y; etc. */
  int nburied;
  double totarea,Narea,Oarea,Carea,Sarea;
  char name[20];
} box_t, *box_p;

#ifdef undefined
typedef struct con_ {
  bndry_p bndries[2];
  SHORT offset[2];
  float direction[4];			/* x, y, z (normalized), r */
} con_t, *con_p;
#endif /* old version */

typedef struct con_ {
  bndry_p bndries[2];
  struct {				/* rever to as con->flags[1].offset */
    SHORT seen:1, done:1, offset:14;
  } flags[2];
  float direction[4];			/* x, y, z (normalized), r */
} con_t, *con_p;

typedef struct seam_ {			/* old type; had to keep it small */
  con_p con;
  short end, length; 
  short from, till;
  bndry_p bndry;			/* always bndry==con->bndries[end] */
  SHORT *first, *last;
} seam_t, *seam_p;

typedef struct br_ {
  struct br_ * next, *prev;
  list_p *where;
  bndry_p bndry;
  shell_p shell;
  int offset;
  const float * coord, *normal;
  uint flags;
  list_p *conlist;
  list_p **conlists;
  ushort *pointflags;
  double atomradius;
  double dir[4]				/* unit vector along established bit */
    ,newdir[4], 			/* unit vector along candidate con */
    phi,				/* angle the seam makes here */
    mindotp, dotp[2];			/*  */
} br_t, *br_p;

/*** global numbers, data ***/
/* in main.args: (options) */
extern int opt_careful, opt_maximize, opt_warn, opt_connect_touch, debug,
	opt_general;
/* extern atom_p root_of_atoms;		/# for calculating index #/
 * extern shell_p root_of_shells;
 */

extern FILE * connections_file;

typedef struct radius_param_ {
  double 
/*    radius, OBSOLETE */
    extension,
    reduction,
    polar_extension,
    polar_reduction;
} radius_param_t;

extern radius_param_t radius_param;

/* #define opt_radius          (radius_param.radius) OBSOLETE */
#define opt_extension 	    (radius_param.extension)
#define opt_reduction 	    (radius_param.reduction)
#define opt_polar_extension (radius_param.polar_extension)
#define opt_polar_reduction (radius_param.polar_reduction)

/* in misc.c: */
extern FILE  *logfile;

extern int overall_npatches;
extern int maxnpoints;			/* set in choose_(uni)_sphere */
extern int ** fast_mod;			/* matrix for speeding up modulo-op. */

/* in templates.c: */
extern sphere_p Template_spheres[ N_TMPLT_SPHERES ];
extern int Template_npoints[ N_TMPLT_SPHERES ];

#define TEMPLATE_NPOINTS_CONTENTS  { 12, 32, 42, 92, 122, 162, 252  }
#define TEMPLATE_NPOINTS_HALFWAY     { 22, 37, 64, 102, 142, 204 }

extern int Template_based_on[ N_TMPLT_SPHERES ];
#define TEMPLATE_BASED_ON_CONTENTS { \
  SPH_ICOS, SPH_PDOD, SPH_ICOS, SPH_ICOS, SPH_PDOD, SPH_ICOS, SPH_ICOS  }

extern int Template_order[ N_TMPLT_SPHERES ];
#define TEMPLATE_ORDER_CONTENTS  { 1,  1,  2,  3,   2,   4,   5  }

/*** prototypes: ***/
/* in io.c: */
box_p get_atoms(FILE * pointsfile, char* protein, char*binspec, uint pdbflags);
void get_template_spheres(FILE *file);

void print_atom(FILE*, atom_p, char*atomcolor, char*conn_color);
void write_atom(char*filename, atom_p, char*atomcolor, char*conn_color);
void print_cons(FILE*, atom_p, char *color);
void write_cons(char*filename, atom_p, char *color);
void print_sphere(FILE*, atom_p, char *color);
void write_sphere(char *filename, atom_p, char *color);
void print_triangulines(FILE*, box_p);

int print_atom_dots(FILE*, atom_p, char*prepend, char*append);
int print_dots(FILE*, box_p box);
void print_areas(FILE*, box_p box);

/* in templates.c: */
void generate_template_spheres(void);
void randomize_points(sphere_p Sphere);

sphere_p tesselate(int type, int order);
void triangulate_sphere(int opt_c, sphere_p);

sphere_p recursion_tessel(int type, int freq);
void free_templates(void);

/* in acc.c: */
int accessability(box_p box);
void free_accessability(box_p box);


/* in top.c: */
int atom_topology(box_p, int save_internal_points);
void clear_top_flags(shell_p);
int find_patches(shell_p, int opt_recover, int *nbndriesp, int firsttime);
int delete_bndry(bndry_p);
int remove_widows(shell_p shell, int *ndeletedp);
int free_bndry(bndry_p);
int free_patch(patch_p patch);
int free_patches(patch_p patches);
void free_topology(box_p);
void free_neibordata(shell_p);
void free_all_neibordata(box_p);

/* in conn.c: */
con_p alloc_con(const con_t *oldcon);
void free_con(con_p);

/* con_p save_con(bndry_p a, int offsa, bndry_p b, int offsb,
 * 	       double *dir, double d2, double rb);
 */
/* int remove_bndry_con(bndry_p); */
/* int delete_conref(con_p con, int end); *//* one side; just the refenrence*/

int remove_con(con_p con);		/* both */
int remove_con_halfknown(list_p, bndry_p); /* if one side known  */

list_p * con_place(con_p toadd, bndry_p bndry, int offset);
int connections(box_p box, int mode);
int widow_check(shell_p, int mode);
void free_connections(box_p box);

/* in seam.c: */
const float * bndry_coordinates(bndry_p bndry);
int find_seam(con_p firstcon, int firstend, int * nseamsp, 
	      int maxnseams, seam_p seams);
/* void set_seams(int nseams, seam_p seams); */

/* int do_pair(seam_p seams, int sort); */
/* int do_simple(shell_p shell); */

/* in stitch.c: */
int stitch(box_p);

/* in misc.c: */
/* int ** fast_modulo_matrix(int from, int to); */

sphere_p choose_uni_sphere(int pts_per_atom);
sphere_p choose_sphere(atom_p atomp, double density);
const char * atom_spec(box_p box, int i);

/* in surf.c: */
int print_faces(FILE *file, box_p box);
int print_surfaces(FILE*, surf_p Surfaces, const char *nature);
surf_p get_surfaces(box_p box, const char * nature);
  /* returns a list of 'surfaces', each of which is a contiguous piece of */
  /* surface of specified nature. For easy access, a surf */
int surftop(surf_p surf);
void free_surface(surf_p);
int recover(FILE*, box_p, int patdots);

/* in randomize.c: */
int ran_atoms(box_p box, int percentage);

#endif /* _QUILT_H_ */
