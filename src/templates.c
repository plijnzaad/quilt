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
 * Philip Lijnzaad, plijnzaad@gmail.com			      *
 *								      *
 * ------------------------------------------------------------------ */

/* support routines for triangulation of one sphere of radius 1 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "extra-math.h"

#include "utils.h"
#include "alloc.h"
/* #include "files.h" */
#include "vecmat.h"
#include "hash.h"

#include "quilt.h"

/* mask use of Floating Point Exception catcher: */
#define catch_fpe(X)

/* program global variables */
sphere_p Template_spheres[ N_TMPLT_SPHERES ];

int Template_npoints[ N_TMPLT_SPHERES ] = TEMPLATE_NPOINTS_CONTENTS ;
int Template_based_on[ N_TMPLT_SPHERES ] = TEMPLATE_BASED_ON_CONTENTS ;
int Template_order[ N_TMPLT_SPHERES ] = TEMPLATE_ORDER_CONTENTS ;


/* file global variables */
static float ** Distances2;		/* upper diagonal matrix */
static double cutoff, cutoff2, arccutoff;
static double cutoff_val[4]= { 1.1, 1.2, 2.0, 1.5 }; 
/* for regular icosahedral, pentakisdodecahedral, Zielenkiewicz, unknown */
#define NEDGES(N) (3*(N) - 6)
#define NFACES(N) (2*(N) - 4)

#define ICOS_NPOINTS 12			/* icoshedral points */
#define ICOS_NEDGES  NEDGES(ICOS_NPOINTS) /* 30 */
#define ICOS_NFACES  NFACES(ICOS_NPOINTS) /* 20 */

#define PDOD_NPOINTS 32			/* pentakisdodecahedral points */
#define PDOD_NEDGES  NEDGES(PDOD_NPOINTS) /* 90 */
#define PDOD_NFACES  NFACES(PDOD_NPOINTS) /* 60 */

/* icosahedral points according to 
   (+-1, 0, +-t)
   (+-t, +-1, 0)
   (0, +-t, +-1)
   where t is golden ratio; scaled down so as to make radius 1.
*/
#define tS 0.5257311
#define tT 0.8506508

static  float icos_points[ICOS_NPOINTS][3]= { /* coordinates */
    {  tS, 0.0,  tT, },
    {  tS, 0.0, -tT, },
    { -tS, 0.0,  tT, },
    { -tS, 0.0, -tT, },
    
    {  tT,  tS, 0.0, },
    { -tT,  tS, 0.0, },
    {  tT, -tS, 0.0, },
    { -tT, -tS, 0.0, },
    
    { 0.0,  tT,  tS, },
    { 0.0,  tT, -tS, },
    { 0.0, -tT,  tS, },
    { 0.0, -tT, -tS, } };

static int icos_edges[ICOS_NEDGES][2]= { /* adjacent pairs of above points */
    { 3, 7}, { 1, 4}, { 7, 11}, { 0, 4}, { 7, 10}, { 0, 6}, { 6, 11}, 
    { 0, 8}, { 6, 10}, { 0, 10}, { 5, 9}, { 5, 8}, { 4, 9}, { 3, 9}, { 1, 11}, 
    { 1, 6}, { 4, 8}, { 3, 11}, { 1, 9}, { 3, 5}, 
    { 2, 5}, { 2, 7}, { 2, 10}, { 2, 8}, { 0, 2}, { 1, 3}, { 10, 11}, { 5, 7},
    { 8, 9}, { 4, 6 } }; 

static int icos_faces[ICOS_NFACES][3]={ /* adjacent triples of above edges */
    { 0, 6, 10}, { 0, 2, 10}, { 0, 2, 8}, { 0, 4, 8}, { 0, 4, 6}, { 1, 3, 9}, 
    { 1, 3, 11}, { 1, 6, 11}, { 1, 4, 6}, { 1, 4, 9}, { 2, 7, 10}, { 2, 5, 7},
    { 2, 5, 8}, { 3, 5, 9},  { 3, 5, 7}, { 3, 7, 11}, { 4, 8, 9}, { 5, 8, 9}, 
    { 6, 10, 11}, { 7, 10, 11 } }; 

/* and now the pentakisdodecahedron (1st tesselation of the dodecahedron)  */
static float pdod_points[PDOD_NPOINTS][3] = { /* coordinates */ 
  { 0.5257311, 0., 0.8506508},
  { 0.5257311, 0., -0.8506508},
  {-0.5257311, 0., 0.8506508},
  {-0.5257311, 0., -0.8506508},
  { 0.8506508, 0.5257311, 0.},
  {-0.8506508, 0.5257311, 0.},
  { 0.8506508, -0.5257311, 0.},
  {-0.8506508, -0.5257311, 0.},
  { 0., 0.8506508, 0.5257311},
  { 0., 0.8506508, -0.5257311},
  { 0., -0.8506508, 0.5257311},
  { 0., -0.8506508, -0.5257311},
  { 0.577350268, -0.577350268, 0.577350268},
  { 0., -0.356822091, 0.934172359},
  { 0., 0.356822091, 0.934172359},
  { 0.577350268, 0.577350268, 0.577350268},
  { 0.934172359, 0., 0.356822091},
  { 0., 0.356822091, -0.934172359},
  { 0., -0.356822091, -0.934172359},
  { 0.577350268, -0.577350268, -0.577350268},
  { 0.934172359, 0., -0.356822091},
  { 0.577350268, 0.577350268, -0.577350268},
  {-0.577350268, -0.577350268, 0.577350268},
  {-0.934172359, 0., 0.356822091},
  {-0.577350268, 0.577350268, 0.577350268},
  {-0.577350268, 0.577350268, -0.577350268},
  {-0.934172359, 0., -0.356822091},
  {-0.577350268, -0.577350268, -0.577350268},
  { 0.356822091, 0.934172359, 0.},
  {-0.356822091, 0.934172359, 0.},
  { 0.356822091, -0.934172359, 0.},
  {-0.356822091, -0.934172359, 0.}
};

static int pdod_edges [PDOD_NEDGES][2] = {
  {9, 17}, {5, 29}, {6, 16}, {5, 26}, {5, 23}, {1, 17}, { 0, 13}, { 0, 14}, 
  {1, 20}, { 0, 16}, {6, 20}, {6, 30}, {7, 23}, {7, 26}, {4, 28}, {4, 20}, 
  {7, 31}, {8, 14}, {4, 16}, {8, 28}, {8, 29}, {9, 28}, {1, 18}, {10, 30}, 
  {10, 31}, {11, 30}, {3, 17}, {3, 26}, {9, 29}, {10, 13}, {3, 18}, 
  {11, 31}, {11, 18}, {2, 23}, {2, 13}, {2, 14}, { 0, 12}, {1, 19}, {2, 22}, 
  { 0, 15}, {2, 24}, {1, 21}, {11, 27}, {11, 19}, {4, 15}, {8, 15}, {5, 24}, 
  {10, 22}, {8, 24}, {6, 19}, {9, 21}, {10, 12}, {4, 21}, {3, 25}, {7, 22}, 
  {3, 27}, {5, 25}, {9, 25}, {7, 27}, {6, 12}, {14, 24}, {17, 21}, 
  {15, 16}, {15, 28}, {18, 19}, {19, 20}, {19, 30}, {20, 21}, {17, 25}, 
  {21, 28}, {13, 22}, {22, 23}, {22, 31}, {14, 15}, {23, 24}, {24, 29}, 
  {12, 16}, {25, 26}, {12, 30}, {12, 13}, {25, 29}, {26, 27}, {18, 27}, 
  {27, 31}, {30, 31}, {23, 26}, {17, 18}, {28, 29}, {13, 14}, {16, 20}
};

static int pdod_faces[PDOD_NFACES][3] ={
  { 0, 14, 15}, { 0, 15, 16}, { 0, 16, 12}, { 0, 12, 13}, { 0, 13, 14}, 
  {1, 18, 19}, {1, 19, 20}, {1, 20, 21}, {1, 21, 17}, {1, 17, 18}, 
  {2, 24, 14}, {2, 14, 13}, {2, 13, 22}, {2, 22, 23}, {2, 23, 24}, 
  {3, 17, 25}, {3, 25, 26}, {3, 26, 27}, {3, 27, 18}, {3, 18, 17}, 
  {4, 21, 20}, {4, 20, 16}, {4, 16, 15}, {4, 15, 28}, {4, 28, 21}, 
  {5, 24, 23}, {5, 23, 26}, {5, 26, 25}, {5, 25, 29}, {5, 29, 24}, 
  {6, 20, 19}, {6, 19, 30}, {6, 30, 12}, {6, 12, 16}, {6, 16, 20}, 
  {7, 23, 22}, {7, 22, 31}, {7, 31, 27}, {7, 27, 26}, {7, 26, 23}, 
  {8, 24, 29}, {8, 29, 28}, {8, 28, 15}, {8, 15, 14}, {8, 14, 24}, 
  {9, 21, 28}, {9, 28, 29}, {9, 29, 25}, {9, 25, 17}, {9, 17, 21}, 
  {10, 31, 22}, {10, 22, 13}, {10, 13, 12}, {10, 12, 30}, {10, 30, 31}, 
  {11, 30, 19}, {11, 19, 18}, {11, 18, 27}, {11, 27, 31}, {11, 31, 30}
};

/* macros: */
#define CNT(X)				/* mask definition */

/* constants: */
#define DFLT_DIM 4			/* unit sphere grid dimension */
#define DFLT_SPACING (2.0/DFLT_DIM)

#define MIDPOINTS_TOO_CLOSE_TOL 0.01	/* if edge midpoints closer, discard */
#define ANGLE_TOL 0.1			/* tolerance for grid mismatch */
#define AREA_TOL  0.02			/* area of triangles should add up */
					/* within this tolereance */

#define BOXING_TOL 0.02			/* needed during hashing */

#define COLLINEAR_TOL 0.01		/* triangle collinear */

/* prototypes: */
/* global functions: */
sphere_p tesselate(int type, int order); /* global for export to sph3.c */
sphere_p recursion_tessel(int type, int freq); /* likewise */

void triangulate_sphere(int opt_c, sphere_p); /* 'main' routine */
void generate_template_spheres(void);
void randomize_points(sphere_p Sphere);
void triangulate_sphere(int opt_c, sphere_p);
void free_templates(void);

/* local functions: */

/* generating points: */

static void subdiv_edges(sphere_p Sphere, float *a, float *b);
static void subdiv_faces(sphere_p Sphere, float *a, float *b, float *c);

/* static void convert_to_final(sphere_p Sphere); */

/*** 'old' (makes use of global hash table; for recursive tesselation) ***/
static void tesselate_face(int freq, float * a, float *b, float *c);
float * subdiv_triangle(int freq, float * a, float *b, float *c);
float * interior_point(float * a, float * b, float *c, 
			   float u, float v);
static void store_point(float*xyz);
static int purge(htable_p point_hash, float edge_length);
/*** end of old stuff ***/


static int edge_crossing(sphere_p Sphere, edge_p a, edge_p b);
static void distances(sphere_p);

static void find_edges(int opt_c, sphere_p);
static int edge_cmp(const void *, const void *);

static void find_faces(sphere_p);
static void areas(sphere_p);
double lHuilier(double a, double b, double c); /* make static ? */

static void edges_misc(sphere_p);
static void set_whoamis(sphere_p);
static void sort_rays(sphere_p, int);
static int ray_cmp(const void *, const void *);

static int make_Delaunay(sphere_p Sphere);
static int spherical_delaunay(float * a, float *b, float *c, float *d, 
			      int *colp);
static double spherical_angle(double a, double b, double c, int*colp);
static int replace_edge(sphere_p Sphere, int e, int a, int b, int c, int d);

static int replace_edge(sphere_p Sphere, int e, int a, int b, int c, int d);

/* static void boxify_points(sphere_p sphere); */

static double estimate_arclen(int npoints);
void randomize_points(sphere_p);	/* also in *.h */

static void points_misc(sphere_p Sphere) { 
  /*
   * does    : sets some extra info on points, mainly more extensive info
   *   on neigbours
   * gets    : 
   *
   * affects : 
   *
   * returns : 
   *
   * warns   : 
   *
   * comment : used to be convert_to_final
   */


  point_p Points, Pp;
  edge_p Edges=Sphere->Edges, Ep;
  int i, j, npoints, nngbs, o, q, e;
  float * coordinates;

/*  Points=(point_p)CALLOC(npoints, sizeof(point_t)); */
  Points=Sphere->Points;
  npoints=Sphere->npoints;
  coordinates=Sphere->coordinates;
  Pp=Points;

  for (i=0; i<npoints; i++, Pp++) {
    nngbs=Pp->nngbs;
    for (j=0; j<nngbs; j++) {
      e= Pp->edges[j];			/* edge index */
      Ep= Edges+e;
      o= (Ep->ends[0]==i);		/* orientation of this edge */
      q= Ep->ends[o];			/* the other point */
      Pp->ngbs[j]=q;
      Pp->Ngbs[j]=Points+q;
      Pp->ngbs_xyz[j]= Points[q].xyz;
      Pp->whoamis[j]=Ep->whoami[o];
    }
  }
  return;
} /* points_misc */

void triangulate_sphere(int opt_c, sphere_p Sphere) {
  /*
   * does    : Delaunay-triangulates sphere->Points to sphere->Edges and
   *   sphere->Faces
   *
   * gets    : 
   *
   * affects : 
   *
   * returns : a triangulated sphere 
   *
   * warns   : about Delaunay-hood of sphere, and when it can't make it.
   *
   * comment : 
   */
  int i, npoints=Sphere->npoints;
  double arclen;

  catch_fpe(FPE_DEFAULT);		/* doubleing point error handling */
  distances(Sphere);			/* squared distances, actually */

  if(Sphere->ideal_arclen==0.0) {	/* may have be done before, in  */
    arclen=estimate_arclen(npoints);	/* randomize_points */
    Sphere->ideal_arclen=arclen;
    Sphere->edgelengths[2]=2.0*sin(0.5*arclen);
  }

  arclen=Sphere->ideal_arclen;
  arccutoff = cutoff_val[Sphere->type]*arclen;
  /* this is a relic from the time  when it also could be an irregular */
  /* distribution. Could dispense with it */

  cutoff= cutoff_val[Sphere->type]*2*sin( 0.5*arclen); 

  cutoff2= SQR(cutoff);

  catch_fpe(FPE_DEFAULT | FPE_UNDERFLOW ); /* doubleing point error handling */

  find_edges(opt_c, Sphere); 
  FREE(Distances2); 

  if ( Sphere->nedges != (3*npoints -6))
    DIE("\
couldn't accomplish triangulation; points too irregular? try running\n\
with '-c' (takes longer), or increasing MAX_NGBS_PER_PNT (currently %d),",
	  MAX_NGBS_PER_PNT); 

  for (i=0; i<npoints; i++)		/* sort edges, per point */
    sort_rays(Sphere, i);

  set_whoamis(Sphere);			/* topology info on the edges */

  if (Sphere->type >= SPH_ZIEL &&	/* it's an irregular distribution */
      Sphere->npoints > 8) {		/* otherwise non-sensical answers */
    if((i=make_Delaunay(Sphere))) {
      warn("%d edges re-assigned to get Delaunay triangulation\n", i);
    } else { 
      warn("first try already a Delaunay triangulation\n");
    }
  }

  find_faces(Sphere);
  assert(Sphere->nfaces == 2*npoints-4);

  areas(Sphere);

  edges_misc(Sphere);			/* directions, lengths, statistics */
  set_whoamis(Sphere);			/* again, destroyed by find_faces */

  points_misc(Sphere);		/* convert t_point_t to point_t */
  return;
} /* triangulate_sphere */


static void edges_misc(sphere_p Sphere) {
  /* sets edgelengths, and some numbers needed during searching over the */
  /* sphere, specifying the min/max fineness of the mesh */
  int i, k, a, b; 
  edge_p ep;
  double f, g ,minlen=100.0, maxlen= -1.0, minarclen=100.0, maxarclen= -1.0
    , maxangle= -100.0, maxgamma, direction[3];
  
  point_p Points= Sphere->Points;
  edge_p Edges=Sphere->Edges; 
  int nedges=Sphere->nedges;

/*   face_p Faces=Sphere->Faces; don't need yet
 *   int nfaces=Sphere->nfaces;
 */

  for (i=0; i<nedges; i++) {
    ep=Edges+i;
    a=ep->ends[0];
    b=ep->ends[1];
    VEC3MIN(direction, Points[b].xyz, Points[a].xyz); /* direction */

    f=ep->length2;
    f=SQRT(f); 
    ep->length2=f;			/* now an unsquared length!  */

    for (k=0; k<3; k++) 
      ep->direction[k]=direction[k]/f;	/* normalized direction vector */

    g=2.0*ASIN(f*0.5);
    ep->arclen=g;			/* arclength */
    
    if (f < minlen) {
      minlen=f; 
/*      assert(g < minarclen); */
      minarclen=g; 
    }
    if (f > maxlen) { 
      maxlen=f; 
/*      assert(g > maxarclen); */
      maxarclen=g; 
    }
  } /* for i < nedges */
  
  /* calculate min,max of internal angles assuming an equilateral triangle: */ 

  /* minangle=2.0*atan(SQRT(sin(0.5*minarclen)/sin(1.5*minarclen))); */
					/* NOT NEEDED */

  maxangle=2.0*atan(SQRT(sin(0.5*maxarclen)/sin(1.5*maxarclen)));

  /* from SQR(tan(A/2) = sin(s-b)*sin(s-c)/(sin(s)*sin(s-a)), where A */
  /* (etc.) is angle on surface of sphere, a etc is (arc)length of */
  /* an elliptical triangle's side , and s is (a+b+c)/2 */


  /* calculate maximal radius (as (arc)length of an elliptical triangle's */
  /* side) that might not contain points, from tan c = tan a / cos B, */
  /* where c is the hypothenuse of the rectangular elliptical triangle, */
  /* and B is 0.5* above angle  */

  /* mingamma = atan(tan(0.5*minarclen)/cos(0.5*minangle)); */
					/* NOT NEEDED */
  maxgamma = atan(tan(0.5*maxarclen)/cos(0.5*maxangle));

  /* transfer relevant bits to structure: */

  Sphere->edgelengths[0]=minlen;
  Sphere->edgelengths[1]=maxlen;
  Sphere->arclengths[0]=minarclen;
  Sphere->arclengths[1]=maxarclen;
/* ideal values of Sphere->edgelengths[2] and Sphere->arclengths[2] */
/* are already set in triangulate_sphere */ 

  f=Sphere->arclengths[2];		/* ideal value */
  Sphere->halfalpha=0.5*f;
/* was:  Sphere->halfalpha=0.5*maxarclen; to maximize halfalpha */

  Sphere->beta = ACOS (cos(f)/cos(0.5*f) );   /* called 'gamma' in logbook */
  Sphere->sinbeta= sin(Sphere->beta);
/* was: Sphere->beta = cos(minarclen)/cos(0.5*maxarclen);, 
 * but should have been:
 * Sphere->beta = ACOS  !!! ( cos(maxarclen !!)/cos(0.5*minarclen !! ) );
 * (so as to 'maximize' beta)
 */
  return;
} /* edges_misc */


static void set_whoamis(sphere_p Sphere) {
  /* for every edge, puts the index a particular edge has in the ray of */
  /* edges emanating out of a point, in edge.sw.w[0,1]. This is the index */
  /* in the array point.edges[] */
  int i, a;
  EDGE_ID *edge_array;
  uchar j;
  edge_p  ep; 
  point_p pp;

  point_p Points= Sphere->Points;
  int npoints=Sphere->npoints;
  edge_p Edges=Sphere->Edges; 

  for (i=0; i<npoints; i++) {
    pp=Points+i;
    edge_array=pp->edges;
    for (j=0; j<pp->nngbs; j++) {
      ep=Edges + edge_array[j];
      a=ep->ends[0]; /* b=ep->ends[1]; */
      if (i==a) 
	ep->whoami[0]=j;
      else
	ep->whoami[1]=j;
    }
  }
} /* set_whoamis */

#ifdef PRELIM				/* debugging stuff */
  FILE * pa, *pb, *pc;
#  define OPEN_PRELIM \
pa=fopen("pa","w"); pb=fopen("pb", "w"); pc=fopen("pc", "w")
#  define PRINT_PRELIM(F) fprintf(F,"%10.4f %10.4f %10.4f\n",p[0],p[1],p[2])
#  define CLOSE_PRELIM fclose(pa); fclose(pb); fclose(pc)
#else
#  define OPEN_PRELIM			/* void */
#  define PRINT_PRELIM(F)		/* void */
#  define CLOSE_PRELIM			/* void */
#endif

sphere_p tesselate(int type, int freq) {
  /*
   * does    : generates a spherical tesselation by division into freq^2 of 
   *   the faces of an icosahedron (type==SPH_ICOS) or pentakisdodecahedron
   *   (type== SPH_PDOD)
   *
   * gets    : frequency
   *
   * affects : nothing
   *
   * returns : a sphere, whose Points is loaded with tPoints
   *
   * warns   : 
   *
   * comment : 
   */

  sphere_p Sphere;
  int  n,i, n_base_points, n_base_edges, n_base_faces;
  int *base_edges, *base_faces;
  float *base_points, *a, *b, *c, *coordinates;

  OPEN_PRELIM;

  Sphere=(sphere_p)CALLOC(1, sizeof(sphere_t)); 
  Sphere->type= type;
  Sphere->order=freq;

  switch(type) {
  case SPH_ICOS:   
    Sphere->maxringlen = 8*freq;
    n = 10 * SQR(freq) + 2;		/* nr of points exspected */ 
    n_base_points = ICOS_NPOINTS;
    n_base_edges = ICOS_NEDGES;
    n_base_faces = ICOS_NFACES;
    base_points = (float*)icos_points;
    base_edges = (int*) icos_edges;
    base_faces  = (int *)icos_faces;
    break;
  case SPH_PDOD:
    Sphere->maxringlen = 13*freq;
    n = 30 * SQR(freq) + 2;		/* nr of points exspected */
    n_base_points= PDOD_NPOINTS;
    n_base_edges = PDOD_NEDGES;
    n_base_faces = PDOD_NFACES;
    base_points = (float*)pdod_points;
    base_edges = (int*) pdod_edges;
    base_faces  = (int *)pdod_faces;
    break;
  default: die("unknown tesselation type: %d\n", type);
    break;
  }

  if (n > MAXNROFPOINTS)
    DIE("\
types SHORT/EDGE_ID don't allow more than %d points; trying to generate %d points,",
	      MAXNROFPOINTS,n);

  Sphere->coordinates=coordinates=(float*)CALLOC(n, sizeof(float[3]));
/*  Points=CALLOC(n, sizeof(t_point_t)); */
/*  Sphere->Points=(point_p)Points; */

  for (i=0; i<n_base_points; i++)		/* transfer base points */
    VEC3ASS(coordinates+3*i, base_points+3*i);

  Sphere->npoints=n_base_points;	/* to start with (may increase) */

  if ( freq == 1 )			/* no subdivision: we're done */
    return Sphere;

  for (i=0; i< n_base_edges; i++) {	/* fill in new edge points */
    a=base_points + 3*base_edges[i * 2];
    b=base_points + 3*base_edges[i * 2 + 1];
    subdiv_edges(Sphere, a, b);
  }

  if ( freq == 2 )
    return Sphere;			/* since doesn't have new facepoints */

  for (i=0; i< n_base_faces; i++) {	/* new faces points */
    a=base_points + 3* base_faces[i*3] ;
    b=base_points + 3* base_faces[i*3+1];
    c=base_points + 3* base_faces[i*3+2];
    subdiv_faces(Sphere, a,b,c);
  }

  assert(Sphere->npoints == n);

/*   for (i=0; i<n; i++) 
 *     Points[i].xyz=coordinates+3*i;
 */

  CLOSE_PRELIM;
  return Sphere;
} /* tesselation */

static void subdiv_edges(sphere_p Sphere, float *a, float *b) {
  int npoints=Sphere->npoints;
  int freq=Sphere->order;
  float * coordinates=Sphere->coordinates;
/*  t_point_p Points= Sphere->Points; */
  int i;
  double p[3], u[3], v[3], angle, mat[3][3];


  if (freq==2) {			/* very simple: just midle */
    for (i=0; i<3; i++)
      u[i]=0.5*(a[i]+b[i]);
    NRMIZE3(u);
    VEC3ASS(coordinates+3*npoints, u);
    Sphere->npoints++;
    return;
  }

  /* following may not be the right way to do it, if the faces are */
  /* subdivided from a 'flat' grid on the original triangle */
  angle=ANGLE(a,b);			/* now the more difficult case */
  angle /= freq; 
  CROSS3P(p, a, b);			/* the vector about which to rotate */
  NRMIZE3(p);
  ARBANGROTMAT(mat, p, angle);		/* calculate rotation matrix */

  VEC3ASS(u, a);			/* and now apply the matrix to x */
  for (i=1; i<freq; i++) {		/* skip first one! */
    VEC3MAT(v,u, mat);			/* y now new point, lying on edge */
    VEC3ASS(coordinates+3*npoints, v);
    VEC3ASS(u, v);			/* for next round */
    npoints++;
  }
  Sphere->npoints=npoints;
  return;
} /* subdiv_edges */

static void subdiv_faces(sphere_p Sphere, float *a, float *b, float *c) {
  /*
   * does    : directly subdivide one triangle into freq^2 smaller triangles
   *
   * gets    : edge points of triangle
   *
   * affects : contents of Sphere->Points, and Sphere->npoints
   *
   * returns : 
   *
   * warns   : 
   *
   * comment : for triangle abc, uses i,j as 'coordinates' during
   *   subdivision. For the triangles bca and cab, first maps i, and j to
   *   their own coordinates u and v, so as to keep same order.
   */

  int freq=Sphere->order ,i,j,u,v,n;
/*  t_point_p Points= Sphere->Points; */
  int offset, npoints;
  float *where, *p, *coordinates;
  double m[3], rf, uf, vf;
  
  coordinates=Sphere->coordinates;

  if (freq==2) 
    return;				/* no interior points */

  npoints=Sphere->npoints;

  if (freq==3) {			/* one interior point: just midle */
    for (i=0; i<3; i++)
      m[i]=(a[i]+b[i]+c[i]);
    NRMIZE3(m);
    VEC3ASS(coordinates+3*npoints, m);
    Sphere->npoints++;
    return;
  }

#define FILL_IN(I,J, A,B,C) u=I; v=J;					\
      			    if (u!=0 && u != freq) {			\
                              uf= u*rf; vf= v*rf;			\
			      p = interior_point(A, B, C, uf, vf);	\
			      VEC3INC(where, p);	 		\
			    }

  offset=Sphere->npoints;		/* nr of points already done */
#define INDEX(I,J) (((I-2)*(I-1))/2 + J - 1)

  rf = 1.0/(double)freq;
  for (i=2; i<freq; i++)
    for (j=1; j<i; j++) {
      n=offset+INDEX(i,j);		/* which point are we working on */
      where = coordinates+3*n;		/* where to put result */


      /* ab = x-axis, ac = y-axis: abc triangle */
      FILL_IN(i, j, a, b, c);
      PRINT_PRELIM(pa);
      
      /* bc = x-axis, ba = y-axis: bca triangle */
      FILL_IN(freq-i+j, freq-i, b, c, a);
      PRINT_PRELIM(pb);

      /* ca = x-axis, cb = y-axis: cab triangle */
      FILL_IN( freq-j, i-j, c, a, b);
      PRINT_PRELIM(pc);
    }

  /* normalize : */
  n = ((freq-1)*(freq-2))/2;		/* the nr to be generated here */
  for (i=0; i<n; i++) {
    where = coordinates + 3*(offset + i);
    NRMIZE3(where);
  }
  Sphere->npoints += n;
  return; 

#undef INDEX
#undef FILL_IN

} /* subdiv_faces  */

void randomize_points(sphere_p Sphere) { 
  /*
   * does    : rotates points such that all points->xyz[i] lie within <-1, +1>
   *
   * gets    : sphere
   *
   * affects : xyz's
   *
   * returns : nothing
   *
   * warns   : 
   *
   * comment : does it by rotating about vector towards 111, over an angle
   *   of 1/3 * sphere->ideal_arclen 
   */

  int i;
  int npoints=Sphere->npoints;
/*  t_point_p Points=Sphere->Points; */
  float * coordinates;
  double p[3], v[3], angle, mat[3][3], arclen;

  arclen=estimate_arclen(npoints);
  Sphere->ideal_arclen=arclen;
  Sphere->edgelengths[2]=2.0*sin(0.5*arclen);
  coordinates=Sphere->coordinates;

  angle= 0.3333*Sphere->ideal_arclen;

  p[0]=p[1]=p[2]=M_1SQRT3;		/* rotate about axis pointing to 111 */

  ARBANGROTMAT(mat,p, angle);

  for (i=0; i<npoints; i++) { 
    VEC3MAT(v, coordinates+3*i, mat);
    VEC3ASS(coordinates+3*i, v);
  }
  return;
} /* randomize_points */

static void distances(sphere_p Sphere) {
  int i,j;
  float *fp, *coordinates;
/*  t_point_p Points=Sphere->Points; */
  int npoints=Sphere->npoints;

  coordinates=Sphere->coordinates;
  Distances2=(float**)CALLOC( (npoints-1),  sizeof(float*) );
  fp=(float*)CALLOC( ((npoints * (npoints-1) )/2), sizeof(float));

  for (i=0; i<npoints-1; i++) {
    Distances2[i]= fp			/* beginning of this row */
      - (i+1);				/* and make look like we start at 0 */
    fp += (npoints -i -1);		/* advance to next row */
  }

  for(i=0; i<npoints; i++)
    for(j=i+1; j<npoints; j++)
      Distances2[i][j] = DIST3_2(coordinates+3*i, coordinates+3*j);

/*    was:   Distances2[i][j] = DIST3_2(Points[i].xyz,Points[j].xyz);
 */

  return;
} /* distances */

static void find_edges(int careful, sphere_p Sphere) {
  /*
   * does    : allocate Sphere->Edges and fills them, given Sphere->Points
   *
   * gets    : careful, a flag, and Sphere that has Points already initialized
   *
   * affects : 
   *
   * returns : 
   *
   * warns   : 
   *
   * comment : the careful flag is in fact not needed, as we're usually
   *   dealing with 'regular' polyhedra. I messed with it a little; may
   *   not work anymore for irregular polyhedra (when it still worked, it
   *   was not impressive anyway ...)
   */

  SHORT  *sp;
  /* NOTE: some of these things (like ncand,i,j) index into arrays */
  /* potentially longer than TYPE_MAX, so should not be of type SHORT ! */
  int i,oldi=0, j, ncand, l,h, maxspread, ncross;
  int maxn;
  float d2;
  edge_p cand_edges, candp, *candp_array, edgei, edgej;
  int nngbs_histo[MAX_NGBS_PER_PNT+1]; 
  
  int current_max;

  point_p Points=Sphere->Points;
  int npoints=Sphere->npoints;
  edge_p Edges=Sphere->Edges; 
  int nedges=Sphere->nedges;
/*   face_p Faces=Sphere->Faces; 
 *   int nfaces=Sphere->nfaces;
 */  

  if (Sphere->type==SPH_ZIEL) { 
    maxn= 10*npoints;			/* was: npoints*(npoints-1)/4 ... */
    if (maxn <=200 )			/* to avoid trouble with the lower */
      maxn=200;				/* triangulations */
  }
  else 
    maxn=5*npoints;			/* wild guess */

  current_max=MAX_NGBS_PER_PNT;		/* to start with */
  maxspread= 2 ;
  /* initialize 'histogram' keep track of balance: */
  nngbs_histo[0]=npoints;		/* all points have 0 neighbours */
  for (i=1; i<=MAX_NGBS_PER_PNT; i++)
    nngbs_histo[i]=0;


#ifndef NO_ALLOCA
  cand_edges=(edge_p)alloca(i=maxn*sizeof(edge_t)); /* the data */
  bzero((char*)cand_edges, i);
  candp_array=(edge_p*)alloca(i=maxn*sizeof(edge_p)); /* pointers to it */
  bzero((char*)candp_array, i);
#else
  cand_edges=(edge_p)CALLOC(maxn, sizeof(edge_t)); /* the data */
  candp_array=(edge_p*)CALLOC(maxn, sizeof(edge_p)); /* pointers to it */
#endif

  Edges=(edge_p)CALLOC(3*npoints-6, sizeof(edge_t)); /* final list of edges */

  ncand=0;
  candp=cand_edges;
  for (i=0; i<npoints; i++)		/* construct array of candidates */
    for (j=i+1; j<npoints; j++) {
      if ( (d2=Distances2[i][j]) < cutoff2) {
	CNT(INSIDE);
	candp->ends[0]=(SHORT)i;
	candp->ends[1]=(SHORT)j;
	candp->length2=d2;
	candp_array[ ncand ]= candp;
	candp++; 
	if( ++ncand > maxn)
	  die("\
finding too many candidate points; decrease cut_off_val[] for distribution\n");
      } else { CNT(OUTSIDE); }
    }
  assert(ncand <= maxn);

  if ( ncand < (3*npoints -6))
    die("\
couldn't accomplish triangulation; found too few candidate points within\n\
cutoff; try increasing cutoff_val[type] (relating used cutoff to that needed\n\
for an equilateral triangle; currently %f). Can also be caused by overflow\n\
of the index; in that case make type bigger\n", cutoff_val[Sphere->type]);
  
  qsort(candp_array, ncand, sizeof(edge_p), edge_cmp); /* sort this */


  for(i=0; i< ncand; i++) {		/* and now triangulate, starting */
					/* with shortest edges */

    if(careful) {
      /* Try to do it carefully */
      /* The idea of this is that the coordination number, the number of */
      /* neigbours per point, should be not too skew. When it is, */
      /* 'current_max', the number of neighbours per point, is lowered, so */
      /* that new candidate edges are discarded if they would exceed it. */
      /* When the skewness disappears, the discarded candidate edges are */
      /* tried again, if they haven't been accepted */

      /* determine spread in the histogram: lowest nr of neighbours found: */
      for (j=0; j<MAX_NGBS_PER_PNT -1 ; j++)
	if (nngbs_histo[j] ==0 && nngbs_histo[j+1]!=0) 
	  l=j;
      /* and highest */
      for (j=MAX_NGBS_PER_PNT - 1; j; j--)
	if (nngbs_histo[j] != 0 && nngbs_histo[j-1]==0) 
	  h=j;

      if( !oldi && ((h-l) > maxspread) ) { /* spread too big  */
	current_max=h;			/* adjust it; h-1 too restrictive */
	oldi=i;				/* save i for restoring later  */
      }
      if (oldi && ((h-l) <= maxspread)) { /* relax again */
	current_max=MAX_NGBS_PER_PNT;
	i=oldi;				/* restore old i, and see if */
	oldi=0;				/* anything is left */
      }
    }
    
    edgei=candp_array[i];

    if (edgei->face_with == UNAVAIL)
      continue;

    if (Points[edgei->ends[0]].nngbs >= current_max) { 
      CNT(TOO_MANY_A); continue;}
    if (Points[edgei->ends[1]].nngbs >= current_max) { 
      CNT(TOO_MANY_B); continue;}


    VEC3PLUS(edgei->midpoint,
	     Points[edgei->ends[0]].xyz, Points[edgei->ends[1]].xyz);
    NRMIZE3(edgei->midpoint);

    ncross=0;

    for (j=0; j<nedges; j++) {
      edgej=Edges+j;
      if ( (ncross += edge_crossing(Sphere,edgei, edgej)) ) 
	break;
    }

    if (ncross == 0) {			/* no crossing: found real new edge */
      Edges[nedges]= *edgei;		/* copy it to final list */
      edgei->face_with=UNAVAIL;		/* and take out of candidate list */
      edgei=Edges+nedges;

      sp = &Points[edgei->ends[0]].nngbs; /* pointer to nr of neighbours */
      Points[edgei->ends[0]].edges[ *sp] = nedges; 
      ++*sp;				/* update the actual number */
      nngbs_histo[*sp]++;		/* update the histogram */
      nngbs_histo[(*sp)-1]--;		/* <- don't forget !  */

      /* same for other end : */
      sp = &Points[edgei->ends[1]].nngbs;
      Points[edgei->ends[1]].edges[ *sp] = nedges; 
      ++*sp; 
      nngbs_histo[*sp]++; 
      nngbs_histo[(*sp)-1]--;

      nedges++;
      if (nedges == 3*npoints -6)
	break;
    }
  }

  //  AFREEA(cand_edges);     @@@@ bug : segfault here
  //  AFREEA(candp_array);    @@@@ bug : segfault here too

  Sphere->nedges=nedges; 
  Sphere->Edges=Edges; 
  return;
} /* find_edges */

/* static int edge_cmp(edge_p *a, edge_p *b) {
 *   if ((*a)->length2 < (*b)->length2) return -1;
 *   if ((*a)->length2 > (*b)->length2) return 1;
 *   return 0;
 * }
 */

static int edge_cmp(const void *a, const void *b) {
  if ((*(edge_p*)a)->length2 < (*(edge_p*)b)->length2) return -1;
  if ((*(edge_p*)a)->length2 > (*(edge_p*)b)->length2) return 1;
  return 0;
}

static int edge_crossing(sphere_p Sphere, edge_p a, edge_p b) {
  /*
   * does    : checks if edges a and b are crossing
   *
   * gets    : edges
   *
   * affects : nothing
   *
   * returns : 0 if all's ok, 1 if crossing
   *
   * warns   : 
   *
   * comment : returns as soon as possible
   */


  int i=(int)a->ends[0], j=(int)a->ends[1]
    ,k=(int)b->ends[0], l=(int)b->ends[1];
  int nc=0;

  point_p Points=Sphere->Points;

  float * aa=Points[i].xyz;
  float * ab=Points[j].xyz;
  float * ba=Points[k].xyz;
  float * bb=Points[l].xyz;
  double am[3], a1[3], a2[3], bm[3], b1[3], b2[3];
  double d2;

#ifdef UNDEFINED
#define BREAK nc++
#define CROSSES 
#else
#define BREAK return 0
#define CROSSES return 1
#endif

  /* first do the cheap tests: */
  if ( (i == k) || (j == l) || (i == l) || (j == k) ) {
    CNT(COM_NGB); BREAK;
  } else { CNT(NO_COM_NGB); }


  /* make sure we always index properly, since matrix is upperdiagonal: */
#define DISTANCES2(I,J) (Distances2[ MIN(I,J)  ][ MAX(I,J) ])

  if (DISTANCES2(i, k) > 1.5 &&
      DISTANCES2(i, l) > 1.5 &&
      DISTANCES2(j, k) > 1.5 &&
      DISTANCES2(j, l) > 1.5) { CNT(END_POINT_OTHER_SIDE); BREAK; }
					/* too far apart anyway */
  if (DISTANCES2(i, k) > cutoff2 &&
      DISTANCES2(i, l) > cutoff2 &&
      DISTANCES2(j, k) > cutoff2 &&
      DISTANCES2(j, l) > cutoff2) { CNT(TOO_FAR); BREAK; }


  if ( (d2=DIST3_2( a->midpoint, b->midpoint)) >= 1.0) {
      CNT(MIDPOINT_OTHER_SIDE); BREAK; }

  if (d2 < MIDPOINTS_TOO_CLOSE_TOL*cutoff2) {
    CNT(MIDPOINTS_TOO_CLOSE); CROSSES; }

  if ( (d2=ANGLE( a->midpoint, b->midpoint)) >= arccutoff) {
      CNT(MIDPOINT_ANGLE_TOO_LARGE); BREAK; }

					/* too far apart anyway */


#define CROSSDOT(A,B,M) (DOT3P(A,B)*DOT3P(M,M) - DOT3P(A,M)*DOT3P(B,M))
  /* comes from ( A x M )  . ( B x M ), which should be > 0 */
  VEC3MIN(am,ab,aa);
  VEC3MIN(bm,bb,ba);
  VEC3MIN(a1,ba,aa);
  VEC3MIN(a2,bb,aa);
  VEC3MIN(b1,ab,ba);
  VEC3MIN(b2,aa,ba);

  if ( (d2=CROSSDOT(a1,a2,am)) > 0.0 ) {
    CNT(NOCROSS1); BREAK; }
  else { CNT(CROSS1); }
  if ( (d2=CROSSDOT(b1,b2,bm)) > 0.0 ) {
    CNT(NOCROSS2); BREAK; }
  else { CNT(CROSS2); }

  if(nc)
    return 0;
  else
    return 1;
} /* edge_crossing */

double estimate_arclen(int npoints) {

/* surface of elliptical triangle (one on surface of sphere), in arclengths :
 * tan^2(epsilon/4) = tan(s/2)*tan(s-a/2)*tan(s-b/2)*tan(s-c/2), where s =
 * (a+b+c)/2 (cf. case of flat triangles). Pretend they are all equilateral
 * spherical triangles, and solve: (tan(3x) = 3tan(x)-tan^3(x)/(1-3tan^2(x))
 */
  int i, nfaces_pre;			/* nr of solutions */
  double R=1.0;
  double sol[3], x, c;

  /* == tan(075x)*tan^3(0.25x): */
  /* solve. That is, solve -t^6 + 3t^4 + 3ct^2 - c==0, t = tan(0.25*x), */
  /* c=tan^2(eps/4),  or p^3 - 3p^2 - 3cp^2 + c == 0 with p==t^2 */

  nfaces_pre= 2*npoints-4;
  c=SQR(tan(M_PI*SQR(R)/nfaces_pre));
  i=root3(1.0, -3.0, -3.0*c, c, sol);	/* solve cubic equation */
  assert( i==3); /* must be three solutions (empirically)  */
  x = sol[1];				/* (empirically so) */
  assert( x > 0.0 && x < 1.0);		/* (empirically so) */

  x=4*atan( SQRT(x));			/* approximate arclength */
  return  x ;
/*   x = 2.0* R * sin(0.5*x);		/# approximate edge length #/
 *   return  x ;
 */
} /* estimate_arclen */


/* enough neighbours can be found within CONV/npoints (unit sphere) */
/*   very roughly: nv=npoints; ne=0.5*(6 edges per point) since counted
 *   double; nv-ne+nf = 2 (Euler), so nf =~ 2n. 2n*A = 4piR^2, so area per
 *   triangle A=2pi/n.  unit radius, roughly equilateral triangles:
 *   A^2=1/16*(6x^4 - 3x^4); A=1/4V3*x^2=2pi/n => x^2 ~=~ 8V3pi/3n => x ~
 *   SQRT(14.51/n). Take that times 2 or so.
 *
 */

#define OTHEREND(EP, W) ( (EP)->ends[ ( (EP)->ends[0] == W ) ] )

static void find_faces(sphere_p Sphere) {
  point_p Pi, Pj;
  int i,j,k, ngbs, p1, p2, p3
    ,e1,e2,e3,e4;
  SHORT *sp;
  edge_p Edgep, Nedgep;
  int nfaces_pre;

  point_p Points=Sphere->Points;
  int npoints=Sphere->npoints;
  edge_p Edges=Sphere->Edges; 
  face_p Faces=Sphere->Faces; 
  int nfaces=Sphere->nfaces;


/*
 * p3 - e3 - p2
 *  \  	    /
 *   e2	   e1
 *     \  /
 * ---	p1 ---
 */

  nfaces_pre=2*npoints-4;
  Faces=(face_p)CALLOC( nfaces_pre, sizeof(face_t));

  for (i=0; i<Sphere->nedges; i++) 
    Sphere->Edges[i].face_with=AVAIL;	/* a number that never can be that */
					/* of another face */

  for (i=0; i<npoints; i++) {
    p1=i; 
    Pi=Points+i;
    ngbs=Pi->nngbs;
    for (j=0; j < ngbs; j++) {
      Edgep=Edges+( (e1=Pi->edges[j]) );
      if ( Edgep->face_with == UNAVAIL ) /* edge already part of 2 faces */
	continue;
      p2 = OTHEREND(Edgep, p1); 

      e2=Pi->edges[(j + 1)%ngbs] ;	/* cyclic! */

      Nedgep=Edges+e2;
      p3 = OTHEREND(Nedgep, p1);

      if (p3==Edgep->face_with)		/* found same face a 2nd time ! */
	continue; 
      /* find the 'missing' edge from the rays around p2 and p3  */
      Pj=Points+p2;
      for (k=0; k<=Pj->nngbs; k++)
	if (Pj->edges[k]==e1)		/* found original edge */
	  break;
      if (k == (Pj->nngbs+1))
	DIE("bug: didn't find edge e1 back in edge list of neibor point,");

      e3= Pj->edges[ (k+Pj->nngbs-1)%Pj->nngbs ]; /* the clock wise edge */

      /* check that with p3 and e2: */
      Pj=Points+p3;
      for (k=0; k<=Pj->nngbs; k++)
	if (Pj->edges[k]==e2)		/* found original edge */
	  break;
      if(k == (Pj->nngbs+1))
	DIE("bug: didn't find edge e2 back in edge list of neibor point,");
      e4 = Pj->edges[ (k+1)%Pj->nngbs];
      if ( e4 != e3)
	DIE("\
topology error: edge %d (%d - %d), %d (%d - %d), and %d (%d - %d)\n\
exspected to form one triangle; instead extra edge involved: %d (%d - %d),",
			   (int)e1, (int)Edges[e1].ends[0], 
			   (int)Edges[e1].ends[1],
			   (int)e2, (int)Edges[e2].ends[0], 
			   (int)Edges[e2].ends[1],
			   (int)e3, (int)Edges[e3].ends[0], 
			   (int)Edges[e3].ends[1],
			   (int)e4, (int)Edges[e4].ends[0], 
			   (int)Edges[e4].ends[1]);
      /* OK, so now we've found a new face; add it, and update: */
      Faces[nfaces].edges[0]=e1;
      Faces[nfaces].edges[1]=e2;
      Faces[nfaces].edges[2]=e3;
      Faces[nfaces].points[0]=p1;
      Faces[nfaces].points[1]=p2;
      Faces[nfaces].points[2]=p3;
      nfaces++; 
      sp= &Edges[e1].face_with;      *sp = (*sp==AVAIL) ? p3 : UNAVAIL; 
      sp= &Edges[e2].face_with;      *sp = (*sp==AVAIL) ? p2 : UNAVAIL; 
      sp= &Edges[e3].face_with;      *sp = (*sp==AVAIL) ? p1 : UNAVAIL; 
      /* keep the point with which the face is formed, such that each edge */
      /* can only be used twice */
    }
  }
  Sphere->nfaces=nfaces;
  Sphere->Faces=Faces;
  return;
} /* find_faces */

static double sort_rays_angles[ MAX_NGBS_PER_PNT ];

static void sort_rays(sphere_p Sphere, int p) {
  point_p Point;
  edge_p Edgep;
  int i, ngbs, q, e, array[MAX_NGBS_PER_PNT]
    , tmp[MAX_NGBS_PER_PNT], maxq;
  double e1[3], e2[3];			/* local xy vectors */
  float *Vec, *fp, d; 
  double x,y, f, z, maxz= -100.0;

  point_p Points=Sphere->Points;
  edge_p Edges=Sphere->Edges; 
  
  Point = Points+p;
  ngbs=Point->nngbs;
  Vec=Point->xyz;
  z = Vec[2];

  for (i=0; i<ngbs; i++) {		/* find the edge which has */
    e=Point->edges[i];			/* smallest angle with +z axis */
    Edgep=Edges+e; 
    q= OTHEREND(Edgep, p);
    f = Points[q].xyz[2] - z; 
    if ( f > maxz ) {
      maxz=f;
      maxq=q;
    }
  }
  fp=Points[maxq].xyz;			/* now form local xy system */

  d=DOT3P(Vec,fp);			/* project on plane normal to p */
  e1[0]=fp[0]-d*Vec[0];
  e1[1]=fp[1]-d*Vec[1];
  e1[2]=fp[2]-d*Vec[2];
  NRMIZE3(e1);				/* that's x */
  CROSS3P(e2, Vec, e1);			/* that's y */

  for (i=0; i<ngbs; i++) {		/* angles of vectors to neibors */
    e=Point->edges[i];			/* form local xy system */
    Edgep=Edges+e;

    q= OTHEREND(Edgep, p);

    fp=Points[q].xyz;
    x=DOT3P(fp, e1);
    y=DOT3P(fp, e2);
    f=atan2(y,x);
    if (f < -0.01)			/* take care of roundoff */
      f += M_2PI;			/* because -pi <atan2< pi, */
    sort_rays_angles[i]= f;
    array[i]=i;				/* indirect acces via static array */
  }
  qsort(array, ngbs, sizeof(int), ray_cmp);
  /*   assert( array[0] == 0); no!! */

  for (i=0; i<ngbs; i++)
    tmp[i]=Point->edges[i];

  for (i=0; i<ngbs; i++)		/* rearrange them */
    Point->edges[i]=tmp[array[i]];

  return;
} /* sort_rays */

static int ray_cmp(const void * a, const void *b) {
  if ( sort_rays_angles[*(int*)a] < sort_rays_angles[*(int*)b])
    return -1;
  if ( sort_rays_angles[*(int*)a] > sort_rays_angles[*(int*)b])
    return 1;
  return 0;
} /* ray_cmp */

static void areas(sphere_p Sphere) {
  int i,j, p0, p1, p2;
  float *x0, *x1, *x2;
  double s,e, a0,a1,a2, m[3];
  face_p facep;
  edge_p e0, e1, e2;

  point_p Points=Sphere->Points;
  edge_p Edges=Sphere->Edges; 
  face_p Faces=Sphere->Faces; 
  int nfaces=Sphere->nfaces;
  
  for (i=0; i<nfaces; i++) {
    facep=Faces+i;
    e0=Edges+(facep->edges[0]);
    e1=Edges+(facep->edges[1]);
    e2=Edges+(facep->edges[2]);

    p0=e0->ends[0];
    p1=e0->ends[1];
    j=e1->ends[0];
    p2= (j==p0 || j == p1) ? e1->ends[1] : j ;
    x0=Points[p0].xyz;
    x1=Points[p1].xyz;
    x2=Points[p2].xyz;
    VEC3ASS(m, x0);			/* first the midpoint */
    VEC3INC(m, x1);
    VEC3INC(m, x2);
    NRMIZE3(m);
    VEC3ASS(facep->midpoint, m);

#define ARC(I) (2.0*ASIN(0.5*SQRT(Edges[facep->edges[I]].length2)))

    a0=ARC(0);				/* the arc lengths */
    a1=ARC(1);
    a2=ARC(2);
    facep->area=(float)lHuilier(a0,a1,a2);
		
    x2=facep->midpoint;
    s=0.0; 
    for (j=0; j<3; j++) {		/* now the point-area thing*/
      p0=Edges[facep->edges[j]].ends[0];
      x0=Points[p0].xyz;
      x1=Edges[facep->edges[j]].midpoint;

      a0=DOT3P(x0, x1);
      a0=ACOS(a0);
      a1=DOT3P(x1,x2); 
      a1=ACOS(a1);
      a2=DOT3P(x2, x0); 
      a2=ACOS(a2);

      e=lHuilier(a0,a1,a2);
      s+=e;				/* to test later on */
      Points[p0].area += (float)e;

      p0=Edges[facep->edges[j]].ends[1]; /* and the same for the other half */
      x0=Points[p0].xyz;
      a0=DOT3P(x0, x1);
      a0=ACOS( a0 );
      a2=DOT3P(x2, x0);
      a2=ACOS( a2 );

      e=lHuilier(a0,a1,a2);
      s+=e;				/* to test later on */
      Points[p0].area += (float)e;
    } /* for j < 3 */
    assert(FABS(facep->area  - s ) < AREA_TOL);
  } /* for i < nfaces */
  return;
} /* areas */

double lHuilier(double a, double b, double c) { 
  /* returns area of triangle with given arclengths */
  double s, t, e;			

  s=0.5*(a+b+c);
  
  t=tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c));
  if (t < 0.0)
    DIE("impossible triangle: triangle inequality doesn't hold,");
  e=4.0*atan(SQRT(t));
  return e;
} /* lHuilier */

static int sort_points_dim=DFLT_DIM;
static double sort_points_spacing=DFLT_SPACING;

#define OFFSET(DIM, X,Y,Z) DIM*(DIM*X + Y)+Z

static int point_cmp(const void * a, const void * b) {
  int ax, ay, az, bx, by, bz, offseta, offsetb;

#define INTDIM(X, I) ((int)floor((((float*)(X))[(I)] + 1.0)/\
				 sort_points_spacing))

  ax=INTDIM(a, 0);
  ay=INTDIM(a, 1);
  az=INTDIM(a, 2);

  bx=INTDIM(b, 0);
  by=INTDIM(b, 1);
  bz=INTDIM(b, 2);

  offseta=OFFSET(sort_points_dim, ax,ay,az);
  offsetb=OFFSET(sort_points_dim, bx,by,bz);

  return offseta - offsetb;
} /* int point_cmp */
#undef OFFSET

static void sort_points(sphere_p Sphere) {
  /*
   * does    : sort points according to boxing, allocates Points and sets
   *   their ->nr and ->xyz.
   *
   * gets    : Sphere
   *
   * affects : Sphere->Points, order of Sphere->coordinates
   *
   * returns : nothing
   *
   * warns   : 
   *
   * comment : sorts into a grid of dimension DFLT_DIM, but can be altered
   */
  int i, npoints;
  point_p Points, point;
  float * coordinates, *x;

  coordinates=Sphere->coordinates;
  npoints=Sphere->npoints;

/*   if(opt_general) {			/# may leave this out #/
 *     sort_points_dim=opt_general;
 *     sort_points_spacing=2.0/opt_general;
 *   }
 */
  qsort(coordinates, npoints, sizeof(float[3]), point_cmp);

  Points=(point_p)CALLOC(npoints, sizeof(point_t));

  for (i=0,   point=Points,   x=coordinates;
       i<npoints; 
       i++, point++, x+=3) {
    point->xyz=x;
    point->nr=i;
  }
  Sphere->Points=Points;
} /* sort_points */

void generate_template_spheres(void) {
  /* generate different subdivisions of the icosahedral and the tesselated */
  /* dodecahedron (pdod) to be used as templates. Sure, it can be done */
  /* more intelligent, but isn't worth the effort */
  int i, n;
  sphere_p Sphere;

  warn("doing subdivision of order points edges faces:\n");
  for (i=0; i<N_TMPLT_SPHERES; i++)  {
    n=Template_npoints[i];

    warn("%2d %3d %3d %3d\n", i, n, 3*n-6, 2*n-4);

    Sphere=tesselate(Template_based_on[i], Template_order[i]);

    randomize_points(Sphere);		/* 'twist' points slightly, to  */
					/* put them within <-1, 1> interval */
    sort_points(Sphere);

    triangulate_sphere(opt_careful, Sphere); /* just do it */

/*    boxify_points(Sphere);		/# put templates point in grid #/
 */

    Template_spheres[i]=Sphere;
  }
  warn("done\n");
  return;
} /* generate_template_spheres */

void free_templates(void) {
  int i;
  sphere_p Sph;

  for (i=0; i<N_TMPLT_SPHERES; i++) {
    Sph = Template_spheres[i];
    FREE(Sph->coordinates);
    FREE(Sph->Points);
    FREE(Sph->Edges);
    FREE(Sph->Faces);
    FREE(Sph);
  }
} /* free_templates */



static float box_length;
static htable_p point_hash;

static void node_free(void *node) {
  if(node) 
    FREE(node);
} /* node_free */

sphere_p recursion_tessel(int type, int freq) {
  /*
   * does    : generates a spherical tesselation by recursive N-division of 
   *   an the faces of an icosahedron (type==SPH_ICOS) or tesselated
   *   dodecahedron (type== SPH_PDOD)
   *
   * gets    : frequency
   *
   * affects : 
   *
   * returns : a sphere
   *
   * warns   : about purging
   *
   * comment : 
   */
  int i, n, npoints, nfaces, p;
  hnode_p node;
  float edgelen;
  sphere_p sphere;
  float  * points, *coordinates;
  int * faces;

  switch(type) {
  case SPH_ICOS:   
    npoints = 10 * SQR(freq) + 2; 
    nfaces = 20;
    points = (float*)icos_points;
    faces  = (int *)icos_faces;
    break;
  case SPH_PDOD:
    npoints = 30 * SQR(freq) + 2; 
    nfaces = 60;
    points = (float *)pdod_points;
    faces  = (int *)pdod_faces;
    break;
  default: die("unknown tesselation type: %d\n", type);
    break;
  }

  edgelen=estimate_arclen(npoints);
  edgelen=2*sin(0.5*edgelen);			/* estimated edge length */
  box_length= M_1SQRT3*edgelen;
  point_hash = Hcreate(NULL, 3*npoints, sizeof(int[3]) ) ;
  
#define FACE_POINT(J) points + 3*(faces[3*i+J])

  for (i=0; i<nfaces; i++) 
    tesselate_face(freq, FACE_POINT(0), FACE_POINT(1), FACE_POINT(2));
  

  n=point_hash->nkeys;
  if (npoints != n) { 
    warn("exspected %d, found %d points; purging ...\n", npoints, n);
    assert(n > npoints);		/* equiv. points in neiboring boxes */
    p=purge(point_hash, edgelen);
    warn("purged %d points ", p);
    if (p == n - npoints)
      warn("succesfully\n");
    else
      die("in vain\n");
  }

  sphere=(sphere_p)CALLOC(1,sizeof(sphere_t));
  coordinates=sphere->coordinates=(float*)CALLOC(npoints, sizeof(float[3]));

  n=0;
  for(node=Hkeys(point_hash); node; node=node->next) {
    if(node->value==NULL)		/* was purged */
      continue;
    VEC3ASS(coordinates +3*n, (float*)node->value);
    n++;
  }
  assert(n==npoints);
  sphere->npoints=n;

  Hfree(point_hash, node_free);
  return sphere;
} /* recursion_tessel */


static void tesselate_face(int freq, float * a, float *b, float *c) {
  /*
   * does    : does a recursive tesselation of triangle abc into freq^2
   *   triangles. 
   *
   * gets    : edge points of triangle
   *
   * affects : 
   *
   * returns : 
   *
   * warns   : 
   *
   * comment : this is a recursive function. Results are collected in a
   *   hash table
   */

  int f,i,j;
  float * triangles;

  if (freq==1) {			/* end of the recursion: store point */
    store_point(a);
    store_point(b);
    store_point(c);
    return;
  }

  for (f=2; freq % f !=0; f++)		/* factorize freq */
    ;

#define INDEX(I,J)  ((I*(I+1))/2 + J ) /* index of point's coordinates */
#define POINT(I,J) (triangles + (3*INDEX((I),(J))))
#define TRIANGLE1(I, J) POINT(I,J), POINT(I, J+1), POINT(I-1,J)
#define TRIANGLE2(I, J) POINT(I,J), POINT(I, J+1), POINT(I+1,J+1)

  /* f is the (current) tesselation frequency. Tesselate with it  */
  triangles=subdiv_triangle(f, a, b, c);

  /* and subsequently tesselate these faces, recursively */
  for (i=1; i<=f; i++) 
    for (j=0; j<i; j++) {
      tesselate_face(freq/f, TRIANGLE1(i,j));
      if (i != f )
	tesselate_face(freq/f, TRIANGLE2(i,j));
    }
  FREE(triangles);
  return;
} /* tesselate_face */

float * subdiv_triangle(int freq, float * a, float *b, float *c) { 
  /*
   * does    : directly subdivide one triangle into freq^2 smaller triangles
   *
   * gets    : edge points of triangle
   *
   * affects : 
   *
   * returns : a dynamically allocated array of float[3]'s, holding
   *   coordinates of the subdivided triangle
   *
   * warns   : 
   *
   * comment : for triangle abc, uses i,j as 'coordinates' during
   *   subdivision. For the triangles bca and cab, first maps i, and j to
   *   their own coordinates u and v, so as to keep same order.
   */

  int i, j, u, v, n;
  int npoints;
  float * triangles, *where, *p;
  double mat[3][3], x[3], y[3], z[3], angle, uf, vf, rf;

  npoints = (freq+1)*(freq+2)/2;	/* number of points to be generated */
  triangles = (float*)CALLOC(npoints, sizeof(float[3]));

  /* corner points: */
  VEC3ASS(triangles,   a);
  VEC3ASS(triangles+3*INDEX(freq,0), b);
  VEC3ASS(triangles+3*INDEX(freq,freq), c);

#define DO_EDGE(A,B,MAIN,SECOND,DERIV)					   \
  VEC3ASS(y, A);			/* starting vector */		   \
  angle = ANGLE(A,B);							   \
  angle /= freq;							   \
  CROSS3P(z, A,B);							   \
  NRMIZE3(z);				/* axis about which to rotate */   \
									   \
  ARBANGROTMAT(mat, z, angle);		/* matrix rotating A towards B */  \
									   \
  for (MAIN=1; MAIN<=freq-1; MAIN++) {	/* main looping variable */	   \
    SECOND DERIV;			/* set other coordinate */	   \
    n=INDEX(i,j);			/* where to leave things */	   \
    where = triangles + 3*n;						   \
    VEC3MAT(x, y, mat);							   \
    VEC3INC(where, x);							   \
    VEC3ASS(y,x);							   \
  }

  /* edge points: */

  DO_EDGE(a,b,i,j, =0);			/* edge ab */
  DO_EDGE(b,c,j,i, =freq);		/* edge bc */
  DO_EDGE(a,c,i,j, =i);			/* edge ac */

#define FILL_IN(I,J, A,B,C) u=I; v=J;					\
      			    if (u!=0 && u != freq) {			\
                              uf= u*rf; vf= v*rf;			\
			      p = interior_point(A, B, C, uf, vf);	\
			      VEC3INC(where, p); 			\
			    }
  /* interior points: */
  rf = 1.0/(double)freq;
  for (i=1; i<freq; i++)
    for (j=1; j<i; j++) {
      n=INDEX(i,j);
      where = triangles + 3*n;

      /* ab = x-axis, ac = y-axis: abc triangle */
      FILL_IN(i, j, a, b, c);

      /* bc = x-axis, ba = y-axis: bca triangle */
      FILL_IN(freq-i+j, freq-i, b, c, a);

      /* ca = x-axis, cb = y-axis: cab triangle */
      FILL_IN( freq-j, i-j, c, a, b);
    }

  /* normalize */
  for (i=0; i<npoints; i++) {
    where = triangles + 3*i;		/* or use INDEX ?  */
    NRMIZE3(where);			/* <== something wrong here  */
  }
  
  return triangles;
} /* subdiv_triangle */


float * interior_point(float * a, float * b, float *c, 
			   float u, float v) {
  /*
   * does    : calculates point lying on v/length part of the arc v1-v2 that's
   *   formed by connecting the u/length of arc a-b and arc a-c.  
   *
			  _ xy axis       
			  /|
			 / 
	                c 
	               / \  
         (angle2)  v2 /   \  
	             `     \ 
                    / `     \ 
       y-axis _    /   `     \ 
	     |\   /     `res. \ (angle3)
               \ /       `     \ 
                a---------`---- b -> x axis
         	          v1 (angle1)

   * gets    : 'origin' a, 'x-axis' pointing to b, 'yx-axis' pointing to
   *   c, and x-coordinate i, y-coordinate j
   *
   * affects : nothing
   *
   * returns : pointer to static float, holding coordinates
   *
   * warns   : 
   *
   * comment : 
   */
  double p[3], v1[3], v2[3], angle1, angle2, angle3;
  double  mat[3][3];
  static float result[3];

  VEC3INIT(result, 0.0);		/* reset */

  angle1=ANGLE(a,b);
  angle1 *= u;				/* angle between origin and x coord */

  CROSS3P(p, a,b);
  NRMIZE3(p); 
  ARBANGROTMAT(mat, p, angle1);
  VEC3MAT(v1, a, mat);
  /* v1 a point with right x coordinate lying on x axis */

  angle2=ANGLE(a,c);
  angle2 *=u;				/* angle between origin and */
					/* xy-intercept */

  assert( FABS(angle2 - angle1) < ANGLE_TOL); /* should not be too far apart */

  CROSS3P(p, a,c);
  NRMIZE3(p); 
  ARBANGROTMAT(mat, p, angle2);
  VEC3MAT(v2, a, mat);
  /* v1 a point with right y coordinate, lying on xy axis */

  angle3=ANGLE(v1, v2);
/*   angle /= u ; angle *= v; */
  angle3 *= (v/u);			/* arc length of v1-v2 */

  CROSS3P(p, v1,v2);
  NRMIZE3(p); 
  ARBANGROTMAT(mat, p, angle3);

  VEC3MAT(result,v1,mat);
  NRMIZE3(result);
  return result;
} /* interior_point */

static int redundant;
static int non_redundant;

static void store_point(float*xyz) {
  int i, ia[3];
  float *y, d;
  hnode_p node;

  for (i=0; i<3; i++) 
    ia[i] = (int)floor(xyz[i]/box_length);

  node = Hget(point_hash, ia);
  if (node) {
    y= (float*)node->value;
    d=DIST3_2(y, xyz);
    if (d > BOXING_TOL)
      die("clash: %8.3g %8.3g %8.3g <-> %8.3g %8.3g %8.3g\n",
	  y[0], y[1], y[2], xyz[0], xyz[1], xyz[2]);
    redundant++;
  } else { 
    node = Hset(point_hash, ia);
    node->value = memsav(xyz, sizeof(float[3]));
    non_redundant++;
  }
  return;
} /* store_point */

static int purge(htable_p point_hash, float edge_length) {
  /*
   * does    : tries to delete points that didn't hash nicely to where
   *   belong, due to round off error
   *
   * gets    : the hash_table to be recovered
   *
   * affects : the hash_table: some keys are deleted
   *
   * returns : new number of entries
   *
   * warns   : about keys deleted
   *
   * comment : 
   */

  hnode_p node, node2;
  float *oldx, d, maxd, *newx;
  int i, j, k, *oldkey, newkey[3], purged;

  maxd = 0.1*SQR(edge_length);
  purged=0;

  for(node = Hkeys(point_hash); node; node=node->next) {
    oldkey=node->key;
    oldx = (float*)node->value;
    if (!oldx)
      continue;				/* nulled before */
    /* check the 26 neighbours for occurrence of a guy that doesn't belong */
    /* there: */
    for (i=-1; i<2; i++) 
      for (j=-1; j<2; j++) 
	for (k=-1; k<2; k++) {
	  if (i==0 && j==0 && k==0)	/* central one */
	    continue;
	  newkey[0]=oldkey[0]+i;
	  newkey[1]=oldkey[1]+j;
	  newkey[2]=oldkey[2]+k;
	  node2=Hfeel(point_hash, newkey);
	  if (!node2 || node2->value==NULL)	/* empty, or nulled before */
	    continue;
	  newx=(float*)node2->value;
	  d = DIST3_2(newx, oldx);
	  if (d < maxd)	{		/* have one! */
	    purged++;
	    FREE(node2->value);
	    node2->value=NULL;
	  }
	}
  }
  return purged;
} /* purge */

static int make_Delaunay(sphere_p Sphere) {
  /*
   * does    : checks if the triangulation is a Delaunay-one, and repairs
   *   if needed.
   *
   * gets    : a complete triangulation in Sphere
   *
   * affects : some sphere->Edges[*] and sphere->Points[*]->whoamis
   *
   * returns : nr of edges that were reassigned
   *
   * warns   : about nr of passes that where made
   *
   * comment : works by checking the quadrilateral of which an edge is the
   *   diagonal, for all edges
   */
  int i,nedges=Sphere->nedges, a, b, c, d, wa, wb, nwrong=0, collinear, t
    ,npass=0 , 	nreplacements=0, lowest=999999, l;
  edge_p Edges=Sphere->Edges, edge, ec,ed;
  point_p Points=Sphere->Points, pa, pb, pc, pd;


  do  {
    npass++;
    nwrong=0;				/* reset that */
    for (i=lowest==999999?0:lowest; i<nedges; i++) {
      collinear=0;
      edge=Edges+i;
      a=edge->ends[0];
      b=edge->ends[1];
      pa=Points+a;
      pb=Points+b;
      wa=edge->whoami[0];
      wb=edge->whoami[1];

      /* find point c, from point a: */
      ec=Edges + pa->edges[ (wa + 1)%pa->nngbs ];
      c = ec->ends[ a == ec->ends[0] ];

      /* check this by starting from b: */
      ec=Edges + pb->edges[ (wb + pb->nngbs - 1)%pb->nngbs ];
      assert(c == ec->ends[ b == ec->ends[0] ]);

      /* find point d, starting from b: */
      ed=Edges + pb->edges[ (wb + 1)%pb->nngbs ];
      d= ed->ends[ b == ed->ends[0] ];

      /* check this by starting from a: */
      ed=Edges + pa->edges[ (wa + pa->nngbs - 1)%pa->nngbs ];
      assert(d == ed->ends[ a == ed->ends[0] ]);

      if (c > d) { t=c; c=d; d=t; }
      
      pc = Points+c;
      pd = Points+d;

      if (spherical_delaunay(pa->xyz, pb->xyz, pc->xyz, pd->xyz, &collinear)) {
	nwrong++;
	nreplacements++;
	warn("replacing edge %d: %d - %d by %d - %d ...\n",
	     i, a,b,c,d);
	l=replace_edge(Sphere, i, a, b, c,d);
	if ( l<lowest )
	  lowest=l;
      }
#ifdef UNDEFINED			/* occurs too often ... */
      if (collinear)
	warn("%d collinear triangles around edge %d: %d - %d \n",
	     collinear, i, a,b,c,d);
#endif

    } /* for(i looping over edges) */
  } while(nwrong);

  if (npass != 1)
    warn("needed %d passes to render triangulation Delaunay\n", npass-1);
      
  return nreplacements ;
} /* make_Delaunay */


static int spherical_delaunay(float * a, float *b, float *c, float *d, 
			      int * colp) {
  /*
   * does    : checks if quadrilateral abcd fulfills delaunay property
   *
   * gets    : points lying on unit sphere
   *
   * affects : *colp is incremented if a collinear triangle is found
   *
   * returns : 0 if first pair of points as edge is delaunay, 1 if 2nd pair
   */

/* For edge a b, check if the other diagonal of the quadrilateral it is
 * the diagonal of, would yield a better triangulation:
 * 
 *     c
 *    /|\
 *   / | \
 *  /  |  \ 
 * a---+---b
 *  \  |  /
 *   \ | /
 *    \|/
 *     d
 * 
 */
  int i;
  /* edges: */
  double ab=ACOS(DOT3P(a,b));
  double bc=ACOS(DOT3P(b,c));
  double ac=ACOS(DOT3P(a,c));
	    
  double ad=ACOS(DOT3P(a,d));
  double bd=ACOS(DOT3P(b,d));
	    
  double cd=ACOS(DOT3P(c,d));

  double angles0[6], angles1[6];
  double min0=1000.0, min1=1000.0, max0=-1000.0, max1=-1000.0;

  /* assuming a-b is the proper edge: */
  /* triangles abc: */
  angles0[0]=spherical_angle(ab, bc, ac, colp);
  angles0[1]=spherical_angle(bc, ac, ab, colp);
  angles0[2]=spherical_angle(ac, ab, bc, colp);
  /* triangles abd: */
  angles0[3]=spherical_angle(ab, bd, ad, colp);
  angles0[4]=spherical_angle(bd, ad, ab, colp);
  angles0[5]=spherical_angle(ad, ab, bd, colp);


  /* now assuming cd is the proper edge: */
  /* triangles acd: */
  angles1[0]=spherical_angle(ac, cd, ad, colp);
  angles1[1]=spherical_angle(cd, ad, ac, colp);
  angles1[2]=spherical_angle(ad, ac, cd, colp);
  /* triangles bcd: */
  angles1[3]=spherical_angle(bc, cd, bd, colp);
  angles1[4]=spherical_angle(cd, bd, bc, colp);
  angles1[5]=spherical_angle(bd, bc, cd, colp);
  
  for (i=0; i<6; i++) {
    if (angles0[i] < min0)
      min0=angles0[i];
    if (angles0[i] > max0)
      max0=angles0[i];
    if (angles1[i] < min1)
      min1=angles1[i];
    if (angles1[i] > max1)
      max1=angles0[i];
  }
  return min0  < min1 ; 
} /* spherical_delaunay */

static double spherical_angle(double a, double b, double c, int * colp) {
  /*
   * does    : calculates angle (opposite edge b) on unit sphere, of
   *   spherical triangle with edge (arc) lenghts a, b & c
   *   *colp is incremented for every collinear triangle
   */

  double s = (a+b+c)/2.0;
  double r;

  r = (sin(s)*sin(s-a));

  if (FABS(r) < COLLINEAR_TOL) { 
    (*colp)++;
    return M_PI;
  }
  else
    r= sin(s-b)*sin(s-c)/r;
  if (r < 0.0 ) {
    assert ( -r < COLLINEAR_TOL);
    return 0.0;
  }
  return 2.0*atan(sqrt(r)) ;
  
} /* spherical_angle */

static int replace_edge(sphere_p Sphere, int e, int a, int b, int c, int d) {
  /*
   * does    : replaces edge i by one that has c and d as endpoints,
   *   instead of a and b
   *
   * gets    : Sphere, edge number, and corner points of the
   *   quadrilateral. a < b, and c < d ! 
   *
   * affects : Sphere->Edges[i], and some of the Points->edges[*]
   *
   * returns : the lowest edge number involved, (i.e. where to start
   *   scanning again ...)
   *
   * warns   : 
   *
   * comment : 
   */

  int i,j, al, ar, bl, br;
  edge_p Edges=Sphere->Edges, edge=Edges+e, ep;
  point_p Points = Sphere->Points, pa, pb, pc, pd, pp;
  int new_edges[MAX_NGBS_PER_PNT];
  EDGE_ID *edge_array;

  pa=Points+a;
  pb=Points+b;
  pc=Points+c;
  pd=Points+d;

  if(pc->nngbs > MAX_NGBS_PER_PNT || pd->nngbs > MAX_NGBS_PER_PNT)
    die("\
cannot re-assign edges to obtain Delaunay-hood, since the points to be\n\
connected already have the maximum nr of edges (MAX_NGBS_PER_PNT) per point\n");

  edge->ends[0]= c;
  edge->ends[1]= d;
  edge->length2= DIST3_2(pc->xyz, pd->xyz);

#define NEXT_EDGE(X, O) X->edges[ (i + X->nngbs + O) % X->nngbs ]

#define REDO_EDGES1(P, L, R)						   \
  for (j=0,i=0; i<P->nngbs; i++) {	/* build new ring of neighbours */ \
    if (P->edges[i] == e) {		/* old edge: leave out */	   \
      L=NEXT_EDGE(P, -1);		/* keep for later */		   \
      R=NEXT_EDGE(P, 1 );						   \
    }									   \
    else								   \
      new_edges[j++] = P->edges[i];					   \
  }									   \
  P->nngbs--;								   \
  assert(j==P->nngbs);							   \
  memcpy(P->edges, new_edges, (P->nngbs)*sizeof(SHORT)); 

#define REDO_EDGES2(P, L, R)						   \
  for (j=0,i=0; i<P->nngbs; i++) {	/* build new ring of neighbours */ \
    new_edges[j++] = P->edges[i];					   \
    if (P->edges[i] == L) {		/* place where new edges comes */  \
      assert( NEXT_EDGE(P, 1 )==R);					   \
      new_edges[j++]=e;		/* the new edge */			   \
    }									   \
  }									   \
  P->nngbs++;								   \
  assert(j==P->nngbs);							   \
  memcpy(P->edges, new_edges, (P->nngbs)*sizeof(SHORT));

  /* delete edge e from ngbs list of a and b: */
  REDO_EDGES1(pa, al, ar);		/* al & ar are edges out of a */
  REDO_EDGES1(pb, bl, br);		/* likewise for b */

  /* insert edge e into ngbs list of c and d: */
  REDO_EDGES2(pc, ar, bl);
  REDO_EDGES2(pd, br, al);

  /* reset whoamis, as they may be needed in rest during rest of checks: */

#define RESET_WHOAMIS(P)			\
  pp=Points+P;					\
  edge_array=pp->edges;				\
  for (j=0; j<pp->nngbs; j++) {			\
    ep=Edges + edge_array[j];			\
    if (P==ep->ends[0])				\
      ep->whoami[0]=j;				\
    else					\
      ep->whoami[1]=j;				\
  }

  RESET_WHOAMIS(a);
  RESET_WHOAMIS(b);
  RESET_WHOAMIS(c);
  RESET_WHOAMIS(d);

  return MIN(MIN(al,ar),MIN(bl,br));	/* return the lowest edge involved*/
} /* replace_edge */

#ifdef undefined
static void boxify_points(sphere_p sphere) {
  /*
   * does    : allocates a grid that contains the points in a grid, and
   *   sort points into a linked list rooted in the grid cell
   *
   * gets    : sphere
   *
   * affects : sphere->cells, sphere->pts_per_cell
   *
   * returns : nothing
   *
   * warns   : 
   *
   * comment : nr of cells is simply dim*dim*dim
   */
  int i, j, npoints, dim, ncells, icoord[3], n, offset;
  point_p point, *cells;
  float *x;
  double f, spacing;
  uchar * pts_per_cell;

  f=1.0;				/* frank says 0.5 */
  npoints=sphere->npoints;
/*  dim=(int)floor(pow( ((double)npoints*f), 0.333333333333333)); */
  dim=(int)ceil(pow( ((double)npoints*f), 0.333333333333333));
  spacing = 2.0/(double)dim;

  ncells=dim*dim*dim;

  cells=(point_p*)CALLOC(ncells, sizeof(point_p));
  pts_per_cell=(uchar*)CALLOC(ncells, sizeof(uchar));

#define OFFSET(X,Y,Z) dim*(dim*X + Y)+Z

  point=sphere->Points;
  for (i=0; i<npoints; i++, point++) {
    x=point->xyz;
    for (j=0; j<3; j++) {		/* get integer coordinates */
      f=(x[j] + 1.0)/spacing;
      n=(int)floor(f);
      assert(n<dim);
      icoord[j]=n;
    }
    offset = OFFSET(icoord[0],icoord[1],icoord[2]);
    assert(offset < ncells);
    point->next=cells[offset];
    cells[offset] = point;
    pts_per_cell[offset]++;
  }

  sphere->dim=dim;
  sphere->ncells=ncells;
  sphere->cells=cells;
  sphere->pts_per_cell=pts_per_cell;
  return;
} /* boxify_points */
#endif /* undefined */
