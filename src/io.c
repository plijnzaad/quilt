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

/* support routines for in- and output of triangulation */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "extra-math.h"

#include "utils.h"
#include "alloc.h" 
#include "pdb.h"
#include "vecmat.h" 

#include "quilt.h"

/* colors */
#define M_GRAY 		"gray"
#define M_GREY 		"gray"		/* i always mix them up ... */
#define M_YELLOW 	"yellow"
#define M_RED 		"red"
#define M_BLUE 		"10"
#define M_DARK 		"blue"
#define M_WHITE 	"white"
#define M_GREEN 	"green"
#define M_CYAN	 	"cyan"

static void write_tsphere(FILE* file, sphere_p Sphere) {
  /* write single template sphere to file */

  fwrite(Sphere, sizeof(sphere_t), 1, file);
  fwrite(Sphere->coordinates, sizeof(float[3]), Sphere->npoints, file);
  fwrite(Sphere->Points, sizeof(point_t), Sphere->npoints, file);
  fwrite(Sphere->Edges, sizeof(edge_t), Sphere->nedges, file);
  fwrite(Sphere->Faces, sizeof(face_t), Sphere->nfaces, file);
/*   fwrite(Sphere->cells, sizeof(point_p), Sphere->ncells, file);
 *   fwrite(Sphere->pts_per_cell, sizeof(uchar), Sphere->ncells, file);
 */
  return;
} /* write_tsphere */

static void dump_template_spheres(FILE *file) {
  /* writes all templates to a binary file */
  int i;
  
  for (i=0; i< N_TMPLT_SPHERES; i++) 
    write_tsphere(file, Template_spheres[i]);
  return;
} /* dump_template_spheres */

static sphere_p read_tsphere(FILE*file) {
  /* read one binary template sphere from file */
  int i, j;
  sphere_p Sphere;
  point_p Points, Pp, oldPoints /* , *cells */ ;
  edge_p Edges, oldEdges;
  face_p Faces, oldFaces;
  int cdiff, pdiff, ediff, fdiff, nngbs;
/*  uchar * pts_per_cell; */
  float * fp, *coordinates, *oldcoords;

  Sphere=(sphere_p)CALLOC(1, sizeof(sphere_t));
  fread(Sphere, sizeof(sphere_t), 1, file); 

  coordinates=(float *)MALLOC(Sphere->npoints*sizeof(float[3]));
  fread(coordinates, sizeof(float[3]), Sphere->npoints, file);

  Points=(point_p)MALLOC(Sphere->npoints*sizeof(point_t));
  fread(Points, sizeof(point_t), Sphere->npoints, file);

  Edges=(edge_p)MALLOC(Sphere->nedges*sizeof(edge_t));
  fread(Edges, sizeof(edge_t), Sphere->nedges, file);

  Faces=(face_p)MALLOC(Sphere->nfaces*sizeof(face_t));
  fread(Faces, sizeof(face_t), Sphere->nfaces, file);

  oldcoords = Sphere->coordinates;
  oldPoints= Sphere->Points;		/* with which it was written */
  oldEdges = Sphere->Edges;		/* "" */
  oldFaces = Sphere->Faces;

  Sphere->coordinates=coordinates;	/* new addresses */
  Sphere->Points=Points;		
  Sphere->Edges=Edges;
  Sphere->Faces=Faces;

  /* fix pointers, add difference of offsets: */
  cdiff= (int)((ulong)coordinates - (ulong)oldcoords);
  pdiff= (int)((ulong)Points - (ulong)oldPoints);
  ediff= (int)((ulong)Edges - (ulong)oldEdges);
  fdiff= (int)((ulong)Faces - (ulong)oldFaces);

  for (i=0; i<Sphere->npoints; i++) {
    Pp=Points+i;

/*     if (Pp->next)			/# not if it's end of the list #/
 *       Pp->next = (point_p)( (ulong)(Pp->next) + pdiff);
 */
    nngbs=Pp->nngbs;
    Pp->xyz = (float*)((ulong)(Pp->xyz) + cdiff);
    for (j=0; j<nngbs; j++) {
      Pp->Ngbs[j]=(point_p)( (ulong)(Pp->Ngbs[j]) + pdiff);
/*      Pp->Edges[j]=(edge_p)( (ulong)(Pp->Edges[j]) + ediff); */
      fp= Pp->ngbs_xyz[j];
      fp= (float*)( (ulong)fp + cdiff);	/* is same distance (or?) */
      Pp->ngbs_xyz[j]=fp;
    }
  }

/*   cells = CALLOC(Sphere->ncells, sizeof(point_p));
 *   fread(cells, sizeof(point_p), Sphere->ncells, file);
 *   Sphere->cells=cells;
 * 
 *   for (i=0; i<Sphere->ncells; i++) 
 *     if (cells[i])			/# not if end of list #/
 *       cells[i]=(point_p)( (ulong)(cells[i]) + pdiff);
 * 
 * 
 *   pts_per_cell=CALLOC(Sphere->ncells, sizeof(uchar));
 *   fread(pts_per_cell, sizeof(uchar), Sphere->ncells, file);
 *   Sphere->pts_per_cell=pts_per_cell;
 * 
 */
  return Sphere;
} /*  read_tsphere */

static void read_template_spheres(FILE*file) {
  /* read all binary template spheres from file */
  int i;

/*  Template_spheres=(sphere_p*)CALLOC(N_TMPLT_SPHERES, sizeof(sphere_p));
 *  now static */ 

  for (i=0; i<N_TMPLT_SPHERES; i++) 
    Template_spheres[i]=read_tsphere(file);

  return;
} /*  read_template_spheres */

void print_sphere(FILE*file, atom_p atom, char * color) {
  /* prints all edges of atom to file */
  int i, j, npoints, nedges;
  SHORT *pointflags;
  shell_p sh;
  sphere_p sphere;
  point_p Points;
  edge_p Edges, ep;
  float *xyz, *pts_xyz;

  sh=atom->pointer;
  pointflags=sh->pointflags;
  sphere= sh->sphere;
  Points=sphere->Points;		/* template coordinates */
  npoints=sphere->npoints;		/* their number */
  Edges=sphere->Edges;			/* template edges */
  nedges=sphere->nedges;		/* their number */
  
  pts_xyz=(float*)ALLOCA(npoints*sizeof(float[3])); /* scratch area */

  fprintf(file, ".color %s\n", color);
  for (i=0; i<npoints; i++) /* prepare coordinates  */
    for (j=0; j<3; j++) 
      pts_xyz[ 3*i + j] = atom->RADIUS * Points[i].xyz[j] + atom->xyz[j];

  for (j=0; j<nedges; j++) {		/* loop over edges of this atom */
    ep=Edges+j;
    if ( (pointflags[ep->ends[0]] & (BURIED | UNWANTED) ) 
	|| (pointflags[ep->ends[1]] & (BURIED | UNWANTED) ) )
      continue;
    
    xyz = pts_xyz + 3*(Edges[j].ends[0]);
    fprintf(file, ".m %.3f %.3f %.3f\n", xyz[0], xyz[1], xyz[2]);
    xyz = pts_xyz + 3*(Edges[j].ends[1]);
    fprintf(file, ".d %.3f %.3f %.3f\n", xyz[0], xyz[1], xyz[2]);
  } /* loop over edges of this atom */
  fprintf(file, "\n");
  AFREEA(pts_xyz);
  return;
} /* print_sphere */

void write_sphere(char *filename, atom_p atom, char *color) {
  /* outputs all edges of atom to file */
  FILE *file;
  
  file=fopen(filename, "w");
  if(!file) 
    perror(filename),die("\n");
  print_sphere(file, atom, color); 
  fclose(file);
  return;
} /* write_sphere */

void print_cons(FILE* file, atom_p atom, char *color) {
  /* writes all connections starting at atom to filename */
  int j, k, length, p;
  shell_p shell;
  patch_p patch;
  bndry_p bndry;
  SHORT * path;
  list_p l, *cons;
  con_p con;
  float * coordinates, *x, *atomcenter, radius, *direction, r;
  float xa[3], xb[3];

  shell=atom->pointer;
  assert(shell);
  atomcenter=atom->xyz;
  radius=atom->RADIUS;
  coordinates=shell->sphere->coordinates;
  fprintf(file, ".color %s\n", color);
  for(patch=shell->patches; patch; patch=patch->next) 
    for(bndry=patch->bndries; bndry; bndry=bndry->next) {
      path=bndry->path;
      length=bndry->length;
      cons=bndry->cons;
      if (cons) {
	for (j=0; j<length; j++) {
	  p=path[j];
	  x=coordinates + 3*p;		/* unit sphere coordinates of p */
	  for (k=0; k<3; k++)		/* one end of con */
	    xa[k]=atomcenter[k] + radius*x[k];
	  for(l=cons[j]; l; l=l->next) {
	    con=l->ptr;
	    if (con->bndries[0] != bndry) /* only one emanating from here */
	      continue;
	    direction=con->direction;
	    r=direction[3];		/* length of con */
	    for (k=0; k<3; k++) 
	      xb[k]=xa[k] + r*direction[k];
	    fprintf(file, ".m %.3f %.3f %.3f\n.d %.3f %.3f %.3f\n",
		    xa[0], xa[1], xa[2], xb[0], xb[1], xb[2]);
	  } /* for l */
	} /* for j < length */
      } /* if (cons) */
    } /* for patch, bndry */
  return;
} /* print_cons */

void write_cons(char *filename, atom_p atom, char *color) {
  /* writes all connections starting at atom to filename */
  FILE *file;

  file=fopen("w", filename);
  if (!file)
    perror(filename),exit(1);
  print_cons(file,atom, color);
  fclose(file);
  return;
} /* write_cons */

void print_atom(FILE*file, atom_p atom, char *atom_color, char *con_color) {
  /* writes all edges and connections of atom to file */
  print_sphere(file, atom, atom_color);
  print_cons(file, atom, con_color);
  return;
} /* print_atom */

void write_atom(char*filename, atom_p atom, 
		char *atom_color, char *con_color ) { /* not used */
  FILE *file;

  file=fopen(filename, "w");
  if (!file)
    perror(filename),exit(1);
  print_atom(file,atom, atom_color, con_color);
  fclose(file);
  return;
} /* write_atom */

static char*choose_color(atom_p atom) { 
  switch(atom->name[1]) {		/* set atom color: */
  case 'C': return M_WHITE;
  case 'S': return M_YELLOW;
  case 'O': return M_RED;
  case 'N': return M_DARK;
  default:  return M_GREEN;		/* unknown atoms */
  }
} /* choose_color */

void print_triangulines(FILE*file, box_p box) {
  /* writes lines of both atoms and connections to file */
  int i, natoms;
  shell_p sh;
  atom_p atom, atoms;

  natoms=box->natoms;
  atoms=box->atoms;

  for (i=0,atom=atoms; i<natoms; i++,atom++) { /* loop over all atoms */
    sh=(shell_p)atom->pointer;
    if ( (!sh) || (sh->flags & (ISWIDOW | ISWIDOW)))
      continue;
    print_atom(file, atom, choose_color(atom), "gray");
  } /* loop over all atoms */
  return;
} /* print_triangulines */

int print_atom_dots(FILE*file, atom_p atom, char *pre, char *post) { 
  int j,k, npoints, n;
  shell_p sh;
  sphere_p sph;
  float r, *xyz, *pxyz;
  SHORT * pointflags;
  point_p points;

  n=0;
  sh=atom->pointer;
  if ( (!sh) || (sh->flags & ISWIDOW  ))
    return 0;
  xyz = atom->xyz;			/* translation vector */
  sph=sh->sphere;
  npoints=sph->npoints;
  points=sph->Points;
  pointflags=sh->pointflags;
  for (j=0; j<npoints; j++) {
    if (pointflags &&			/* if NULL, atom is ISOLATED */
	pointflags[j] & BURIED)	/* , so should print all */
      continue;
    fprintf(file, "%s", pre);
    pxyz = points[j].xyz;
    r=xyz[3];				/* radius */
    for (k=0; k<3; k++)		/* coordinates */
      fprintf(file, "%10.3f ", (double)(r*pxyz[k] + xyz[k]) );
    fprintf(file, "\n");
    n++;
  } /* for j < npoints */
  return n;
} /* print_atom_dots */

int print_dots(FILE * file, box_p box) {
  int i, n;

  n=0;
  for (i=0; i<box->natoms; i++) { 
    fprintf(file, "# %d\n", i);
    n+=print_atom_dots( file,  &box->atoms[i], ".dot ", "");
  }
  return n;
} /* print_dots */

static int residue_areas(FILE * file, structure_p str) { 
  /* prints totals per residue */
  int i, j;
  residue_p res;
  atom_p at;
  char ch, ss, *atname;
  static char secstr[3]={ '-', 'a','b', };
  float area, total, hfob, hfyl, mainchain, sidechain;


  fprintf(file, "REMARK residue data follows\n");
  fprintf(file, "\
REMARK res chain nr secstr; total, hydrophobic, hydrophilic, main, side\n");
  for (i=0, res=str->residues; i<str->nresidues; i++,res++) {
    ch=res->chain->name; if (ch==' ')ch='-';
    ss =  secstr [ res->flags.SS ] ;
    total=hfob=hfyl=mainchain=sidechain=0.0;
    for (j=0,at=res->atoms; j<res->natoms; j++,at++) {
      area=at->AREA;
      total += area;
      if (strchr("CS", at->name[1]))
	hfob += area;
      else 
	hfyl += area;
      atname=at->name;
      if (atname[0]==' ' && atname[3]==' ' &&
	  (atname[2]==' ' || atname[2]=='A' )   )
	mainchain += area;
      else
	sidechain += area;
    } /* for j */
    fprintf(file, 
	    "REMARK AREA %.4s %c %4hd%c %c %7.2f %7.2f %7.2f %7.2f %7.2f\n",
	    (char*)&res->name[0],		/* not 0-terminated ! */
	    ch, res->nr, res->inscode, ss, total,hfob, hfyl, mainchain,sidechain);
  } /* for i */
  fprintf(file, "REMARK END OF RESIDUE DATA\n");
  return 0;
} /* residue_areas */

void print_areas(FILE*file, box_p box) {
  int i;
  atom_p a;
  shell_p sh;
  structure_p str;

  str=box->structure;

  /* set areas properly: */
  for (i=0,a=str->atoms; i<str->natoms; i++,a++) {
    sh=a->pointer;
    if (sh==NULL)
      a->AREA=0.0;
    else
      a->AREA=SQR(a->RADIUS)*sh->area;
  } /* i  */

  residue_areas(file, str);

  fprintf(file, "REMARK ATOM DATA FOLLOWS\n");
  fprintf(file, "REMARK \n");
  fprintf(file, "REMARK 4th float field is RADIUS, 5th is AREA\n");

  print_pdb(str, file); 

/*   reset ????
 *   for (i=0,a=str->atoms; i<str->natoms, i++,a++) {
 *     sh=a->pointer;
 *     if (sh)
 *       a->AREA  /=SQR(atom->RADIUS)*sh->area; /# reset #/
 *   } /# for i #/
 */  
  return;
} /* print_areas */

static atom_p read_points_file(FILE*file, int *npointsp) {
  /* reads a file with points and radii as atoms. Format: n \n; n *4 %f\n */
  /* Ignore empty and commented lines; lines beyond ^END$ (regexp) are */
  /* also ignored */

  int n=0, natoms;
  char line[133], *cp;
  float dummy, *fp;
  atom_p Atoms=NULL, atom;

  while(fgets(line, 80, file)) {
    if ( (cp=strchr(line, '#')) )		/* ignore stuff beyond '#' */
      *cp = '\0';
    if ( strspn(line, " \t\n") == strlen(line) )
      continue;				/* ignore empty lines */
    if (!strcmp(line, "END\n"))		/* ignore stuff beyond ^END$  */
      break;
    if (! Atoms) {			/* get n of points, and allocate */
      if(sscanf(line, "%d%f", &natoms, &dummy) !=1)
	die("expected 1 integer (number of points) in this line:\n%s\n",
	      line);

      Atoms=(atom_p)MALLOC(natoms*sizeof(atom_t));
    }
    else {				/* fill array */
      atom=Atoms+n;
      fp = atom->xyz;		/* read coordinates */
      if(sscanf(line, "%f %f %f %f", fp+0,  fp+1, fp+2, fp+3) != 4)
	die("expected %d floats in this line:\n%s\n", 4, line);
      n++;
    }
  }
  if (n  != natoms)
    die("Exspect %d points, found %d\n", (int)natoms, n);
  *npointsp = natoms;
  return Atoms;
} /* read_points_file */

void get_template_spheres(FILE* infile) { 
  /* will load templates into Templates[*], either by reading from infile, */
  /* finding a properly named file, or generating them from scratch */
  char *filename;
  FILE *file;

#ifndef TEMPLATES_DIR
#  define TEMPLATES_DIR DFLT_TEMPLATES_DIR
#else
  /* was defined from cc command line to incorporate system dependencies */
#endif
#ifndef TEMPLATES_FILE
#  define TEMPLATES_FILE DFLT_TEMPLATES_FILE
#else
  /* was defined from cc command line to incorporate system dependencies */
#endif

  if (infile)
    file=infile;
  else
    filename=findfile(TEMPLATES_DIR, TEMPLATES_FILE, "r" , &file);

  if(file) {				/* file given or found */
    if(filename)			/* was found itself, so tell */
      warn("reading file %s ...\n", filename );
    read_template_spheres(file);	/* read in to Template_spheres[] */
    fclose(file);
    return;
  }

  warn("template spheres not found; generating them ...\n"); 
  generate_template_spheres();	/* are put in Template_spheres[] */
  warn("done; now writing template spheres to %s ... ", TEMPLATES_FILE); 
  file=fopen(TEMPLATES_FILE, "w");
  if(!file)
    perror(TEMPLATES_FILE), die("\n");
  dump_template_spheres(file);
  fclose(file);
  warn("done\n");
  return;
} /* get_template_spheres */

box_p get_atoms(FILE * pointsfile, char* protein, char*bin, uint pdbflags) {
  /*
   * does    : reads in a set of atoms 
   *
   * gets    : a file with points, or the name of a pdb structure or file
   *  (to be read with pdbflags to findstruc)
   *
   * affects : contents of box
   *
   * returns : a box, with box->natoms, box->structure, box->atoms set.
   *
   * warns   : not
   *
   * comment : this is the place where box is allocated
   */
  int n;
  box_p box;
  structure_p str;
  char *cp;

  box = (box_p)CALLOC(1, sizeof(box_t));

  if (pointsfile)
    strcpy(box->name, "file");
  else  { 
    cp= (protein[0]) ? protein : bin;
    n=strlen(cp)+1;
    n= sizeof(box->name) - n;
    if (n < 0) {			/* too long: truncate at front */
      n= -n;
      box->name[0]='*';
      strcpy(box->name+1, cp+n);
    } else { 
      strcpy(box->name, cp);
    }
  }

  if (pointsfile) {			/* we're given a pointsfile */
    box->atoms=read_points_file(pointsfile, &box->natoms);
    box->structure = NULL;		/* there's no structure */
    return box;
  }

/*  assert(protein[0]); */
  if (protein[0]) {
    str=findpdb(protein, pdbflags);
    n=finddssp(protein, str);
    if(n>1)
      warn("error reading dssp file of %s\n", protein);
    /*     ((err==1)?(warn):(die))("error reading dssp file of %s\n", protein);
     */
  } else
    str=findpdb_bin(bin);

  box->structure=str;
  box->atoms = str->atoms;
  box->natoms = str->natoms;
  
  return box;
} /* get_atoms */

