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

#ifndef _PDB_H_
#define _PDB_H_ 1

#include "compat.h"
/* constants & flags: */
#define PDBPATH "PDBPATH"		/* environment variable containing */
					/* colon separated list of dirs to */
					/* look for pdb files */
#define DSSPPATH "DSSPPATH"		/* environment variable containing */
					/* colon separated list of dirs to */
					/* look for pdb files */

#define PDBSUFFIX ".pdb"		/* file extension of pdb-files */
#define PDBSUFFIX2 ".brk"		/* 2nd one to try */

#define PDBCHAINSEP '-'		/* allowed separators for chains */

#define DSSPSUFFIX ".dssp"
#define DSSPSUFFIX2 ".dssp_pre"

/* residue and atom-templates:  */
#define TEMPLATE_PATH ""
#define ATOMS_FILE ".atoms"
#define RESIDUES_FILE ".residues"

#define A_ACID 0
#define N_ACID 1
#define COFACTOR 2
#define ION 3
#define SOLVENT 4

/* flags to readpdb: */
#define IGNORE_A_ACID  	BIT(A_ACID)	/* 0 */
#define IGNORE_N_ACID  	BIT(N_ACID)	/* 1 */
#define IGNORE_COFACTOR BIT(COFACTOR)	/* 2 */
#define IGNORE_ION      BIT(ION)	/* 3 */
#define IGNORE_SOLVENT  BIT(SOLVENT)	/* 4 */
#define IGNORE_HS 	BIT(5)
/* #define SET_RADII 	BIT(6) */
#define SURFACE		BIT(7)
#define ADD_PROBE_RADIUS BIT(8)		/* does what it says to atom radii */

/* #define READ_DSSP	BIT(7) */

#define PDB_DEFAULT_FLAGS ( IGNORE_SOLVENT | IGNORE_HS /* | SET_RADII */  \
			    | IGNORE_ION )

/* selecting atoms by atom->flags: */
#define SELECTED 	BIT(0)		/* general selection flag */
#define CALC_AREA  	BIT(0)

/* internal constants: */

#define NCATEGORIES 5			/* nr of categories */
#define CAT_IGNORE_MASK (BIT(NCATEGORIES)-1)

#define ATNAMLEN    4
#define RESNAMLEN   3 

#include "utils.h"			/* for typedefs of uchar etc. */

typedef struct residue_template_ {
  int type				/* one letter code */
    , category,				/* a_acid, n_acid, cofactor, ion */
    natoms;				/* , solvent */ 
  char name[RESNAMLEN+1];
/*  atom_p Atoms; */
} residue_template_t, * residue_template_p;

/* extern data: */
extern residue_template_p * ResidueTemplates;
#include "hash.h" 
extern htable_p AtomTypes, ResidueTypes;

/* typedefs: */
typedef struct atom_  { 
  void * pointer;			/* general purpose pointer */
  struct residue_ * residue;		/* to which it belongs */
  float xyz[5];				/* x, y, z, occupancy, B-factor */
  char name[ATNAMLEN+1];		/* (B-factor often accessability) */
  char altloc;				/* '\0' if atom coordinates  missing */
  short flags; 
} atom_t, * atom_p;
#define OCCUPANCY xyz[3]
#define BFACTOR xyz[4]
#define RADIUS xyz[3]
#define CHARGE xyz[4]
#define AREA   xyz[4]

typedef struct residue_ {
  void * pointer;
  atom_p atoms;				/* pointer to its atoms */
  struct chain_ * chain;		/* 'parent' of residue */
  char name[RESNAMLEN+1];

  short nr;				/* as found in pdbfile */
  char inscode; 

  uchar natoms;
  uchar type; /* equals the onelettercode in case of amino acid residues ! */

#define SS_NONE  0
#define SS_ALPHA 1
#define SS_BETA  2
  struct { 
    uchar SS:2				/* see constants above,  */
      , bb_missing:1			/* 1 if backbone not complete */
	, rest:5;			/* misc. */
  } flags;
  float phi, psi, area;
} residue_t, * residue_p;

typedef struct chain_ { 
  residue_p residues;
  int nresidues;
  char name;
  char * sequence;
} chain_t, *chain_p;
 
typedef struct structure_ {
  chain_p chains;
  int nchains;
  char name;
  char * allchains;

  int natoms, nresidues;
  residue_p residues;
  atom_p atoms;
  double center[3], bounding_box[6], maxradius, xyz_covariance[3][3];
} structure_t, * structure_p;

/* prototypes */
structure_p findpdb(const char *name, uint flags);
  /* reads protein from name, which is a code (-<chain>) or a filename */
structure_p readpdb(FILE * file, char *chainids, uint flags);
  /* reads protein from file. Use _ for a blank chain, NULL or "*" for all */

int finddssp(const char *name, structure_p str);
int readdssp(FILE * file, structure_p str);

const char * atom_name(atom_p);		/* ca, cb, etc. */
const char * residue_name(residue_p);	/* VA133 (A for chain a) */
const char * full_atom_name(atom_p atom); /* eg. VA133@cg */

int print_pdb(const structure_t *str, FILE *file);

/* binary i/o stuff: */
#define BINSEP '%'			/* separates template-filename from
					 * coordinate-filename */
#define TPL_EXT "tpl"			/* filename extension for template */
#define COORD_EXT "x"			/* filename extension for coordinates*/
#define BINDATA_ELTS 3			/* nr of floats per atom to be
					 * read/written with pdb_*bin() */
structure_p findpdb_bin(const char *spec);
  /* if SPEC is A%B, reads file from binary A.tpl and B.x */
int print_pdb_bin(const structure_t *str, const char *filename);

/* cleaning up: */
void free_pdb(structure_p);		/* one structure */
void free_pdb_related(void);		/* templates, scratch, etc. */

/* misc.: */
void set_radii(structure_p, htable_p atom_radii, float probe_radius);

hnode_p find_atomtype(atom_p atom, htable_p atomtypes);
/* finds ATOM's type in  ATOMTYPES; if latter is NULL, looks in AtomTypes */
#endif /* _PDB_H_ */
