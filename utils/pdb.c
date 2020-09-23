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

/* module to read in pdb files */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "extra-math.h" 
#include <ctype.h>
#include <string.h> 

#include "vecmat.h"			/* macros for vector things */
#include "utils.h"			/* some basic definitions/prototypes */
#include "alloc.h"
#include "hash.h"

#include "pdb.h"			/* the data structure. */

/* about the worst you can get in a brookhaven/pdb file:
 * 
ATOM    436 1DG2ATHR 21339B     16.597  49.096 -14.477  0.50 17.98
01234567890123456789012345678901234567890123456789012345678901234567890
0         1         2         3         4         5         6         7
 * 
 * atom nr 436 is first deuterium on heavy atom *G2, altloc A, on residue
 *  Threonine 1339B of chain '2'
 */

/* constants: */
#define MAXRESTYPES 128			/* for fast lookup of templates */


#define MAXNCHAINS 8			/* only when specifying explicitly */
#define LINLEN	 256			/* length of pdb-file lines */
#define NRESUNIQ 10			/* nr of chars that makes res unique */

#define atoms_BLOCK 2048		/* by which to extend allocation */
#define residues_BLOCK 256		/* for these things. SHOULD BE  */
#define chains_BLOCK 4			/* POWERS OF TWO !! */


/* file global stuff:  */
static struct {
  int templates_read;
} local_once;			/* initialized once */

static struct { 
  int category_ignored [ NCATEGORIES ] ; /* ignore_category */
  int model, Q_ig, H_ig, altloc_ig	/* check_line */
    , bb_order;				/* check_backbone */
  char residue_id[ NRESUNIQ ];
} local;


int CatNrOf(const char *sym) {
/* a_acid
 * n_acid
 * cofactor
 * ion
 * solvent
 */
  switch(sym[0]) {
  case 'a': return A_ACID;
  case 'n': return N_ACID;
  case 'c': return COFACTOR;
  case 'i': return ION;
  case 's': return SOLVENT;
  default: return -1;
  }
} /* CatNrOf */

static hnode_p make_entry(htable_p table, char *name, 
			  int type, int category, int natoms) {
  hnode_p node;
  residue_template_p res;

  node=Hget(table, name);
  if (node)				/* already defined */
    return node; 

  node=Hset(table, name);		/* create new entry */

  res=(residue_template_p)MALLOC(sizeof(residue_template_t));
  node->value=res;

  strcpy(res->name, name);
  res->type = type;
  res->category=category;
  res->natoms=natoms;
  return NULL;				/* here means: succes */
} /* make_entry */

static htable_p residue_types(FILE*file, char*filename, char *line, 
			 htable_p atom_data /* not used */ ) { 
  /*
   * does    : reads residue templates from file into hash table, and
   *   return returns it
   *
   * gets    : open file, filename for error messages, scratch line to
   *   use, hash table with atom data
   *
   * affects : ResidueTemplates is allocated and set here
   *
   * returns : nr of residue types found
   *
   * warns   : 
   *
   * comment : the hash table is indexed by residue name (0-terminated),
   *  and has (pointers to) struct residue_template' as value. 
   */

  int linenr, l, nfound;
  int type, cat, natoms;
  char name[RESNAMLEN+1], catname[12];
  char *cp;
  htable_p residue_data;

#define GETLINE								\
  if(fgets(line, LINLEN, file)==NULL) /* end of file */			\
    break;								\
  linenr++;								\
  if((cp=strchr(line, '#')))		/* get rid of comment */	\
    *cp='\0';								\
  l=strlen(line); if(l > LINLEN)					\
    die("%s line %d: line length longer than %d\n",			\
	filename, linenr, LINLEN);					\
  if(strspn(line, " \t\n")==l )						\
    continue

#undef ERROR
#define ERROR(S) die("%s (%s line %d):\n%s\n", S, filename, linenr, line)

  residue_data=Hcreate(NULL, 50, 0);
  linenr=0;
  nfound=0;				/* nr of types found so far */
  while(1) {
    GETLINE;
    
    /* scan residue information */
    if ( sscanf(line, " \"%[a-zA-Z0-9_ ]\" %d %s %d", 
		name, &type, catname, &natoms) < 4 )
      ERROR("exspecting residue specification: name type-nr category natoms");

    if  ( type > MAXRESTYPES )
      ERROR("type nr > MAXRESTYPES");

    cat=CatNrOf(catname);
    if (cat < 0)
      ERROR("unknown category");

    if ( make_entry(residue_data, name, type, cat, natoms) )
      die("residue '%s' used more than once\n", name);
      
  } /* while (1) */
  make_entry(residue_data, "UNK", 'X', A_ACID, 0);
  make_entry(residue_data, "XXX", 'x', COFACTOR, 0);
  return residue_data;
} /* residue_types */


#define ATOMDATA_ELTS 3
static void * parse_atom_data(const char * string) { 
  int i;
  float values[ ATOMDATA_ELTS ];
  
  bzero((char*)values, sizeof(values));
  i=sscanf(string, "%f %f %f", values, values+1, values+2 );
			  /* radius charge solvation-energy */
  if (i < 2  )				/* may omit third value */
    return NULL;
  
  return memsav(values, sizeof(float[ ATOMDATA_ELTS ]));
} /* parse_atom_data */

static htable_p atom_types(FILE * file, char *filename, char *line) {
  /*
   * does    : read in atom types  
   *
   * gets    : an opened file, a filename (for the error messages), and
   *   storage to read the lines into
   *   
   * affects : sets pointer arguments with hash tbale. line is used
   *
   * returns : nr of atom types found
   *
   */

  htable_p atom_data;
  atom_data = Hcreate(NULL,		/* don't load other table */
		      50,		/* nr of buckets */
		      0);		/* keys are strings: atom names */
  Hread(file, atom_data, "%s ", "%[-+ \t0-9.e]", 0, &parse_atom_data);
  return atom_data;
  /* that's all folks */

#ifdef undefined
  while(fgets(line, LINLEN, file)) {
    if((cp=strchr(line, '#')))		/* get rid of comment */
      *cp='\0';
    l=strlen(line); 
    if(l > LINLEN)
      die("%s line %d: line length longer than %d\n",
	  filename, linenr, LINLEN);
    if(strspn(line, " \t\n")==l )	/* emty line: skip  */
      continue;

    if(sscanf(line, "%s %s %s" , name, &s_radius[0], 
	      &s_charge[0], ) != 3)
      die("%s\n: exspected <name> <radius> <charge> (%s line %d)\n",
	 line, filename, linenr);

    if(Hget(atom_data, name))
      die("defining atom type %s multiply (%s line %d)\n",
	  name, filename, linenr);

    values[0]=strtod(s_radius, &cp);
    if (cp[0]) {			/* found a symbol, not a value */
      node=Hget(atom_data, s_radius); /* address of value of hash entry  */
      if(node)
	values[0] = ((float*)node->value)[0];
      else
	die("'%s' not defined (%s line %d)\n", 
	    name, filename, linenr);
    }
    values[1]=strtod(s_charge, &cp);
    if (cp[0]) {			/* found a symbol, not a value */
      node=Hget(atom_data, s_charge);
      if(node)
	values[1] = ((float*)node->value)[1];
      else
	die("'%s' not defined (%s line %d)\n", 
	    name, filename, linenr);
    }

    /* now store */
    Hset(atom_data, name)->value = memsav(values, sizeof(float[3]));

  } /* while(fgets) */
#endif
} /* atom_types */

static int read_templates(char * path, char*atoms, char*residues) {
  /*
   * does    : reads residue and atom templates
   *
   * gets    : colon separated path to search through, and names of atom
   *   and residue template files
   *
   * affects : AtomTypes and ResidueTypes are set here, and
   *   ResidueTemplates is set to contain pointers into ResidueTemplates
   *   for faster lookup
   *
   * returns : 
   *
   * warns   : about atom and residue types
   *
   * comment : 
   */

  int type;
  char * filename, line[ LINLEN ];
  FILE *file;
  htable_p atom_data, residue_data;
  hnode_p node;
  residue_template_p t;
  
#define LINTOOLONG "length of line longer than LINLEN"

  if (local_once.templates_read++) /* already there */
    return 0;

  filename=findfile(path, atoms, "r", &file);
  if(!filename)
    die("no atom types file\n");
  atom_data=atom_types(file, filename, line);
  fclose(file);
  warn("read %d atom types from %s\n", atom_data->nkeys, filename);

  filename=findfile(path, residues, "r", &file);
  if(!filename)
    die("no residue template file\n");
  residue_data=residue_types(file, filename, line, atom_data);
  fclose(file);
  warn("read %d residue types from %s\n", residue_data->nkeys, filename);
  

  ResidueTemplates=(residue_template_p *)
    CALLOC(MAXRESTYPES, sizeof(residue_template_p));

  /* build fast lookup table  */
  for(node = Hkeys(residue_data); node; node=node->next) {
    t=node->value;
    type=t->type;
    if ( ResidueTemplates[type] )	/* already tested before ? */
      die("type %d defined as %s and as %s\n", type, t->name,
	  ResidueTemplates[type]->name);
    ResidueTemplates[type]=node->value;
  }

  AtomTypes=atom_data;			/* set global variables */
  ResidueTypes=residue_data;

  return 0;
} /* read_templates */


static residue_template_p type_of(const char * res) {
  /*
   * does    : find out residue_template of residue name
   *
   * gets    : start of residue name in pdb line
   *
   * affects : nothing
   *
   * returns : pointer to residue_template 
   *
   * warns   : if unknown residue
   *
   * comment : 
   */

  char name[4];
  hnode_p node;

  memcpy(name, res, 3);
  name[3]='\0';
  
  node=Hget(ResidueTypes, name);
  
  if (!node) {
    node=Hget(ResidueTypes, "XXX");	/* NOT UNK: would be a_acid */
    warn("residue %s unknown, will be considered a COFACTOR\n", name);
  }

  return node->value;
} /* type_of */


static int ignore_category(int cat, uint flags) { 
/*  static int category_ignored [ NCATEGORIES ] ; */
  static char * names_of_category[ NCATEGORIES ] = {
    "amino acid",
    "nucleic acid",
    "cofactor",
    "ion",
    "solvent",
  };

  if ( BIT(cat) & (flags & CAT_IGNORE_MASK) ) {
    if ( ! local.category_ignored[cat]++ )
      warn("%s's ignored\n", names_of_category[cat]);
    return 1;
  }

  return 0;
} /* ignore_category */


#define R_end -1
#define R_ok 0
#define R_new 1
#define R_dont_want 2

static int check_line(const char * line, char *chainids, uint flags, 
		      residue_template_p * templatep) {
  /*
   * does    : reads one line of pdb, and returns NULL if not suitable,
   *   otherwise 'parses' it and returns pointer to static structure
   *
   * gets    : 
   *
   * affects : *templatep is set to the template(pointer) of the new
   *   residue  (also if a new residue)
   *
   * returns : R_dont_want if not wanted, R_end if done, R_new if new
   *   residue, R_ok if still old residue
   *
   * warns   : about every things ignored (one time per 'kind' of thing)
   *
   * comment : very simple, should be either ^ATOM  or ^HETATM, rest
   *   depends on flags
   */
  int mod;
/*   static int model;
 *   static int Q_ig, H_ig, altloc_ig;
 */
  char *residue_id = local.residue_id;
  const char * atnam;
  residue_template_p template;

  mod=0;

  if (   !memcmp(line, "ATOM", 4)	/* ATOM */
      || !memcmp(line, "HETA", 4)	/* HETATM */
      || (mod= !memcmp(line, "MODE", 4)) ) /* MODEL */
    ;					/* line may be interesting */
  else
    return R_dont_want;
  /* could do it faster by making integers out of 4 byte entities, and */
  /* comparing these ... have to be sure they're aligned properly though ... */
  
  if(mod) {				/* having a MODEL record: NMR thing */
/*     sscanf(&line[6], "%d", &nr); 
 *     if (nr==1)				/# first model #/
 */
    if(!local.model++)		/* first model; can also have nr 0 */
      return R_dont_want;
    /* else: other models */
    warn("ignoring further models...\n");
    return R_end;
  }

  if (chainids && strchr(chainids, line[21])==NULL)
    return R_dont_want;

  atnam=&(line[12]);
  if ( memchr(atnam, (uchar)'Q', ATNAMLEN) ) {
    if (!local.Q_ig++)
      warn("ignoring Q's...\n");
    return R_dont_want;
  }

  if ( (flags & IGNORE_HS) && ( atnam[1]=='H' || atnam[1]=='D' ) ) {
    if (!local.H_ig++)
      warn("ignoring H's...\n");
    return R_dont_want;
  }

  if ( line[16] != ' ') {		/* finding an altloc */
/*    atomflags |= ALTLOC;  */
    if  ( ! (line[16] =='A' || line[16] == '1') ) { /* not first of altloc */
      if (!local.altloc_ig++)
	warn("finding altloc's: using first conformer, ignoring rest\n");
      return R_dont_want;
    }
  }

  if (memcmp(residue_id, &line[17], NRESUNIQ)) { /* new residue */
    template=type_of(&line[17]);
    if (ignore_category(template->category, flags)) /* does warnings itself */
      return R_dont_want;

    /* else: seem to want it */
    memcpy(residue_id, &line[17], NRESUNIQ); /* keep its identity */
    *templatep=template;
    return R_new;
  }
  return R_ok;
} /* check_line */

static int enter_residue(residue_p residue, ulong atom_offset, int natoms, 
			 const char * line, 
			 residue_template_p template) { 
  /*
   * does    : enters last natoms atoms as a new residue
   *
   * gets    : residue_p to write to, atom offset of first atom, total nr
   *   of them, and template the residue is supposed to look like
   *
   * affects : *residue
   *
   * returns : deviation of expected nr of atoms
   */

  residue->atoms=(void*)atom_offset;
  residue->natoms=natoms;
  residue->type=(uchar)template->type;

  sscanf(&line[22], "%hd", &residue->nr);
  residue->inscode=line[26];

  residue->name[0]=line[17];
  residue->name[1]=line[18];
  residue->name[2]=line[19];
  residue->name[3]='\0';

  *(uchar*)&residue->flags=(uchar)'\0';
  residue->pointer=NULL;

  if (!natoms) {				/* unknown, no use */
    assert(residue->type=='X');
    return 0;
  }
  return natoms - template->natoms;
} /* enter_residue */

static void enter_chain(chain_p chain, ulong residue_offset, 
			int nresidues, int name, residue_p residues) {
  /*
   * does    : finishes chain read so far
   *
   * gets    : a chain_p to write into, the offset of the first residue,
   *   the nr of residues, and id
   *
   * affects : *chain
   *
   * returns : nothing
   *
   * warns   : 
   *
   * comment : 
   */
  int i;
  char * seq;

  chain->name=(char)name;
  chain->residues=(void*)residue_offset;
  chain->nresidues=nresidues;
  seq=chain->sequence=MALLOC(nresidues+1);
  for (i=0; i<nresidues; i++) 
    seq[i]=residues[residue_offset+i].type;
  seq[nresidues]=0;
  return;
} /* enter_chain */

static void enter_atom(atom_p atom, const char * line) {
  /*
   * does    : enters this line's atom into memory
   *
   * gets    : line to read from
   *
   * affects : *atom is set to whatever it looks like
   *
   * returns : nothing
   *
   * warns   : 
   *
   * comment : 
   */

  sscanf(&(line[30]),"%f %f %f %f %f",
	 &atom->xyz[0], &atom->xyz[1], &atom->xyz[2],
	 &atom->OCCUPANCY, &atom->BFACTOR);
  
  memcpy(atom->name, &line[12], ATNAMLEN);
  atom->name[ ATNAMLEN ]=(char)0;

  /* correct old use of ASN/GLN atoms A[D,E][1,2] */

  if ( line[13]=='A' /*elementname*/ && line[19]=='N' ) 
    switch(line[15]) {
    case '1':  atom->name[1]='O';	/* slightly arbitrary */
      break;
    case '2':  atom->name[1]='N';
      break;
    }
  return;
} /* enter_atom */

/* fill in the bounding box and covariance */
static void bounding_box(structure_p Str, 
                  double * coordstats, double covariance[3][3]) {  
  int i, natoms,j;

  natoms=Str->natoms;

  for (i=0; i<3; i++)
    Str->center[i] = coordstats[i]/ natoms;
  for (i=0; i<3; i++)			/* following maybe not very accurate */
    for (j=0; j<3; j++)
      Str->xyz_covariance[i][j]= covariance[i][j] - 
        coordstats[i]*coordstats[j]/natoms;
  /* this is matrix with elts. sigma{ (Pu - avg(Pu))*(Pv - avg(Pv)) },
     where P are coordinates, and u, v run over {x, y, z}  */
  /* radius of gyration is sqrt(sigma(xyz_covariance[i][i]/natoms)) */
  for (i=0; i<6; i++)
    Str->bounding_box[i] = (coordstats+3)[i];
} /* bounding_box */


#define MAYBE_EXTEND(THING) if ( (n##THING##s & (THING##s_BLOCK-1) ) == 0) \
    THING##s=(THING##_p)AREALLOC(THING##s,				   \
				 (n##THING##s + THING##s_BLOCK)		   \
				 *sizeof(THING##_t));			   \
    THING=THING##s+n##THING##s;						   \
    n##THING##s++;
         
static structure_p relocate(int natoms, int nresidues, int nchains,
		    atom_p atoms, residue_p residues, chain_p chains,
		     double * coordstats, double covariance[3][3]) { 
  /*
   * does    : returns unused memory, allocates a true structure_p, and
   *   sets pointer correspondences
   *
   * gets    : 
   *
   * affects : 
   *
   * returns : allocated structure
   *
   * warns   : dies if no atoms, residues or chains were read
   *
   * comment : 
   */
  int c, r, a;
  structure_p Str;
  chain_p chain;
  residue_p residue;
  atom_p atom;
  char * allchains;

  
  if_not( natoms|nresidues|nchains )
    die("no atoms, residues or chains read\n");
  allchains=(char*)MALLOC( (nchains+1)*sizeof(char) );

  atoms=REALLOC(atoms, natoms*sizeof(atom_t));
  residues=REALLOC(residues, nresidues*sizeof(residue_t));
  chains=REALLOC(chains, nchains*sizeof(chain_t));

  for(c=0, chain=chains; c<nchains; c++, chain++) {
    allchains[c]=chain->name;
    chain->residues=residues+(ulong)chain->residues;
    for(r=0, residue=chain->residues; r<chain->nresidues; r++, residue++) {
      residue->chain=chain;
      residue->atoms=atoms + (ulong)residue->atoms;
      for(a=0, atom=residue->atoms; a < residue->natoms; a++, atom++) {
	atom->residue=residue;
      }	/* for atom */
    } /* for residue */
  } /* for chain */
  allchains[nchains]='\0';

  Str=(structure_p)CALLOC(1, sizeof(structure_t));

  Str->natoms=natoms;
  Str->nresidues=nresidues;
  Str->nchains=nchains;
  Str->atoms=atoms;
  Str->residues=residues;
  Str->chains=chains;
  Str->allchains=allchains;

  bounding_box(Str, coordstats, covariance);

  return Str;
} /* relocate */


/* initialize or update the statistics for the bounding box */
static void add_coord(atom_p atom, double * coordstats, 
		      double covariance[3][3]) {
  int i, j;
  double x, *boundbox=coordstats+3;

  if (atom==NULL) {			/* reset */
    for (i=0; i<3; i++) {
      coordstats[i] = 0.0;
      boundbox[ 2*i ] = HUGE_VAL;
      boundbox[2*i+1] = -HUGE_VAL;
      for (j=0; j<3; j++) 
	covariance[i][j]=0.0;
    }
    return;
  }

  for (i=0; i<3; i++) { 
    x=atom->xyz[i];
    coordstats[i] += x;
    for (j=0; j<3; j++) 
      covariance[i][j] += x*atom->xyz[j];
    if (x < boundbox[2*i]) boundbox[ 2*i ] = x;
    if (x > boundbox[2*i+1]) boundbox[2*i+1] = x;
  }
  return;
} /* add_coord */

static int check_backbone(residue_p res, atom_p atoms) { 
  /* make sure it's N,CA,C,O + rest. return 0 if all OK, 1 if wrong order */
  /* or wrong name, or minus nr of missing backbone atoms */
  /* NOTE: this does not check coordinates, just names !!! */
  ulong offset;
  int i;
  atom_t *bbatom;
  static char *bbnames[4]=  { " N  ", " CA ", " C  ", " O  " };
  
  offset=(ulong)res->atoms;
  for (i=0; i<4; i++) { 
    bbatom=atoms + offset + i;
    if ( strcmp(bbnames[i], bbatom->name) != 0 )
      break;				/* wrong name */
  }
  if (i==4)
    return 0;				/* ok */

/*   NOT FIXING THE ORDER; this needs more intellicence!
 *   for (i=0; i<4; i++) {
 *     bbatom=atoms + offset + i;
 *     strcpy(bbname, bbnames[i]); chardel(bbname, " ");
 *     for (j=0; j<res->natoms; j++) {	/# delete white space #/
 *       atom = atoms + offset + j;
 *       strcpy(name, atom->name); chardel(name, " "); /# delete white space #/
 *       if (strcmp(bbname, name) == 0) {	/# found #/
 * 	if (i != j) {			/# not in proper place: exchange #/
 * 	  memcpy(&tempatom, bbatom,   sizeof(atom_t));
 * 	  memcpy( bbatom,   atom,     sizeof(atom_t));
 * 	  memcpy( atom,     &tempatom, sizeof(atom_t));
 * 	}
 * 	break;
 *       }
 *     } /# for j < natoms #/
 *     if (j==res->natoms)			/# not found #/
 *       return -1;
 *   } /# for i < 4 #/
 * 
 */

  return 1;
} /* check_backbone */

static structure_p readatoms(FILE * file, char * chainids, uint flags) {
  /*
   * does : reads all actual atom data from file. If chainid != NULL or
   * "", only some chains are read ('_' for ' '). Other behaviour
   * specified by flags
   *
   * gets    : open file, possibly a chain id, and flags
   *
   * affects : 
   *
   * returns : allocated structure
   *
   * warns   : about funnies
   *
   * comment : */
  char lines[ 2*LINLEN ], *line, *prevline;
  int which, status, s;
  int linenr, natoms, nresidues, nchains, nresatoms, nchainres, chainname;
  atom_p atom,atoms;
  residue_p residue, residues;
  chain_p chain, chains;
  residue_template_p template, newtemplate;
  double coordstats[9], covariance[3][3]; /* center of mass; xmin,max, etc. */

  linenr=natoms=nresidues=nchains=nresatoms=nchainres=chainname=0;
  lines[ 0 ]=lines[ LINLEN ]=(char)0;
  atoms=NULL;
  residues=NULL;
  chains=NULL;
  template=NULL;
  add_coord(NULL, coordstats, covariance);		/* reset */
  
  line=lines; prevline=lines+LINLEN; which=1;
  while( 1 ) {
    if ( fgets(line, LINLEN, file) ) {	/* if not NULL, ie. eof  */
      linenr++;
      if( (status=check_line(line, chainids, flags, &newtemplate)) 
	 == R_dont_want )
	continue;			/* line & prevline remain same */
    } else
      status = R_end;			/* eof */

/*     if(status==R_end)
 *       break;
 */
    if(status == R_new || status == R_end ) { /* new residue */
      if (natoms) {			/* only if not very first */
	MAYBE_EXTEND(residue);
	s=enter_residue(residue, (uint)natoms - nresatoms, nresatoms, 
			   prevline, template);
	if (s &&			/* too many/few atoms */
	    (char)template->type != 'X' ) /* not unknown residue */
	  warn("residue %s %c%hd%c has %d too %s atoms\n", 
	     residue->name, (char)chainname, residue->nr, residue->inscode,
	     s <0? -s : s, s < 0 ? "few" : "many");

	if ( ResidueTemplates[residue->type]->category == A_ACID ) { 
	  s=check_backbone(residue, atoms);
	  if (s == 1 && ! local.bb_order++)
	    warn("*** wrong order or atom names in backbone\n");
	  if (s < 0) { 
	    residue->flags.bb_missing=1;
	    warn("*** backbone atoms missing for %s\n", residue_name(residue));
	  }
	}
	nchainres++;
      }
      template=newtemplate;
      nresatoms=0;
    }

    if( chainname != (int)line[21] || status== R_end) {	/* finish prev chain */
      if (natoms) {			/* not the very first one */
	MAYBE_EXTEND(chain);
	enter_chain(chain, (uint)nresidues - nchainres, 
		    nchainres, chainname, residues);
      }
      chainname = (int)line[21];
      if( line[25] != '1' && status != R_end )
	warn("first residue of chain '%c' has number %.4s\n", 
	     (char)chainname, &line[22]);
      nchainres=0;
    } /* new chain */

    if (status == R_end)
      break;

    MAYBE_EXTEND(atom);
    enter_atom(atom, line);		/* read it into atoms */
    add_coord(atom, coordstats, covariance);
    nresatoms++;

    /* switch lines: */
    prevline=line;
    line=lines + which*LINLEN;
    which = (!which);
  } /* while (fgets) */

  if(!natoms)
    die("error: no atoms read\n");

  return relocate(natoms, nresidues, nchains,
		  atoms, residues, chains, coordstats, covariance);
} /* readatoms */

static char * chain_names(char *chainspec) {
  /* returns NULL, or translated form of chain specification */
  char *cp;

  if (chainspec == NULL || chainspec[0]=='\0' || strstr(chainspec,"*"))
    return NULL;

  cp=strchr(chainspec, '_');		/* _ means ' '; have to translate */
  if (cp)
    *cp=' ';
  cp=chainspec;
  while(*cp) {
    *cp=toupper(*cp);
    cp++;
  }
  
  return chainspec;
} /* chain_names */

hnode_p find_atomtype(atom_p atom, htable_p table) {
  /*
   * does    : finds ATOM's entry in TABLE; if TABLE==NULL, uses AtomTypes
   *
   * gets    : atom, to look at its name
   *
   * affects : nothing
   *
   * returns : radius
   *
   * warns   : about radius not found (if atom unknown)
   *
   * comment : 
   */
  hnode_p node;
  char elt[3];
  residue_p res;

  elt[0]=toupper(atom->name[0]), elt[1]=toupper(atom->name[1]), elt[2]='\0';
					/* make sure it's 0 terminated */ ;

  if (table==NULL)
    table=AtomTypes;

  node= Hget(table, elt+ (int)(elt[0]==' ') ); /* latter == 0 or 1 */
  if (node)
    return node;

  /* atoms like 'FE  '. But what with wrong formats like 'CB__' (should */
  /* be '_CB_') In case of Calcium / C-alpha, this will be confusing ... */ 
  node=Hget(table, elt);
  if (node) 
    return node;
  
  node=Hget(table, elt+1);	/* funny things like 'AC' meant as C */
  if (node) {			/* seems to exist, but warn */
    res=atom->residue;
    if(res && res->type=='X')	/* else we can trust it */
      warn("atom %s [%s %c%hd%c] treated as if of type '%c'\n",
	   atom->name, res->name, res->chain->name, res->nr, res->inscode, 
	     elt[1]);
    return node;
  }
  node=Hget(table, "x");
  if (!node)
    die("entry 'x' (for unknown atoms) missing from atom type file\n");
  else  { 
    res=atom->residue;
    if (res)
      warn("using atom type `x' for atom `%s' [%s %c%hd%c]\n",
	   atom->name, res->name,
	   res->chain->name, res->nr, res->inscode);
    else 
      warn("using atom type `x' for atom `%s' (unknown residue)\n", 
	   atom->name);
  }
  return node;				/* may be NULL */
} /* find_atomtype */

static float find_radius(atom_p atom, htable_p table) {
  hnode_p node;
  
  node=find_atomtype(atom, table);
  return ((float*)node->value)[0];
} /* find_radius */

static int init(void) { 
  read_templates(TEMPLATE_PATH, ATOMS_FILE, RESIDUES_FILE);
  bzero((char*)&local, sizeof(local)); /* reset all file-globals */
  return 0;
}
     
structure_p readpdb(FILE * file, char *chainspec, uint flags) {
  char * chainids;
  structure_p Str;

  init();
/*  info=readinfo(file); */
  chainids=chain_names(chainspec);
  Str=readatoms(file, chainids, flags);

  warn("read %d chains, %d residues, %d atoms\n", 
       Str->nchains, Str->nresidues, Str->natoms);
  return Str;
} /* readpdb */

static FILE *findstruc(char *name) {
  /*
   * does    : finds file that is readable as a pdb-file
   *
   * gets    : 
   *
   * affects : 
   *
   * returns : opened file
   *
   * warns   : about the file found
   *
   * comment : If name=="", stdin is returned
   */
  char filename[80], *foundname;
  FILE * file=NULL;

  if(name[0]=='\0') {			/* if "", read from stdin */
    file=stdin;
    strcpy(filename, "stdin");
    foundname=filename;
  } else if ( strpbrk(name, "./") ) {	/* check occurrence of filename chars*/
    sprintf(filename, "%s", name);	/* it's just a file name */
    foundname=findfile(PDBPATH, filename, "r", &file);
  } else {     /* a pdb code. Build a file name out of it: */
    sprintf(filename, "%.4s%s", name, PDBSUFFIX);
    foundname=findfile(PDBPATH, filename, "r", &file);
    if(!foundname) {			/* try 2nd suffix */
      sprintf(filename, "%.4s%s", name, PDBSUFFIX2);
      foundname=findfile(PDBPATH, filename, "r", &file);
    }
  }

  if (!file)
    die("structure %s not found\n", name);
  
  warn("reading %s\n", foundname);

  return file;
} /* findstruc */

static residue_p findresnr(chain_p chain, int nr, char inscode) { 
  /* returns residue such that residue->nr==nr */
  int leftOK, rightOK;
  residue_p residues, res, res2, endres;
  static int start=1;			/* most common starting nr */

  residues=chain->residues;
  res= &residues[ nr - start ];		/* obvious guess */
  endres=residues+chain->nresidues;

  if (res->nr == nr 			/* usual case */
      && (res>= residues) && (res <= endres)
      && res->inscode==inscode)		/* if nr==0, residues[-1].nr */
    return res;				/* often will be 0 coincidentally */

  res2=res;

  do {
    leftOK= (res>= residues);
    rightOK= (res2 <= endres);
    res -= leftOK;			/* doesn't go beyond array bounds */
    res2 += rightOK;
    if (res2->nr==nr && res2->inscode==inscode) { 
      start =  nr - (res2-residues);	/* new start */
      return res2;
    }
    if (res->nr==nr && res->inscode==inscode) { 
      start =  nr - (res-residues);	/* new start */
      return res;
    }
  }  while ( leftOK || rightOK );

  return NULL;				/* not found */
} /* findresnr */

int readdssp(FILE * file, structure_p Str) {
  /* reads DSSP secondary structure info from file into Str. returns 0 if */
  /* all's OK, 1 if slightly wrong, 2 if exit prematurely */
  int brknr, n;
  char c, *cp, *allchains, line[ LINLEN ], inscode;
  chain_p chain;
  residue_p res;

  allchains=Str->allchains;

  /* spool to first residue */
  while( fgets(line, LINLEN, file) )
    if(line[2]=='#')
      break;
  if ( feof(file) )
    return 1+warn("Could not locate residue data by '#'\n");

  n=0;
  while( fgets(line, LINLEN, file)) { 
    if(line[13] == '!') /* chain breaks are denoted that way: ingore them */
      continue;
    c=line[11];
    cp=strchr(allchains, c);		/* is c in allchains? */
    if(! cp )				/* chain ID not in string: not read */
      continue;
    chain=Str->chains+ (cp - allchains); /* simple as that */
    if (sscanf(&(line[5]) ,"%d", &brknr) != 1) 
      return 1+warn("Could not locate brknr in dssp file: %s", line);
    inscode=line[10];
    res=findresnr(chain, brknr, inscode);
    if(!res)
      return 1+warn("Could not find residue corresponding to %s\n", line);

    if (islower((int)line[13])) /* cystines (!= cysteines ...) */
      line[13]='C';
    if (res->type=='a')			/* ACE : counted as A_ACID */
      continue;
    if (res->type != (uchar)line[13]) { 
      if (line[13] == 'X') {		/* botched  (eg. 3dfr ...) */
	warn("########################################\n");
	warn("#### HET-group with same residue nr and chain ID id as ");
	warn("residue %s %c%hd%c\n",
	     res->name, res->chain->name, res->nr, res->inscode);
	warn("########################################\n");
	continue;
      } else
	return 1+warn("Residuetype mismatch: \n%sexspected %d %c", 
	  line, res->nr, res->type);
    }
    n++;
    c=line[16];				/* sec. struc. */
    if (c == 'E')
      res->flags.SS=SS_BETA;
    else if (c >= 'G' && c <= 'I')
      res->flags.SS=SS_ALPHA;
/*     else 
 *       res->flags.SS=SS_NONE
 */
  } /* while fgets(file) */
  if (n != Str->nresidues)
    return warn("read dssp of %d residues, exspected %d, for chain%s '%s'\n",
	 n, Str->nresidues, Str->nchains==1?"":"s", Str->allchains);
  return 0;
} /* readdssp */

static FILE *finddssp_file(const char * name) { 
  /* opens dssp file belonging to name */
  char filename[80], *foundname;
  FILE * file;

  if ( strpbrk(name, "./") ) {	/* check occurrence of filename chars*/
    sprintf(filename, "%s", name);	/* it's just a file name */
    foundname=findfile(DSSPPATH, filename, "r", &file);
  } else { 
    sprintf(filename, "%.4s%s", name, DSSPSUFFIX); /* synthesize filename */
    foundname=findfile(DSSPPATH, filename, "r", &file);
  }

  if (!file) { 
    warn("structure %s not found\n", name);
    return NULL;
  }
  
  warn("reading %s\n", foundname);
  return file;
} /* finddssp_file */

int finddssp(const char * protein, structure_p str) { 
  /* finds dssp file, and reads it into structure */
  /* returns 0 if all ok, 1 if not, 2 if premature exit of reading */
  int err;
  FILE * file;

  file=finddssp_file(protein); 
  if (!file) {
    warn("no dssp file for %s\n", protein);
    return 1;
  }
  err=readdssp(file, str);
  fclose(file);
  return err;
} /* finddssp */

structure_p findpdb(const char * inname, uint flags) {
  /*
   * does    : finds and reads a pdb file, with flags controlling the reading
   *  Files are searched for in PDBPATH or in $PDBPATH (colon separated list of
   *   directories), 
   *    (NOT READY: if not found, filenames are also tried as .z, .Z, .gz)
   * gets    : a pdb-code, possibly extended by  -<chainids>, or a simple
   *   filename (if it contains a '.' or '/'). The pdb-filename extension
   *   is presumed to be PDBSUFFIX or PDBSUFFIX2.
   *
   * affects : 
   *
   * returns : a structure
   *
   * warns   : statistics, things ignored, etc.
   *
   * comment : 
   */
  char *chainids=NULL;
  char * name, *cp;
  FILE *file;
  structure_p Str;
  int  l;

  l=strlen(inname);
  MEMSAV(name, inname, l);		/* local copy (on stack!) of inname */

  if (! strpbrk(name, "./")) {		/* no filename; maybe 4letter code*/ 
    cp=strchr(name, PDBCHAINSEP);
    if (cp) {				/* one or more chain IDs, eg 4hhb-a */
      chainids=cp+1;
      *cp='\0';				/* truncate name */
    } 
  }

  file=findstruc(name);			/* just name (without chainids) */
  Str = readpdb(file, chainids, flags);
  fclose(file);

  AFREEA(name);
  return(Str);
} /* findpdb */

static int write_bintemplate(const structure_t *str, FILE* file) {
  /* write struc in binary fashion; coordinates are written, but not used  */
  int n=0;
  
  n += fwrite(str, sizeof(structure_t), 1, file);
  n += fwrite(str->chains, sizeof(chain_t), str->nchains, file);
  n += fwrite(str->residues, sizeof(residue_t), str->nresidues, file);
  n += fwrite(str->atoms, sizeof(atom_t), str->natoms, file);

  return n;
} /* write_bintemplate */

static structure_p read_bintemplate(FILE* file) {
  /* read struc;  */
  int i, chaindiff, resdiff, atomdiff;
  char * vp;
  structure_p str;

  str=(structure_p)readfile(file);	/* all in one contiguous go: */
  vp=(char*)str;

  /* adjust the pointers of str: */
  vp += sizeof(structure_t);
  chaindiff= (int)((ulong)vp - (ulong)str->chains) ;
  str->chains=(chain_p)vp;

  vp += str->nchains * sizeof(chain_t);
  resdiff= (int) ((ulong)vp - (ulong)str->residues) ;
  str->residues=(residue_p)vp;

  vp += str->nresidues * sizeof(residue_t);  
  atomdiff= (int) ((ulong)vp - (ulong)str->atoms);
  str->atoms=(atom_p)vp;

  /* and adjust the pointers in all residues */
  LOOPDOWN(i, str->nresidues) { 
      ulong ul=((ulong)str->residues[i].atoms) + atomdiff;
      str->residues[i].atoms = (atom_p)ul;
      ul=((ulong)str->residues[i].chain) + chaindiff;
      str->residues[i].chain = (chain_p)ul;
  }
  
  LOOPDOWN(i, str->natoms) {
      ulong ul=((ulong)str->atoms[i].residue) + resdiff;
      str->atoms[i].residue = (residue_p)ul ;
  }

  LOOPDOWN(i, str->nchains) {
      ulong ul =((ulong)str->chains[i].residues) + resdiff;
      str->chains[i].residues = (residue_p)ul;
  }

  return str;
} /* read_bintemplate */

#define SIZE (sizeof( float[ BINDATA_ELTS ]))

static int write_bincoords(const structure_t *str, FILE* file) {
  /* write str->natoms * sizeof(float[BINDATA_ELTS]) from file */
  int i,natoms;
  float * coords;

  natoms=str->natoms;
  coords=(float*)ALLOCA(natoms*SIZE);

  for (i=0; i<natoms; i++) 
    memcpy( &coords[3*i], str->atoms[i].xyz, SIZE);

  fwrite(coords, SIZE, natoms, file);

  AFREEA(coords);
  return 0;
} /* write_bincoords */

static int read_bincoords(FILE* file, structure_p str) {
  /* read str->natoms * sizeof(float[BINDATA_ELTS]) from file */
  int i, natoms, n;
  float * coords;

  double coordstats[9], covariance[3][3]; /* center of mass; xmin,max, etc. */
  add_coord(NULL, coordstats, covariance);		/* reset */


  natoms=str->natoms;
  coords=(float*)ALLOCA(natoms*SIZE);
  n=fread(coords, SIZE, natoms, file);
  if (n != natoms) {
    WARN("read only %d of %d atom coordinates, ", n, natoms );
    return 1;
  }
  if (ferror(file)) { 
    WARN("Error during read, ");
    return 1;
  }

  for (i=0; i<natoms; i++) { 
    atom_p atom = &str->atoms[i];
    memcpy(atom->xyz, &coords[ BINDATA_ELTS *i], SIZE);
    add_coord(atom, coordstats, covariance);
  }

  bounding_box(str, coordstats, covariance);

  AFREEA(coords);
  return 0;
} /* read_bincoords */

#undef SIZE

int print_pdb_bin(const structure_t  *str, const char*spec) { 
  /* given a SPEC like TEMPLATE%COORDS, writes STR to <template>.<TPL_EXT>
   * and <frame>.<COORD_EXT> in binary fashion 
   */

  char *copy, *cp, filename[30];
  FILE * file;

  MEMSAV(copy, spec, 0);
  cp=strchr(copy, BINSEP);
  if (!cp)
    die("no '%c' separating '%s' into TEMPLATE."TPL_EXT 
	" and COORDINATES."COORD_EXT"\n", BINSEP, spec);
  cp++;
  cp[-1]='\0';

  sprintf(filename, "%s.%s", copy, TPL_EXT);
  file=fopen(filename, "w");
  if(!file) {
    perror(filename);
    return 1;
  } 
  write_bintemplate(str, file);
  fclose(file);
  
  sprintf(filename, "%s.%s", cp, COORD_EXT);
  file=fopen(filename, "w");
  if(!file) {
    perror(filename);
    return 1;
  }
  write_bincoords(str, file);
  fclose(file);
  
  AFREEA(copy);
  return 0;
} /* print_pdb_bin */

structure_p findpdb_bin(const char *spec) {
  /*
   * does    : reads in a protein structure from two binary files,
   *   indicated by SPEC, under control of FLAGS
   *
   * gets    : a SPEC that looks like TEMPLATE%COORDINATES
   *
   * affects : 
   *
   * returns : newly allocated structure
   *
   * warns   : about names of files and statistics on contents
   *
   * comment : 
   */

  char*copy, *cp, *name, filename[30];
  structure_p str;
  FILE *file;

  init();
  MEMSAV(copy, spec, 0);
  cp=strchr(copy, BINSEP);
  if (!cp)
    die("%s cannot be separated in template and coordinate file\n", spec);
  cp++;
  cp[-1]='\0';

  sprintf(filename, "%s.%s", copy, TPL_EXT);
  name=findfile(PDBPATH, filename, "r", &file);
  if(!name) {
    perror(filename);
    return NULL;
  } else warn("reading binary template from %s", name);
  str=read_bintemplate(file);
  if(!str) { 
    perror(name);
    return NULL;
  }    
  fclose(file);

  sprintf(filename, "%s.%s", cp, COORD_EXT);  
  name=findfile(PDBPATH, filename, "r", &file);
  if(!name) {
    perror(filename);
    return NULL;
  } else warn(" and coordinates from %s\n", name);
  
  if( read_bincoords(file, str) ) {
    perror(name);
    return NULL;
  }    
  fclose(file);
  
  AFREEA(copy);
  return str;
} /* findpdb_bin */

void set_radii(structure_p Str, htable_p radii, float probe_radius) {
  /*
   * does    : sets radii of atoms. probe_radius is added to the radii 
   *
   * gets    : Str whose atoms to treat, and hash table to lookup radii
   *
   * affects : Str->atom[*]->RADIUS
   *
   * returns : nothing
   *
   * warns   : about radii not found (assigns entry of "unknown" in that case)
   *
   * comment : 
   */

  int i, natoms;
  atom_p atom;
  float r, maxrad= 0.0;

  natoms=Str->natoms;
  for (i=0, atom=Str->atoms; i<natoms; i++, atom++) { 
    r=find_radius(atom, radii) + probe_radius;
    atom->RADIUS=r;
    if (r > maxrad) maxrad=r;
  }
  Str->maxradius=maxrad;
  return;
} /* set_radii */

/* void correct_gromos_order(structure_p str) {
 *   /# evil, this  #/
 *   int i, j, k, n;
 *   char chainid, *card, altloc=' ';	/# not dealt with now #/
 *   atom_p a;
 *   residue_p res;
 *   chain_p chain;
 * 
 *   n=1;					/# atom nr #/
 *   for(i=0,chain=str->chains; i<str->nchains; i++, chain++) {
 *     chainid=chain->name;
 *     for (j=0, res=chain->residues; j<chain->nresidues; j++, res++) { 
 *       if ( ResidueTemplates[res->type]->category==A_ACID ) {
 *         float Cxyz[5];
 *         float Oxyz[5];
 *         memcpy(Cxyz, res->atoms 2*sizeof(Cxyz));
 *       }
 *     } /# for j #/
 *     res--;				/# last residue: needed in TER-card #/
 *     fprintf(file, "TER              %3s %c%4i%c\n", res->name, chainid, 
 * 	    res->nr, res->inscode );
 *   } /# for i #/
 *   return;
 * } /# correct_gromos_order #/
 */

const char * residue_name(residue_p res) { 
  static char resnam[10];
  
  resnam[0]=0;
  if (res == NULL)
    return resnam;			/* empty string */

  sprintf(resnam, "%c%c%hd%c", res->type, res->chain ? 
	  res->chain->name : ' ', res->nr, res->inscode);
  chardel(resnam, " ");		/* throw away all spaces */
  return resnam;
} /* residue_name */

const char * atom_name(atom_p atom) { 
  int i;
  static char atnam[5];
  char *a, *b=&atnam[0];

  for ( a=atom->name,i=0; i<4; i++,a++) 
    if ( *a  != ' ')
      *b++=tolower(*a);
  *b='\0';
  return atnam;
}

const char * full_atom_name(atom_p atom) { 
  static char fullname[20];

  sprintf(fullname, "%s.%s", residue_name(atom->residue), atom_name(atom));

  return fullname;
} /* full_atom_name */

int print_pdb(const structure_t *str, FILE *file) { 
  int i, j, k, n;
  char chainid, *card, altloc=' ';	/* not dealt with now */
  atom_p a;
  residue_p res;
  chain_p chain;

  n=1;					/* atom nr */
  for(i=0,chain=str->chains; i<str->nchains; i++, chain++) {
    chainid=chain->name;
    for (j=0, res=chain->residues; j<chain->nresidues; j++, res++) { 
      card = ResidueTemplates[res->type]->category==A_ACID ? 
	"ATOM  " : "HETATM";
      for (k=0,a=res->atoms; k<res->natoms; k++, a++, n++) {
	fprintf(file, 
		"%6s%5i %4s%c%3s %c%4i%c    %7.3f %7.3f %7.3f %5.2f %5.2f\n",
		card, n, a->name, altloc, res->name, chainid, 
		res->nr, res->inscode,
		a->xyz[0], a->xyz[1],a->xyz[2], a->xyz[3], a->xyz[4]);
      }	/* for k */
    } /* for j */
    res--;				/* last residue: needed in TER-card */
    fprintf(file, "TER              %3s %c%4i%c\n", res->name, chainid, 
	    res->nr, res->inscode );
  } /* for i */
  return 0;
} /* print_pdb */

void free_pdb(structure_p str) {
  /*
   * does    : frees everything about str
   *
   * gets    : a structure
   *
   * comment : does not free pdb.c's internal tables, such as templates
   */

  FREE(str->atoms);
  FREE(str->residues);
  FREE(str->chains);
  FREE(str->allchains);			/* string of all chain ID's */
  FREE(str);
  return;
} /* free_pdb */

void free_pdb_related(void) {
  /*
   * does    : frees pdb.c's internal tables, such as templates
   * comment : does not free any structures
   */
  FREE(ResidueTemplates);		/* the array */
  Hfree(AtomTypes, _free);
  Hfree(ResidueTypes, _free);
  return;
} /* free_pdb_related */

/* extern global variables, data: */
residue_template_p * ResidueTemplates;	/* for fast lookup */
htable_p AtomTypes, ResidueTypes;
