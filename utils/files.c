/*  
 *   copyright/warranty notice:
 *  
 *   Although this software is NOT part of GNU, nor copyrighted by the
 *   Free Software Foundation, Inc., (FSF) and also should not be regarded
 *   as such, I declare the "GNU General Public Licence" of the FSF.,
 *   FULLY APPLICABLE to the code in this file. The GNU General Public
 *   Licence can be found in any GNU product, or obtained from the FSF,
 *   675 Mass Ave, Cambridge, MA 02139, USA.
 *   
 *   This means that all of the code included here must be used, modified,
 *   and redistributed completely for free. This notice itself shall
 *   remain intact verbatimly.  There is absolutely no waranty that this
 *   software works, or is fit for its purpose.
 *   
 *   
 *   Software written by Philip lijnzaad@embl-heidelberg.de
 *   
 *   end of copyright/warranty notice
 *    
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>

/* #include <sys/types.h>			/# not needed #/
 * #include <sys/stat.h>
 */

#include "files.h"

#include "utils.h" 
#include "alloc.h" 

#define MAXNFILES 4 /* includes stdin and out */

/* global variable */
char * files_fullfilename;

/* general purpose variables in this file, for speed and ease */
static int i,j;

typedef struct filerec {
  int lnr; /* linenr */
  char line[MAXLINLEN+1];
  char pb, pipe;		/* pipe indicates popen()ed file (pipe-open) */
  FILE * filedes;
} filerec_t;
typedef filerec_t * filerec_p; 

static filerec_t filelist[MAXNFILES];
static filerec_p fr,lastfile;
static int nfiles;

#define FILNAMERR "file name longer than MAXFILNAMLEN"
#define LINLENERR "line longer than MAXLINLEN"
#define FOPENFAIL "old_findfile: could not find file %s for opening as \"%s\", \nin path '[$]%s'\n", inname, inmode, path
#define POPENFAIL  "old_findfile: file '%s' found, but filter '%s' went wrong", foundname, filter

static void addfile(FILE * file, char * pipe);

FILE * old_findfile(char * path, char * inname, char * inmode) {
/* just calls findfilename() to see wether file exist (the tricky bit), then
 * opens it and wires it into linked list of known files
 */
  static  char command[5 + MAXFILNAMLEN + 1];
  char name[MAXFILNAMLEN + 1];
  char * filter=NULL;
#define MODE "r"			/* @@@ for now @@@ */
  char * foundname, *mode=MODE;
  FILE * file;

  if (!(path&&inname&&inmode))
    ABORT("old_findfile: NULL argument passed");

  if(inmode[0]!='r')
    ABORT("old_findfile: as yet only reading allowed");
  if(inmode[1])
    filter=inmode+1;			/* NULL otherwise, 'cos of init. */

  if (strcmp(inname, "<")==0)		/* stdin */
    file=stdin;
  else { 
    strcpy(name, inname);		/* make copy so we can mess around */
    files_fullfilename=
      foundname=findfile(path, name, mode, NULL); /* locate file */
    if(foundname) {
      if (filter) {			/* popen it */
	sprintf(command, "%s < %s", foundname, filter);
	file=POPEN(command, "r");
	if(!file) {
	  WARNING(POPENFAIL);
	  return(NULL);
	}
      }
      else {				/* fopen it */
	file=fopen(foundname, mode);	/* testing already done so dont */
	rewind(file);
      }
    }					/* if no filename found */
    else { 
      WARNING(FOPENFAIL);
      return(NULL);
    }
  }
  /* now we have an open file descriptor (possibly stdin) */
  addfile(file, filter);		/* even if NULL */
  return(file);				/* even if NULL */
}



filerec_p getfilerec(FILE * file) {
  /* returns pointer to corresponding filerecord; this way we still can use */
  /* the normal file descriptors. Returns NULL if the filerec_p not found */
  static filerec_p lastfile;
  if(lastfile)  /* first invocation misses it */
    if ((lastfile->filedes)==file) /* for quicker acces */
      return(lastfile);
  for (i=0; i<MAXNFILES; i++)  {
    fr=filelist+i;
    if ( fr->filedes == file) {
      lastfile=fr; /* for quicker acces */
      return(fr);
    }
  }
  return(NULL); /* must be NULL first time */
}

static void addfile(FILE * file, char * pipe) { /* add file in list */
  if (nfiles++ >= MAXNFILES)
    ABORT("trying to open more than MAXNFILES");
  if (getfilerec(file)) /* should be the first time NULL */
    ABORT("same file pointer (still) in list");  
  for (i=0; i< MAXNFILES; i++) {
    fr=filelist+i;
    if(!(fr->filedes)) /*  if you find one that's NULL */
      break;
  }
  if (i==MAXNFILES)
    ABORT("Internal error: no empty slots but nfiles < MAXNFILES");
  fr->filedes=file;
  fr->pipe=pipe?'\1':'\0';
  lastfile=fr;
  return;
}
  

FILE * mypopen(const char * command, const char * mode ) {
  /* popen does not always return a NULL string if something failed; (bug?) */
  /* therefore this one. Also check if legal mode */ 
  FILE * file;
  int c=0;

  file=popen(command, mode);		/* may give pointer also if failed */
  c=fgetc(file);

  if(c==EOF) {				/* so test it like this */
    pclose(file);
    return(NULL);
  }
  ungetc(c, file);			/* should work ! */

  return file;
}

void  closefile(FILE * file) {
  fr=getfilerec(file);
  if (fr->pipe)
    pclose(file);
  else
    fclose(file);
  /*  makezero((char *)fr, sizeof filerec_t); */
  fr->lnr=
    fr->line[0]=
      fr->pb=0;
  fr->filedes=NULL;
  nfiles--;
  return;
}

int fgetl(char * string, FILE *stream){
  fr=getfilerec( stream);
  if(!fr)
    ABORT("fgetl: stream unknown");
  if (fr->pb)
    { fr->pb=0;
      fr->lnr++;
      return(lstrcpy(string, fr->line));
    }
  else {
    for (i=0; (j=fgetc(stream)) !=EOF && (j != '\n'); i++)
      string[i]=j;
    if (i > MAXLINLEN)
      ABORT(LINLENERR);
    string[i]='\0';
    fr->lnr++;
    if (j == EOF)
      return (EOF);
    else return (i);
  }
}

void fungetl(char * string, FILE * stream) {
  fr=getfilerec(stream);
  if(!fr)
    ABORT("fgetl: stream unknown");
  fr->pb=1;
  fr->lnr--;
  if (lstrcpy(fr->line, string) > MAXLINLEN)
    ABORT(LINLENERR) ;
}

void fnextl(FILE * stream) {
  fr=getfilerec(stream);
  if(!fr)
    ABORT("fgetl: stream unknown");  
  if (fr->pb)
    fr->pb=0;
  else {
    while ((i=fgetc(stream)) != '\n' && i != EOF)
      ;
    /* fgetc(stream); */
  }
  (fr->lnr)++;
}

int linenr(FILE *file) {
  fr=getfilerec(file);
  if(!fr)
    ABORT("fgetl: stream unknown");
  return(fr->lnr);
}

const char * currentline(FILE * file) {
  char dummy[MAXLINLEN+1];
  fr=getfilerec(file);
  if(!fr)
    ABORT("fgetl: stream unknown");  
  if(fr->pb)
    return(fr->line);
  else {
    (void) fgetl(dummy, file);
    fungetl(dummy, file);
  }
  return(fr->line);
}

int readpast(FILE * infile, char * line, const char * match) {
  int i,j;
  while( (i=fgetl(line, infile)) != EOF) {
    j=strlen(match);
    if (strncmp(line, match, j)==0)
      break;
  }
  if(i==EOF)
    return(EOF);
  else
    return(strlen(line));
}
