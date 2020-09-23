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

#ifndef _FILES_H_
#define _FILES_H_ 1


/* #define MAXLINLEN 512 is in basics.h */
#define MAXFILNAMLEN 132 /*  FILENAME_MAX is too large */
/* shared data */
extern char * files_fullfilename;
/* full pathname of the file most recently opened with a succesfull findfile */

/* prototypes */

#define STRIPBLANK  "regrep -v '^[ \t]*$'" /* egrep faster than grep */
#define CUTCOMMENT  "rcut -d# -f1"
#define UNCOMMENT   "rsed 's/[#].*$//; /^[ \t]*$/d'"
#define UNCOMPRESS  "runcompress"

FILE * old_findfile(char *path, char *filename, char *mode);
/* WHEN READING (as yet only "r" mode allowed): if path is an environment
 * variable, it is interpreted as such, otherwise it will be taken as is;
 * likewise for filename. If path is a colon-separated list of pathnames, it
 * will look through this list of pathnames for finding filename.  If path=="",
 * findfile looks only for file in current dir, subsequently $HOME. 
 * As for now, mode should be "rX", where X is empty, or a command that would
 * evaluate to a possibly compound filter. If X is not empty, the file
 * will be opened as a pipe. Above #define's may serve as useful macros
 * to this end.

   (FOLLOWING NOT YET READY)

 * WHEN WRITING (as yet only "w" and "wZ" allowed): if path is an environment
 * variable, it is interpreted as such, otherwise it will be taken as is;
 * likewise for filename. path should not be a colon-separated pathname-list,
 * just one pathname or env. var. that evaluates to that. It tries to open
 * ($)path/($)filename for writing; if this fails, a message is issued and
 * current dir and $HOME are tried. If this also fails, another message is
 * issued, and NULL returned. A 'Z' in mode makes it open files with
 * popen("compress > ($)path/($)filename.Z", "w");
 */
#define TELLWHICHFILE \
fprintf(stderr, "reading %s\n", files_fullfilename)

#define FINDFILE_OR_Z(FP, DIR,NAME) ( ((FP)=old_findfile((DIR), (NAME), "r"))\
       || (strcat((NAME),".Z"),(FP)=old_findfile((DIR), NAME, UNCOMPRESS)))

void closefile(FILE * file);		/* goes with findfile() */

int fgetl(char * string, FILE * stream); /* and */
void fungetl(char * string, FILE * stream);
/* like fgets() with following modifications:
   - fgetl() only works on files that have been opened with
     old_findfile( , , "r") 
   - fgetl() gets input from either file internal pushback buffer
   - fungetl() effectively pushes back string into a per file internal buffer
   - fgetl() always reads until '\n', there is no length argument
   - fgetl() never stores a '\n'; yet it moves the position
     indicator to the first position after a '\n'
   - returns strlen(string) (may be == 0) or EOF;
   - keeps track of linenr (of the line most recently read); this is available
     as linenr(file);
   - by means of currentline(file) returns corresponding line
   - BUG: if file ends without last \n, the string is normally read, but EOF
     is returned. Therefore, it may be wise to have all files end with a \n .
 in future maybe also a fgetw,fungetw for reading and pushing back words */

extern int linenr(FILE *file);  /* returns linenr of line from file being read */
extern void fnextl(FILE *file); /* moves filepointer to beginning of next line */
extern const char * currentline(FILE * file); /* returns p to internally stored line*/
extern int readpast(FILE * infile, char * line, const char * match);
/* does fgetl(line, infile) until it finds match at beginning of file; 
 * returns the line in line, and length as result, or EOF if not found 
 */
#define POPEN mypopen			/* make it look like macro */
extern FILE * mypopen(const char * string, const char * mode );
/* does a safer popen(); */

#endif /*#ifndef _FILES_H_ */

