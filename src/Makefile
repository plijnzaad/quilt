### Makefile for quilt
### $Id$
### Philip lijnzaad@ebi.ac.uk, 22-7-96

### This Makefile require GNU Make.

include $(BUILD_HOME)/makefile.include

### following should be set in environment, but can also do it here
# HOME					# obvious
# SYS					# uname (1)
# src					# where source code resides
# obj					# where object files reside
# bin					# where exectables reside

### compiler flags:
### all flags given are gcc flags; have to twiddle things yourself if you're
### not using  gcc

DEFINES = -D$(SYS)			# announce system by a '#define'

### location of auxiliarly .h files and libraries:
utilsdir = ../utils			# may have to be created as symlink
# utilsdir = .
utilincdir = $(utilsdir)		
utillibdir = $(utilsdir)
utillibs = -lutils

incdirs = -I$(utilincdir)
libdirs = -L$(utillibdir)
libs = $(utillibs)			#libraries to be linked against

CFLAGS = -g -Wall $(incdirs) $(DEFINES)    # -Wall: gcc's maximal nr. of warnings

### highest optimization flags (only for most number-crunchy modules)
optim = $(incdirs) -g -O2 -ffast-math -finline-functions -funroll-all-loops -DNOSAFE -DNDEBUG -UDEBUG -Wall

utils_h = pdb.h utils.h list.h alloc.h hash.h vecmat.h
# utils_c = getargs.c utils.c alloc.c list.c pdb.c extra-math.c jacobi.c \
# 	hash.c compat.c
# utils_o =  $(utils_c:.c=.o)

quilt_h = quilt.h
quilt_c = main.c acc.c templates.c io.c misc.c top.c conn.c stitch.c \
		seam.c surf.c randomize.c 

quilt_o = $(quilt_c:.c=.o)		# all quilt objects files

### The name for the template file of spherical point distributions:
SPHERES_FILE = .quilt-spheres.$(SYS)

targets = all install help clean utils quilt_o debugmake

.PHONY: all install help clean utils quilt_o debugmake

all: quilt

help:
	@echo "Usage: make [ one of: $(targets)  ]"

install: quilt
	$(INSTALL) quilt $(bindir)

# may have to be created as a symlink
utils:
	(cd $(utilsdir); $(MAKE) )

$(quilt_o): $(quilt_h)			# forget about utils_h dependcy
quilt_o: $(quilt_h) $(quilt_o)

quilt: $(quilt_o) utils
	$(CC) -o $@ $(quilt_o) $(libdirs) $(libs) -lm

# rule for helping debug this Makefile
debugmake:
	echo $(SYS)


### 
### following entries deviate from make's rules, and are included separately
### 
main.o: main.args			#additional dependency

acc.o: acc.c
	$(CC) $(optim) -c $< -o $@	#different command
### accessability routines. Must be compiled with highest optimization

io.o: io.c
### i/o routines
	$(RM) -f $(SPHERES_FILE)
	$(CC) $(CFLAGS) -DTEMPLATES_FILE=\"$(SPHERES_FILE)\" -c $< -o $@

templates.o: templates.c
### templates of spherical point distributions
	$(RM) -f $(SPHERES_FILE)
	$(CC) $(CFLAGS) -c $< -o $@
### (delete spheres.bin, since format may have changed)


clean:
	$(RM) -f quilt *.o $(SPHERES_FILE)
