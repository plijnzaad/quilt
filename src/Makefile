### Makefile for quilt
### $Id$
### Philip lijnzaad@ebi.ac.uk, 22-7-96

### This Makefile require GNU Make.

include $(BUILD_HOME)/make.include

# force error, since not ready yet:
$(shell asdfa	;lkjq13lrj;@ `)

### following should be set in environment, but can also do it here
# HOME					# obvious
# SYS					# uname (1)
# src					# where source code resides
# obj					# where object files reside
# bin					# where exectables reside

# variables pointing to various utilities:
MAKE = /usr/local/bin/make		# gnu make
SHELL = /bin/sh
CC = gcc
RM = rm
INSTALL = install


### compiler flags:
### all flags given are gcc flags; have to twiddle things yourself if you're
### not using  gcc

DEFINES = -D$(SYS)			# announce system by a '#define'

### where to find .h files:
INCDIRS = 
CFLAGS = -g -Wall $(INCDIRS) $(DEFINES)    # -Wall: gcc's maximal nr. of warnings

### highest optimization flags (only for most number-crunchy modules)
OPTIM = $(INCDIRS) -g -O2 -ffast-math -finline-functions -funroll-all-loops -DNOSAFE -DNDEBUG -UDEBUG -Wall

support_h = pdb.h utils.h list.h alloc.h hash.h vecmat.h
support_c = getargs.c utils.c alloc.c list.c pdb.c extra-math.c jacobi.c \
	hash.c compat.c
support_o =  $(support_c:.c=.o)

quilt_h = quilt.h
quilt_c = main.c acc.c templates.c io.c misc.c top.c conn.c stitch.c \
		seam.c surf.c randomize.c 

quilt_o = $(quilt_c:.c=.o)		# all quilt objects files

### The name for the template file of spherical point distributions:
SPHERES_FILE = .quilt-spheres.$(SYS)

targets = all install help clean support_o quilt_o debug

.PHONY: all install help clean support_o quilt_o debug

all: quilt

help:
	@echo Usage: make $(targets)

makedebug:
	echo $(support_o)

install: quilt
	$(INSTALL) quilt $(bindir)

debug:
	echo $(SYS)

$(support_o): $(support_h)
support_o: $(support_o)

$(quilt_o): $(support_h) $(quilt_h)
quilt_o: $(support_h) $(quilt_h) $(quilt_o)
quilt: $(support_o) $(quilt_o) 
	$(CC) -o $@ $(LDFLAGS) $(support_o) $(quilt_o)  -lm


### 
### following entries deviate from make's rules, and are included separately
### 
main.o: main.args getargs.h		#additional dependencies

acc.o: acc.c
	$(CC) $(OPTIM) -c $< -o $@	#different command
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
