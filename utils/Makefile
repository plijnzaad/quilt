# $Id$

# requires GNU make

# Philip lijnzaad@ebi.ac.uk

include $(BUILD_HOME)/makefile.include

# where to install stuff:

# prefix = /usr/local
# bindir = $(prefix)/bin
# libdir = $(prefix)/lib

### following should be set in environment, but can also do it here
# HOME					# obvious
# SYS					# uname (1)
# src					# where source code resides
# obj					# where object files reside
# bin					# where exectables reside

### compiler flags:
### all flags given are gcc flags; have to twiddle things yourself if you're
### not using  gcc

defines = -D$(SYS)			# announce system by a '#define'

### where to find .h files:
incdirs = 
libdirs = 
libs = 
optimflags  = -c -g -Wall $(defines) $(incdirs) -O2 -ffast-math -finline-functions -funroll-all-loops -DNOSAFE -DNDEBUG -UDEBUG 
normalflags = -c -g -Wall $(defines) $(incdirs) 
CFLAGS = $(normalflags)

utils_h = utils.h list.h alloc.h hash.h vecmat.h pdb.h pdbutils.h extra-math.h \
          getargs.h jacobi.h files.h
utils_c = getargs.c utils.c alloc.c list.c extra-math.c jacobi.c \
	hash.c compat.c pdb.c pdbutils.c files.c
utils_o =  $(utils_c:.c=.o)

.PHONY: all utils clean debugmake help 

execs = conv usleep Pgamma phyper

all:  $(execs) 

install: $(execs)
	cp $(execs) $(bindir)

conv: conv.c
	$(CC) -g -Wall -O4 $< -o $@
Pgamma: Pgamma.c
	$(CC) $< -o $@ utils.o alloc.o -lm
usleep: usleep.c
	$(CC) -g -Wall -O4 $< -o $@

phyper: phyper.c
	$(CC) -g -Wall -O4 $< -lRmath -lm -o phyper
# (this requires Debian package r-mathlib)

# all/most utils depend on each other, so put all of utils_h  in  there:
$(utils_o): %.o : %.c $(utils_h)
	$(CC) $(CFLAGS) $*.c
	$(AR) sr libutils.a $*.o

help: 
	@echo "make all | utils | clean | help | debugmake"


clean:
	$(RM) *.o *.a conv usleep Pgamma

utils: $(utils_o)

$(utils_o):

# $(utils_o): $(utils_h)
# 
# $(utils_o):
# 	  $(CC) -c $^
# 	  $(AR) sr libutils.a $@

debugmake:
	pwd
	echo $(SYS) $(MAKELEVEL)
