# -*- mode: Makefile; -*-
# $Id$

# This file requires GNU make


include $(BUILD_HOME)/makefile.include

.PHONY:  all quilt utils install clean debugmake

all: quilt utils

quilt:
	(cd src; $(MAKE))

utils:
	(cd utils; $(MAKE)) || true;	#may have been installed elsewhere

# test:
# 	(cd test; $(MAKE))

install:
	(cd src; $(MAKE) install)

clean:
	(cd src; $(MAKE) clean)
	(cd utils; $(MAKE) clean)
	(if [ -d test ] ; then cd test;  $(MAKE) clean; fi)

debugmake:
	pwd
	echo SYS=$(SYS) $(MAKELEVEL)
	(cd src; $(MAKE) debugmake)
	(cd utils; $(MAKE) debugmake)

