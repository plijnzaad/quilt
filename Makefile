# -*- mode: Makefile; -*-
# $Id$

# This file requires GNU make


include $(BUILD_HOME)/makefile.include

.PHONY:  all quilt utils install clean

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
	(cd utils; $(MAKE) clean) || true;
	(cd test; $(MAKE) clean)

debug:
	pwd
	echo $(SYS) $(MAKELEVEL)
	(cd src; $(MAKE) debug)
