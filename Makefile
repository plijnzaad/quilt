# -*- mode: Makefile; -*-

# This file requires GNU make


include $(BUILD_HOME)/makefile.include

.PHONY:  all quilt utils install clean debugmake

all: quilt utils

quilt:
	(cd src; $(MAKE))

utils:
	(cd utils; $(MAKE) utils) || true;	#may have been installed elsewhere

test:
	(env PDBPATH=examples ./src/quilt -n 252 -ep 1.4 -R  -p 8lyz.pdb > 8lyz-test.pat 2> 8lyz-test.log; if diff -q 8lyz-test.pat examples/8lyz.pat ; then echo "test OK"; exit 0; else echo "test not OK: 8lyz-test.pat and examples/8lyz.pat differ"; exit 9; fi)

install:
	(cd src; $(MAKE) install)
	$(CP) -p scripts/* $(bindir)

clean:
	(cd src; $(MAKE) -k clean)
	(cd utils; $(MAKE) -k clean)
	(if [ -d test ] ; then cd test;  $(MAKE) -k clean; fi)
	rm $(bindir)/{patresidues,histo,one-structure.sh,Rpatsel}

debugmake:
	pwd
	echo SYS=$(SYS) $(MAKELEVEL)
	(cd src; $(MAKE) debugmake)
	(cd utils; $(MAKE) debugmake)

