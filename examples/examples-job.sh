#!/bin/sh

# script to produce the output in the example directory

tmpdir=t
mkdir $tmpdir 2>/dev/null
cd $tmpdir
echo "will do work in `pwd` to avoid overwriting the distribution's examples"

quilt_home=`pwd`/../..
quilt_scripts=$quilt_home/scripts
PATH=$quilt_scripts:$PATH
export PDBPATH=..:.

one-structure.sh 8lyz






