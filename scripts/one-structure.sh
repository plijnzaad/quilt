#! /bin/sh -x
## Written by Philip lijnzaad@ebi.ac.uk
## This script produces all the stuff necessary for analyzing one particular protein
## structure. The pdb filename must be given as a parameter.



quilt_home=$HOME/quilt                  # change to needs
# quilt_path=$quilt_home/scripts:$quilt_home # if binary not installed in $PATH
quilt_path=$quilt_home/scripts          # just for scripts

PATH=$quilt_path:$PATH

pdbcode=$1                              # assumes something.pdb
pdbfile=$pdbcode.pdb
## It is searched for in PDBPATH, so set that to at least include "."
export PDBPATH=.

areafile=$pdbcode.area
patfile=$pdbcode.pat
logfile=$pdbcode.log
histofile=$pdbcode.histo
ranhistofile=$pdbcode.ranhisto
tmpdir=./tmp                            # where all the randomized pat files go
textfile=$pdbcode.txt

polar_expansion=1.4
npoints=252

# run it on structure
quilt  -n $npoints -ep $polar_expansion -R -a $areafile $pdbfile \
    > $patfile 2> $logfile

# parse it into something more readable
patresidues < $patfile > $textfile

# produce rasmol script for the 5 largest patches
for i in 0 1 2 3 4 
do 
  Rpatsel $patfile $i  
done > $pdbcode.rasmol


# produce histograms and randomized histograms to assess significance:
awk '/^# [0-9]/{print $10}' $patfile | histo -l 0 -s 10 > $histofile

for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 
do 
  quilt -ran -n 252 -ep 1.4 -R -p $areafile > $tmpdir/$pdbcode.ran$i 2> /dev/null
done

cat $tmpdir/$pdbcode.ran[01][0-9] | awk '/^# [0-9]/{print $10}' | histo -l 0 -s 10 \
  | awk '{print $1, $2/20.0}' > $ranhistofile

