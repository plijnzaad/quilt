#! /bin/sh -x
## Written by plijnzaad@gmail.com This script produces all the stuff
## necessary for analyzing one particular protein structure. The pdb filename
## must be given as a parameter.

if [ $# -ne 1 ]; then
    echo "Usage: one-structure.sh PDBCODE"
    exit 21
fi

### Note: quilt and the scripts are assumed to be on the $PATH. If this is
### not the case, uncomment the following:

# quilt_home=$HOME/quilt                  # change to needs
# # quilt_path=$quilt_home/scripts:$quilt_home # if binary not installed in $PATH
# quilt_path=$quilt_home/scripts          # just for scripts
# PATH=$quilt_path:$PATH

pdbcode=$1
pdbfile=$pdbcode.pdb
## It is searched for in $PDBPATH.  Provide sensible default
if [ ${PDBPATH-notdefined} = notdefined ]; then
    export PDBPATH=.:..:$HOME
fi

areafile=$pdbcode.area
patfile=$pdbcode.pat
logfile=$pdbcode.log
histofile=$pdbcode.histo
ranhistofile=$pdbcode.ranhisto
tmpdir=./ran                            # where all the randomized pat files go
rm -fr $tmpdir                          # has to be clean first
mkdir $tmpdir 
textfile=$pdbcode.txt

polar_expansion=1.4

# run it on structure
quilt=quilt # must be in path, e.g. ~/Linux/bin/quilt
$quilt -ep $polar_expansion -R -a $areafile  -p $pdbfile \
    > $patfile 2> $logfile

# parse it into something more readable
patresidues < $patfile > $textfile

# produce rasmol script for the 5 largest patches, called p0, p1 .. p4:
for i in 0 1 2 3 4 
do 
  Rpatsel $patfile $i  
done > $pdbcode.rasmol

# produce histograms and randomized histograms to assess significance:
awk '/^# [0-9]/{print $10}' $patfile | histo  -l 0 -s 10 > $histofile

# now randomize the surface n times
ntimes=100

echo "now randomizing $ntimes times ..."
for i in $(seq 1 $ntimes); do 
  $quilt -ran  -n $npoints -ep $polar_expansion -R -p $areafile \
     > $tmpdir/$pdbcode.ran$i.pat 
done 2>> $logfile

# see if none failed (this still happens occasionally ...)
find $tmpdir -size 0 -exec echo 'randomization failed' \;  -exec  rm -v {} \; 

#adjust ntimes so histo's end up at correct height:
ntimes=$(ls ran/*.ran*|wc -l)

cat $tmpdir/$pdbcode.ran* | awk '/^# [0-9]/{print $10}' | histo -l 0 -s 10 \
  | awk '{print $1, $2/'$ntimes'}' > $ranhistofile

set +x

echo "Done. Patches are in $pdbcode.pat (more readeable version: $pdbcode.txt); 
Rasmol script is in $pdbcode.rasmol. Logs are in $pdbcode.log.

The patch size distribution of $pdbcode is in $pdbcode.histo
The patch size distribution of $ntimes surface-randomized versions of $pdbcode 
is in $pdbcode.ranhisto. This may help decide which patches are 'really large'.
 
 "
