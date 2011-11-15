#! /bin/bash
# $Id$

# Construct tar file for distribution. Should be invoked from within the
# quilt toplevel directory, and puts resulting tar-file  one level up.
# The version number is contained within the file VERSION.

# where to find stuff:
privcvs=":ext:login.genomicscenter.nl:/home/gen/philip/cvsroot"
# where is the histo script:
histo=/home/philip/perl/histo

echo "Need to fix templates.c AFREEA bug at some point first" >&2 
basedir=`pwd`

if [ $(basename $basedir) != 'quilt' ]; then
    echo "invoke from within quilt directory!"
    exit 3
fi

manifest=$basedir/MANIFEST
version=`cat VERSION`

# check if we can build:

if cvs -q up 2>&1 | egrep '^[?M]|nothing known'; then  
# if false; then                # debugging purposes
  echo '
***************************************************************
Current CVS working copy not completely up-to-date or committed. 
Not building tar file; please solve/commit first! (or specify env force_release=TRUE)
***************************************************************
' >&2
  if [ x$force_release != xTRUE ]; then
    exit 2
  fi
else 
  echo 'CVS ok' >&2
fi

cd ..

source_dir=quilt
release_dir=quilt-$version
tarfile=$release_dir.tgz

## get rid of CVS cruft:
rm -fr $release_dir
cvs -d $privcvs export -d $release_dir -D today quilt || exit 8

# missing stuff: just copy it from me:
cp $histo $release_dir/scripts/
mkdir $release_dir/utils
utilfiles=$(cat $manifest | sed -n  's/[ 	]*#.*//;/[a-z]/p;' | grep 'utils/')
for file in $utilfiles; do
     cp $source_dir/$file $release_dir/utils/
done

# get the examples in the right place:
cd $release_dir/examples
./examples-job.sh || exit 1
mv t/8lyz* t/ran/8lyz.ran* ./
cd ..
rm -fr examples/t
pwd

# process the templates:
cd $release_dir
sed "s/<<version>>/$version/" README.tmpl  > README
sed "s/<<version>>/$version/" INSTALL.tmpl  > INSTALL
echo Version $version built by `whoami`@`hostname` on `date` > BUILDDATE

rm README.tmpl INSTALL.tmpl

# Check if all files are there:
pwd
files=$(cat $manifest | sed -n  's/[ 	]*#.*//;/[a-z]/p;')
for file in $files; do
  if [ -f $file -o -L $file ]; then
    :
  else
    echo "$manifest says we need file '$file'; not found" >&2
    exit 3
  fi
done
cd ..

tar  --exclude-from=$source_dir/EXCLUDE -zhcvf $tarfile  $release_dir || exit 4
echo "done creating $tarfile"

# rm -fr $release_dir
