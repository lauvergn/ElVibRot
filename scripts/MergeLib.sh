#!/bin/bash

nArg=$#
nLib=$((nArg-1))

echo Script name: $0
echo The script has $# arguments 
echo Number of libraries: $nLib

if test "$nLib" -le 1 ; then 
  echo "The script has only "$nLib" library."
  echo "It must have 2 or more libraries."
  exit 1
fi

NewLib=$1
rm -f $NewLib
rm -rf tempOBJ

mkdir tempOBJ || exit 1
cd tempOBJ
pwd

for Lib in "$@"
do
  test "$Lib" = "$NewLib" && continue
  echo $Lib "to be merged"

  if test ! -e "../$Lib" ; then
    echo The library, $Lib, does not exist! 
    exit 1
  fi
  ar -x ../$Lib
  ar -qc ../$NewLib *.o
  rm *.o
done

cd ..
ar -s $NewLib
rm -rf tempOBJ
