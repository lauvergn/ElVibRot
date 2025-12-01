#!/bin/bash

here=`pwd`
DIR_vib=$here/../..
DIR_pot=$DIR_vib/sub_pot

if [ -e $DIR_vib/EVR-T.exe ] 
then
  echo " No compilation EVR-T.exe exist"
  echo " EVR OK" #for the test in calc and calc_micro
  exit
fi

cd $DIR_vib
 make OPT=0
cd $here
compOK=`grep -c 'EVR OK' comp.log `
if [ $compOK = '1' ]
then
 echo " ElVibRot compilation: OK"
else
 echo " ElVibRot compilation: ERROR"
 exit 1
fi