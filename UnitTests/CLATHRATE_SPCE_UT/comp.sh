#!/bin/bash

here=`pwd`
DIR_vib=$here/../..
DIR_pot=$DIR_vib/sub_pot

ls $DIR_pot
cp sub_system_sc_potSPCE-LF.f $DIR_pot/sub_system.f
cp calc_f2_f1Q.f90            $DIR_vib/sub_operator_T/

grep "Alavi" $DIR_pot/sub_system.f


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