#!/bin/bash

name_dep=dependencies.mk
SRCFile=fortranlist.mk

list=`ls Source_ElVibRot/*.f90 Source_ElVibRot/*/*.f90 Source_ElVibRot/*/*/*.f90`
ExcludeList='Gauss_numlH.f90 spectre.f90 sub_rotation.f90 sub_CRP-jo.f90 temp.f90'

echo "#===============================================" > $name_dep
echo "#===============================================" > $SRCFile
echo "SRCFILES= \\" >> $SRCFile

for ff90 in $list
do
   ff=`awk '{name=$1 ; n=split(name,tab,"/") ; if (n > 0) {l=length(tab[n]) ; print tab[n]}}' <<< $ff90`
   #echo $ff
   if grep -vq $ff <<< $ExcludeList;  then
     echo $ff " \\" >> $SRCFile
     awk -f scripts/mod2file.awk $ff90 >> $name_dep
   fi

   #awk -f scripts/mod2file.awk $ff90 >> $name_dep
   #echo $ff90 | awk '{name=$1
   #n=split(name,tab,"/")
  #if (n > 0) {
  #  l=length(tab[n])
  #  print tab[n] " \\"
  #}
  #}' >> $SRCFile
done
echo "#===============================================" >> $name_dep
for ff90 in $list
do
   ff=`awk '{name=$1 ; n=split(name,tab,"/") ; if (n > 0) {l=length(tab[n]) ; print tab[n]}}' <<< $ff90`
   if grep -vq $ff <<< $ExcludeList;  then
     awk -f scripts/dep2.awk $ff90 >> $name_dep
   fi
   #awk -f scripts/dep2.awk $ff90 >> $name_dep
done