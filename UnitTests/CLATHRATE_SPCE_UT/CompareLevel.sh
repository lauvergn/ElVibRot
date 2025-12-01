#!/bin/bash

RES=$1

echo "TEST for: " $RES

if [ -e comp.log ]
then
  compOK=`grep -c 'EVR OK' comp.log `
else
  compOK=1
fi

if [ "$compOK" = '1' ]
then
 echo " ElVibRot compilation: OK"
else
 echo " ElVibRot compilation: ERROR"
 exit
fi

epsi=$2
if [ -z $epsi ]
then
  epsi=0.1
fi


EVROK=`grep -c "ElVibRot-Tnum AU REVOIR" $RES `
echo EVROK: $EVROK
if [ "$EVROK" -ne "1" ]
then
 echo "ElVibRot execution: ERROR"
 #exit
fi

LANG=C awk '{
if (NR == 1) {nbE=NF ; for (i=1;i<=nbE;i++) Eref[i]=$i}
else {for (i=1;i<=nbE;i++) DE[i]=sqrt((Eref[i]-$i)^2)}
}
END {
DEmax=0.
for (i=1;i<=nbE;i++) if (DE[i] > DEmax) DEmax=DE[i]
print "DEmax: " DEmax
}' lev_comp

