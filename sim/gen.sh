#!/bin/bash

nx=256
ny=256

I0=1

zd_muln="0.02"
zd_sigln="2"

lam_muln="0.005"
lam_sigln="0.5"

apl_muln="0.005"
apl_sigln="0.5"

base_nlog="-8"
base_zd="0.00001"
base_zlog=$(echo "l($base_zd)/l(10)" | bc -l)
base_zdiff=$(echo "$base_nlog - $base_zlog" | bc -l)

cutoff_time=6000 # maximum time given for a pattern to be generated (in seconds)

num=10000 # number of patterns to generate
digs=$(echo "val=l($num - 1)/l(10);scale=0;val/1 + 1;" | bc -l) # number of digits in file name
counter=1

gen_pat() {
  local pref=$(printf "%0${digs}d" $counter)
  timeout $cutoff_time dffrcc --nx $nx --ny $ny --lam $1 --xa $2 --xb $3 --ya $2 --yb $3 --z $4 --w $5 --l $5 --I0 $I0 \
  --dffr $pref.dffr --bmp $pref.bmp --ptime 0 -v --cmap bgr > log$pref 2>&1
  return $?
}

while [ $counter -le $num ]; do
  lam=$(logn $lam_muln $lam_sigln)
  apl=$(logn $apl_muln $apl_sigln)
  if [ $(echo "$lam > $apl" | bc -l) == 1 ] || [ $(echo "$lam < $apl/20" | bc -l) == 1 ]; then
    continue;
  fi
  zd=$(logn $zd_muln $zd_sigln)
  wl=$(echo "1.5*$zd" | bc -l) # physical width and length of detector (since the detector is square)
  if [ $(echo "$wl <= 2*$apl" | bc -l) == 1 ]; then
    continue;
  fi
  xyb=$(echo "$apl/2" | bc -l) # upper x- and y-limits of the aperture (since the aperture is square)
  gen_pat $lam "-$xyb" $xyb $zd $wl
  if [ $? != 0 ]; then
    continue;
  fi
  echo -en "Pattern $counter generated with:\n\tz-dist: $zd m\n\tap_len: $apl m\n\tlambda: $lam m\n\n"
  counter=$((counter + 1))
done
