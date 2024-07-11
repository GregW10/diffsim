#!/bin/bash

nx=256
ny=256

I0=1

apl=()

zdist=("0.00001" "0.00002" "0.00005" "0.00010" "0.00020" "0.00050" "0.00100" "0.00200" "0.00500" "0.01000" "0.02000" \
"0.05000" "0.10000" "0.20000" "0.50000" "1.00000") # distance to detector

lam_frac=20

echo -e "Starting wavelength is $(echo "scale=2;$lam_i*1000000000" | bc -l) nm.\nFinal wavelength is \
$(echo "scale=2;$lam_f*1000000000" | bc -l) nm.\nStep is $(echo "scale=2;$lam_d*1000000000" | bc -l) nm."

num=10000 # number of patterns to generate
digs=$(echo "val=l($num - 1)/l(10);scale=0;val/1 + 1;" | bc -l) # number of digits in file name
counter=0

gen_pat() {
  local pref=$(printf "%0${digs}d" $counter)
  dffrcg --nx $nx --ny $ny --lam $1 --xa $2 --xb $3 --ya $2 --yb $3 --z $4 --w $5 --l $5 --I0 $I0 --dffr $pref.dffr \
  --bmp $pref.bmp
}

for zd in "${zdist[@]}"; do
  apl+=( $zd )
  for aplen in "${apl[@]}"; do
    lam_s=$(echo "$aplen/$lam_frac" | bc -l)
    lam_val=$lam_s
    upper_lim=$(echo "$aplen + $aplen/($lam_frac*2)" | bc -l) # to account for f.p. errors
    while [ $(echo "$lam_val <= $upper_lim" | bc -l) == 1 ]; do
      xyb=$(echo "$aplen/2" | bc -l) # upper x- and y-limits of the aperture (since the aperture is square)
      wl=$(echo "2*$zd" | bc -l) # physical width and length of detector (since the detector is square)
      # gen_pat $lamv "-$xyb" $xyb $zd $wl
      # echo $lamv $xyb $zd
      echo -en "Pattern $((counter + 1)) generated with:\n\tz-dist: $zd m\n\tap_len: $aplen m\n\tlambda: $lam_val m\n\n"
      lam_val=$(echo "$lam_val + $lam_s" | bc -l)
      counter=$((counter + 1))
    done
  done
done
