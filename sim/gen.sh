#!/bin/bash

nx=256
ny=256

I0=1

lam_i="0.000000380" # 380 nm - approx. bottom of human perception
lam_f="0.000000740" # 740 nm - approx. top of human perception
lam_n=21 # number of wavelength values to use
lam_d=$(echo "($lam_f - $lam_i)/($lam_n - 1)" | bc -l) # step in wavelength

apsl=("0.00001" "0.00002" "0.00005" "0.00010" "0.00020" "0.00050" "0.00100" "0.00200" "0.00500" "0.01000" "0.02000" \
"0.05000") # aperture side-length

zdist=("0.00001" "0.00002" "0.00005" "0.00010" "0.00020" "0.00050" "0.00100" "0.00200" "0.00500" "0.01000" "0.02000" \
"0.05000" "0.10000" "0.20000" "0.50000" "1.00000") # distance to detector

echo -e "Starting wavelength is $(echo "scale=2;$lam_i*1000000000" | bc -l) nm.\nFinal wavelength is \
$(echo "scale=2;$lam_f*1000000000" | bc -l) nm.\nStep is $(echo "scale=2;$lam_d*1000000000" | bc -l) nm."

num=10000 # number of patterns to generate
digs=$(echo "val=l($num - 1)/l(10);scale=0;val/1 + 1;" | bc -l) # number of digits in file name
counter=0

gen_pat() {
  local pref=$(printf "%0${digs}d" $counter)
  dffrcg --nx $nx --ny $ny --lam $1 --xa $2 --xb $3 --ya $4 --yb $5 --z $6 --w $7 --l $8 --I0 $I0 --dffr $pref.dffr \
  --bmp $pref.bmp
}
