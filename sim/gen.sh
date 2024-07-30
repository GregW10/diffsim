#!/bin/bash

nx=256 # detector x-resolution
ny=256 # detector y-resolution

I0=10000 # incident light intensity in W/m^2

ptol="1e-15" # periodic tolerance

zd_muln="0.02" # mean distance to background screen
zd_sigln="4" # standard deviation in distance to background screen

lam_muln="0.005" # mean wavelength of light
lam_sigln="3" # standard deviation in wavelength

apl_muln="0.005" # mean aperture side-length
apl_sigln="3" # standard deviation in aperture side-length

l_max="-8"
l_b="-5"
z_b=$(echo "e($l_b*l(10))" | bc -l)
min_tol=$(echo "e($l_max*l(10))" | bc -l)

tol_cut="0.000000000001" # cut-off tolerance below which the simulation switches to using `long double` values
# tol_cutlog=$(echo "l($tol_cut)/l(10)" | bc -l)

lams="-8.5"
# lamb=$(echo "$lams + l($2)/l(10) - $l_b" | bc -l)

gen_lt() {
  logz=$(echo "l($1)/l(10)" | bc -l)
  lamdiff=$(echo "$lams + l($2)/l(10) - $logz" | bc -l)
  if [ $(echo "$1 <= $z_b" | bc -l) == 1 ] && [ $(echo "$lamdiff >= $l_max" | bc -l) == 1 ]; then
    echo $min_tol
    return 0
  fi
  zdiff=$(echo "$l_max + 1.5*($l_b - $logz)" | bc -l)
  if [ $(echo "$zdiff <= $lamdiff" | bc -l) == 1 ]; then
    echo "e($zdiff*l(10))" | bc -l
  else
    echo "e($lamdiff*l(10))" | bc -l
  fi
  return 0
}

base_nlog="-8"
base_zd="0.00001"
base_zlog=$(echo "l($base_zd)/l(10)" | bc -l)
base_zdiff=$(echo "$base_nlog - $base_zlog" | bc -l)

cutoff_time=300 # maximum time given for a pattern to be generated (in seconds)
cutoff_sleep=10
min_prog=$(echo "($cutoff_sleep/$cutoff_time)*100.0" | bc -l)

num_pats=100000 # number of patterns to generate
digs=$(echo "val=l($num_pats - 1)/l(10);scale=0;val/1 + 1;" | bc -l) # number of digits in file name
counter=1

gen_pat() {
  local pref=$(printf "%0${digs}d" $counter)
  if [ $(echo "$6 < $tol_cut" | bc -l) == 1 ]; then
    ftype="Lf"
  else
    ftype="lf"
  fi
  timeout $cutoff_time dffrcc --nx $nx --ny $ny --lam $1 --xa $2 --xb $3 --ya $2 --yb $3 --z $4 --w $5 --l $5 --I0 $I0 \
  --ptol_x $ptol --ptol_y $ptol --atol_x $6 --atol_y $6 --rtol_x $6 --rtol_y $6 -f $ftype \
  --dffr $pref.dffr --bmp $pref.bmp --ptime 1 -v --cmap grayscale > log$pref 2>&1 &
  pid=$!
  sleep $cutoff_sleep
  if ! [ -z "$(grep "100.00% complete." log$pref)" ]; then
    return 0
  fi
  num=$(tail -1 log$pref | awk -F '%' '{ print $1 }')
  # echo "num: $num"
  if [ $(echo "$num >= $min_prog" | bc -l) == 1 ]; then
    wait $pid
    return $?
  fi
  kill $pid
  echo "The simulation with the following parameters would have taken approx. $(printf "%.2Lf" $(echo "$cutoff_sleep*(100/$num)" | bc -l)) s to complete:"
  cat log$pref && echo
  return 1
}

while [ $counter -le $num_pats ]; do
  lam=$(logn $lam_muln $lam_sigln)
  apl=$(logn $apl_muln $apl_sigln)
  if [ $(echo "$lam > $apl" | bc -l) == 1 ] || [ $(echo "$lam < $apl/20" | bc -l) == 1 ]; then
    continue;
  fi
  zd=$(logn $zd_muln $zd_sigln)
  wl=$(echo "1.5*$zd" | bc -l) # physical width and length of detector (since the detector is square)
  one_p5_apl=$(echo "1.5*$apl" | bc -l)
  if [ $(echo "$wl < $one_p5_apl" | bc -l) == 1 ]; then
    wl=$one_p5_apl # have added this cut-off (rather than just disallowing it) to permit very near-field patterns
  fi
  xyb=$(echo "$apl/2" | bc -l) # upper x- and y-limits of the aperture (since the aperture is square)
  tol=$(gen_lt $zd $lam)
  if ! gen_pat $lam "-$xyb" $xyb $zd $wl $tol; then
    continue;
  fi
  echo -en "Pattern $counter generated with:\n\tz-dist: $zd m\n\tap_len: $apl m\n\tlambda: $lam m\n\t   tol: 0$tol\n\n"
  counter=$((counter + 1))
done
