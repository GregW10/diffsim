#!/bin/bash

nx=256
ny=256

num=10000
digs=$(echo "val=l($num - 1)/l(10);scale=0;val/1 + 1;" | bc -l)
counter=0

gen_pat() {
  local pref=$(printf "%0${digs}d")
  dffrcg --nx $nx --ny $ny --lam $1 --xa $2 --xb $3 --ya $4 --yb $5 --z $6 --w $7 --l $8 --I0 $9 --dffr $pref.dffr \
  --bmp $pref.bmp
}


