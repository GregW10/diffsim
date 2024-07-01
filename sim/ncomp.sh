#!/bin/bash

cp diffim_rct.cpp temp.cu
nvcc -std=c++20 -O3 temp.cu -o ~/coms/dffrcg
ret=$?
rm temp.cu
exit $ret
