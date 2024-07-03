#!/bin/bash

cp diffim_rct.cpp temp.cu
nvcc -std=c++20 -O3 -rdc=true temp.cu -o ~/coms/dffrcg
ret=$?
rm temp.cu
exit $ret
