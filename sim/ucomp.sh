#!/bin/bash

cp diffim_rct.cpp diffim_rct.cu
nvcc -std=c++20 -O3 diffim_rct.cu -o ~/coms/dffrc
ret=$?
rm diffim_rct.cu
exit $ret
