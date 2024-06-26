#!/bin/bash

cp diffim_rct.cpp diffim_rct.cu
nvcc -std=c++20 -O3 diffim_rct.cu -o ~/coms/dffrc
rm diffim_rct.cu
