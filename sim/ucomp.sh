#!/bin/bash

cp diffim.cpp diffim.cu
nvcc -std=c++20 -O3 diffim.cu -o ~/coms/diffimg
rm diffim.cu
