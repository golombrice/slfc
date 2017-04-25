#!/bin/bash

gcc -c arithmetic_decode.c
gcc -c arithmetic_encode.c
gcc -c bit_output.c
gcc -c bit_input.c
gcc -c adaptive_model.c
gcc -c cerv.c
#gcc -o depth_compressor bit_input.o main.o adaptive_model.o bit_output.o arithmetic_encode.o arithmetic_decode.o -lm
ar rcs libcerv.a *.o
