#!/bin/bash

# patch .ptt
cd .ptt
make clean 
../../generator/pthreads/patch_pthreads.sh nbody.Pthreads.c 8000
make
cd ..

# patch .cpntrace
cd .cpntrace
CFILES=`find . -name '*.c'`
for i in $CFILES; do ../../generator/pthreads/patch_pthreads.sh $i 8000; done
