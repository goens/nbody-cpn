#!/bin/bash

declare -a mapping_algos=("gbm-100" "gbm-10" "load-balancer")

for algo in "${mapping_algos[@]}" 
do
  for i in `seq 1 6`;
  do
    echo "Generate binary with $i worker and mapping obtained by $algo."
    sed -i "s/^__PNgeneric(NUM_WORKERS,[0-9]\+/__PNgeneric(NUM_WORKERS,$i/" nbody.cpnm4 
    grep "_PNgeneric(NUM_WORKERS" nbody.cpnm4
    make cpnm4
    make exynos_clean
    MAPPING=../mapper/mappings/$algo/$i/default.mapping make exynos
    cd ../generator/exynos
    ./patch_pthreads.sh nbody.exynos.c 8000
    CC_EXYNOS_USER_FLAGS=-O3 make 
    mkdir -p ../../mapper/mappings/bins
    mv nbody.exynos ../../mapper/mappings/bins/nbody.exynos.$algo.$i.O3
    cd ../../spec
  done
done 

