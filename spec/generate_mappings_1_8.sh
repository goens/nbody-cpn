#!/bin/bash

for i in `seq 1 8`;
do
  echo "Starting mapping generation with $i worker."
  sed -i "s/^__PNgeneric(NUM_WORKERS,[0-9]\+/__PNgeneric(NUM_WORKERS,$i/" nbody.cpnm4 
  grep "_PNgeneric(NUM_WORKERS" nbody.cpnm4
  make cpnm4
  make mapper_clean
  make mapper
  cd ../mapper
  ./patch_mapper.sh
  make mapper_bb
  mkdir -p ../mapper/mappings/$i
  mv `find ../mapper/results/ -maxdepth 1 -mindepth 1 -type d` ../mapper/mappings/$i
  cd ../spec
done
