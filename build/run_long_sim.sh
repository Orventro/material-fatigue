#!/bin/bash


if [ -z "$2"]
  then
    mesh=../../meshes/mesh_100_35_plate_small.mesh
else
    mesh=$2
fi

if [ -z "$1" ]
  then
    echo "Arguments:"
    echo "    1. simulation directory number"
else
    mkdir -p ../sims/sim$1
    make plate_load || break
    mv plate_load ../sims/sim$1/

    cd ../sims/sim$1

    # for simNum in  $(seq 0 500)
    # do
    simNum=0
    echo "./plate_load $simNum $mesh ."
    ./plate_load $simNum $mesh . || break
    # done
fi
