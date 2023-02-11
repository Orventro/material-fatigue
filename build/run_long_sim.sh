#!/bin/bash

make plate_load
for i in $(seq 0 10)
do
    echo "./plate_load $i ../meshes/mesh_96_50_numbered.mesh ../sims/sim3"
    ./plate_load $i ../meshes/mesh_96_50_numbered.mesh ../sims/sim3
done
