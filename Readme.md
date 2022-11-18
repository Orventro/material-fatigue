# Material fatigue calculation

## How to make mesh

```
python3 ./util/gen_mesh.py 56 30 > ./meshes/mesh_56_30_circle_dmg.mesh
```

## How to run simulation

You might need to edit `CMakeLists.txt` to include directory of mfem lib.

```
cd build
cmake ..
make plate_load
./plate_load ../meshes/mesh_56_30_circle_dmg.mesh ../sims/sim1
```

The output will be saved in `./sims/sim1`. Simulation will take about 50k steps to converge.

## How to debug

```
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

then

```
./plate_load ../meshes/mesh_56_30_circle_dmg.mesh ../sims/sim1
```

or


```
valgrind ../meshes/mesh_basic_16_5.mesh ../sims/sim1
```

In line above valgrind is used and also the mesh is much smaller (since valgrind slows program down by a lot)
