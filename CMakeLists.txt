cmake_minimum_required(VERSION 3.18)
project(material-fatigue)

set(LIB_LOCATION "/home/orventro/Programms")
set(CMAKE_CXX_COMPILER g++)
# add_compile_options(-O3)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

set(CMAKE_LDFKA)

set(CMAKE_CUDA_COMPILER nvcc)

add_library(mfem STATIC IMPORTED GLOBAL)
set_target_properties(mfem PROPERTIES IMPORTED_LOCATION ${LIB_LOCATION}/mfem-4.4/build/libmfem.a)
add_library(mfem_cuda STATIC IMPORTED GLOBAL)
set_target_properties(mfem_cuda PROPERTIES IMPORTED_LOCATION ${LIB_LOCATION}/mfem-4.4/build_cuda/libmfem.a)
add_library(mfem_mpi STATIC IMPORTED GLOBAL)
set_target_properties(mfem_mpi PROPERTIES IMPORTED_LOCATION ${LIB_LOCATION}/mfem-4.4/build_parallel/libmfem.a)

add_executable(plate_load plate_loading.cpp ElasticWaveOperator.cpp)# FatigueCalculator.cpp)
target_include_directories(plate_load PRIVATE "../mfem-4.4/build")
target_link_libraries(plate_load mfem)

add_executable(plate_load_cuda plate_loading.cpp ElasticWaveOperator.cpp)
target_include_directories(plate_load_cuda PRIVATE "../mfem-4.4/build_cuda")
target_include_directories(plate_load_cuda PRIVATE "/usr/local/cuda-11.4/include")
target_link_directories(plate_load_cuda PRIVATE "/usr/local/cuda-11.4/targets/x86_64-linux/lib/")
target_link_libraries(plate_load_cuda mfem_cuda cudart cusparse)
target_compile_definitions(plate_load_cuda PRIVATE MFEM_USE_CUDA)

add_executable(plate_load_mpi plate_loading.cpp ElasticWaveOperator.cpp)
target_link_libraries(plate_load_mpi mfem_mpi mpi mpi_cxx metis HYPRE)
target_include_directories(plate_load_mpi PRIVATE "/usr/include/suitesparse/")
target_link_directories(plate_load_mpi PRIVATE "../hypre/src/lib")
target_link_directories(plate_load_mpi PRIVATE "../metis-4.0/")
# target_include_directories(plate_load_mpi PRIVATE "/usr/local/cuda-11.8/targets/x86_64-linux/include/")
target_include_directories(plate_load_mpi PRIVATE "../mfem-4.4/build_parallel")
target_include_directories(plate_load_mpi PRIVATE "/usr/lib/x86_64-linux-gnu/openmpi/include/")
target_include_directories(plate_load_mpi PRIVATE "../hypre/src/hypre/")
target_include_directories(plate_load_mpi PRIVATE "../hypre/src/hypre/include/")