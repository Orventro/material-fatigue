cmake_minimum_required(VERSION 3.18)
project(plate)

set(LIB_LOCATION "/home/orventro/Programms")
set(CMAKE_CXX_COMPILER g++)
add_compile_options(-O3)

add_library(mfem STATIC IMPORTED GLOBAL)
set_target_properties(mfem PROPERTIES IMPORTED_LOCATION ${LIB_LOCATION}/mfem-4.4/build/libmfem.a)
add_library(mfem_cuda STATIC IMPORTED GLOBAL)
set_target_properties(mfem_cuda PROPERTIES IMPORTED_LOCATION ${LIB_LOCATION}/mfem-4.4/build_cuda/libmfem.a)

add_executable(plate_load plate_loading.cpp ElasticWaveOperator.cpp)
target_include_directories(plate_load PRIVATE "../mfem-4.4/build")
target_link_libraries(plate_load mfem)

add_executable(plate_load_cuda plate_loading.cpp ElasticWaveOperator.cpp)
target_include_directories(plate_load_cuda PRIVATE "../mfem-4.4/build_cuda")
target_link_libraries(plate_load_cuda mfem_cuda)
target_compile_definitions(plate_load_cuda PRIVATE MFEM_CUDA)
