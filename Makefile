cpu:
	g++ -O3 plate_loading.cpp -I ../mfem-4.4/build -L ../mfem-4.4/build -lmfem -o plate_load

par:
	mpicxx -O3 -std=c++11 plate_loading.cpp -I ../mfem-4.4/build_parallel -I ../hypre/src/hypre -I ../hypre/src/seq_mv -I ../hypre/src/hypre/include \
		-L ../hypre/src/hypre -lHYPRE \
		-L ../metis -lmetis \
		-L ../mfem-4.4/build_parallel -lmfem \
		-o plate_load_par

cuda:
	nvcc -O3 plate_loading.cpp -g -I ../mfem-4.4/build_cuda -L ../mfem-4.4/build_cuda -lmfem -lcusparse -o plate_load_cuda