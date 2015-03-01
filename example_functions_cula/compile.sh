#!/bin/bash

mkdir logs/

MATLAB_INCLUDE="/usr/local/MATLAB/R2013b/extern/include"
MATLAB_LIB="/usr/local/MATLAB/R2013b/bin/glnxa64"
MATLAB_MAP_FILE="/usr/local/MATLAB/R2013b/extern/lib/glnxa64/mexFunction.map"

SOCKET_INCLUDE="../socket_functions"

gcc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I $MATLAB_INCLUDE -I $SOCKET_INCLUDE svd_mex_client.c "$SOCKET_INCLUDE"/socket_functions.c -L"$MATLAB_LIB" -Wl,--version-script,"$MATLAB_MAP_FILE"  -lmex -lmx -lmat -lm -o svd_mex_client.mexa64

nvcc -I "$SOCKET_INCLUDE" svd_with_socket_server_nvidia_cula.c "$SOCKET_INCLUDE"/socket_functions.c matrix_vector_functions_nvidia_cula.c  -o svd_with_socket_server_nvidia_cula -Xcompiler -fopenmp -lcula_lapack -lcublas -lcudart -liomp5 -L$CULA_LIB_PATH_64 -I$CULA_INC_PATH

