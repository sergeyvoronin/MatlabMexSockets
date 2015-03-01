#!/bin/bash

mkdir logs/

MATLAB_INCLUDE="/usr/local/MATLAB/R2013b/extern/include"
MATLAB_LIB="/usr/local/MATLAB/R2013b/bin/glnxa64"
MATLAB_MAP_FILE="/usr/local/MATLAB/R2013b/extern/lib/glnxa64/mexFunction.map"

SOCKET_INCLUDE="../socket_functions"

icc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I $MATLAB_INCLUDE -I $SOCKET_INCLUDE svd_mex_client.c "$SOCKET_INCLUDE"/socket_functions.c -L"$MATLAB_LIB" -Wl,--version-script,"$MATLAB_MAP_FILE"  -lmex -lmx -lmat -lm -o svd_mex_client.mexa64

icc -mkl -openmp -I "$SOCKET_INCLUDE" svd_with_socket_server_intel_mkl.c "$SOCKET_INCLUDE"/socket_functions.c matrix_vector_functions_intel_mkl.c  -o svd_with_socket_server_intel_mkl


icc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I $MATLAB_INCLUDE -I $SOCKET_INCLUDE qr_full_rank_mex_client.c "$SOCKET_INCLUDE"/socket_functions.c -L"$MATLAB_LIB" -Wl,--version-script,"$MATLAB_MAP_FILE"  -lmex -lmx -lmat -lm -o qr_full_rank_mex_client.mexa64

icc -mkl -openmp -I "$SOCKET_INCLUDE" qr_full_rank_with_socket_server_intel_mkl.c "$SOCKET_INCLUDE"/socket_functions.c matrix_vector_functions_intel_mkl.c  -o qr_full_rank_with_socket_server_intel_mkl


