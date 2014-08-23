#!/bin/bash

icc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I "/usr/local/apps/matlab2014a/extern/include" svd_mex_client.c socket_functions.c -L"/usr/local/apps/matlab2014a/bin/glnxa64" -Wl,--version-script,"/usr/local/apps/matlab2014a/extern/lib/glnxa64/mexFunction.map"  -lmex -lmx -lmat -lm -o svd_mex_client.mexa64

icc -mkl -openmp svd_with_socket_server_intel_mkl.c socket_functions.c matrix_vector_functions_intel_mkl.c  -o svd_with_socket_server_intel_mkl


icc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I "/usr/local/apps/matlab2014a/extern/include" qr_full_rank_mex_client.c socket_functions.c -L"/usr/local/apps/matlab2014a/bin/glnxa64" -Wl,--version-script,"/usr/local/apps/matlab2014a/extern/lib/glnxa64/mexFunction.map"  -lmex -lmx -lmat -lm -o qr_full_rank_mex_client.mexa64

icc -mkl -openmp qr_full_rank_with_socket_server_intel_mkl.c socket_functions.c matrix_vector_functions_intel_mkl.c  -o qr_full_rank_with_socket_server_intel_mkl
