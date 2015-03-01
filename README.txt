Socket Code for Mex Interface
Includes example of calling from Matlab the SVD from NVIDIA CULA and 
SVD and pivoted QR from Intel MKL library
Written by Sergey Voronin, 2014 - 2015. 
Should work on any recent Linux distribution.
Tested on Ubuntu 14.04.2 with cuda 6, cula R16, icc 14.03 and Matlab 2013b and 2014a.

========= what is this ==============================

This is a socket interface for passing data between Matlab and an external program. 
This is useful for performing a computation in Matlab which is programmed 
with an external library (probably with parallelization - such as on GPU). 
In the examples here, we create 
a code with the NVIDIA CULA library and the Intel MKL library to do the SVD 
and pivoted QR decompositions and are able 
to call these functions from Matlab using the socket interface. This is just an example, 
the socket interface can be used for other libraries too. 
This is useful because it is often difficult to build Matlab 
mex files linked with external libraries (like CULA and MKL) and running 
in parallel. This code 
completely separates the mex file and the compute engine (the compute engine can even 
run on a different machine, so for example you can run a Matlab code on a laptop and 
have the SVD computation inside be done for you on a local server).


========= how it works ==============================

Matlab passes data (matrix) to the wrapper. Wrapper m-file calls up the socket server 
with a randomly chosen port number between 1000 and 2000 (note: if you have a firewall, 
please adjust port numbers accordingly in the wrapper). Then the wrapper 
passes the data to the mex client and the mex client communicates with the socket server 
via sockets. The server does the computation and sends the results back to the mex client 
via sockets. The mex client receives results from the compute server and passes the results 
to the Matlab wrapper which returns them to the Matlab session or script from which it was called.
The reads and writes between the socket server and client are performed in a synchronous 
manner using a small amount of disk i/o (i.e. server sends client data when client is ready 
and vice versa). This is done to minimize possibility of errors with transfers.

========= summary of installation and usage =========

The following steps are for the use of the included SVD examples. The examples use 
the Intel MKL library and the NVIDIA CULA library and are callable 
from Matlab. The socket code can be used with different libraries also.

1) Install Intel C Compiler (icc) and Math Kernel Library (MKL) 
OR gcc and the CULA library

For Intel MKL:
Follow instructions from Intel website (basically download install package and unpack)
Afterwards, to set paths, use (in bash shell):

source /opt/intel/bin/compilervars.sh intel64

For NVIDIA CULA:
install CUDA 
get and install CULA package from http://www.culatools.com/downloads/dense/ 
change parameters accordingly in setup_vars.sh and then run:
source setup_vars.sh 

2) Compile code (svd and qr codes from mkl library, mex files, socket client/server) - use the included compile.sh script

The compilation looks something like below for the SVD example, similar for QR. 
You must put the proper path to your Matlab installation.

mkdir logs/

MATLAB_INCLUDE="/usr/local/MATLAB/R2013b/extern/include"
MATLAB_LIB="/usr/local/MATLAB/R2013b/bin/glnxa64"
MATLAB_MAP_FILE="/usr/local/MATLAB/R2013b/extern/lib/glnxa64/mexFunction.map"

SOCKET_INCLUDE="../socket_functions"

icc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I $MATLAB_INCLUDE -I $SOCKET_INCLUDE svd_mex_client.c "$SOCKET_INCLUDE"/socket_functions.c -L"$MATLAB_LIB" -Wl,--version-script,"$MATLAB_MAP_FILE"  -lmex -lmx -lmat -lm -o svd_mex_client.mexa64

icc -mkl -openmp -I "$SOCKET_INCLUDE" svd_with_socket_server_intel_mkl.c "$SOCKET_INCLUDE"/socket_functions.c matrix_vector_functions_intel_mkl.c  -o svd_with_socket_server_intel_mkl

For the CULA code it looks something like this:

gcc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I $MATLAB_INCLUDE -I $SOCKET_INCLUDE svd_mex_client.c "$SOCKET_INCLUDE"/socket_functions.c -L"$MATLAB_LIB" -Wl,--version-script,"$MATLAB_MAP_FILE"  -lmex -lmx -lmat -lm -o svd_mex_client.mexa64

nvcc -I "$SOCKET_INCLUDE" svd_with_socket_server_nvidia_cula.c "$SOCKET_INCLUDE"/socket_functions.c matrix_vector_functions_nvidia_cula.c  -o svd_with_socket_server_nvidia_cula -Xcompiler -fopenmp -lcula_lapack -lcublas -lcudart -liomp5 -L$CULA_LIB_PATH_64 -I$CULA_INC_PATH

You can modify the compile.sh script with proper paths and run that to compile all the examples.

3) Run the code

First make a directory called logs (this was done in the compilation 
step above). This is needed for i/o communication between 
server and client.

mkdir logs/

Open Matlab in same directory as the compiled code.
make sure intel mkl paths are set properly in the terminal shell 
where Matlab is launched from
i.e. source /opt/intel/bin/compilervars.sh intel64
If using CULA, make sure you source the CULA vars.

Create matrix:

>> A = rand(2000,3000);

Run SVD CULA wrapper:

>> [U,S,Vt] = svd_wrapper_nvidia_cula(A); % computation done on GPU

Check result:

>> norm(A - U*S*Vt)

ans =

   2.6487e-12

>> 

Run SVD MKL wrapper:

>> [U,S,Vt] = svd_wrapper_intel_mkl(A); % computation done on CPU with MKL lib

Check result:

>> norm(A - U*S*Vt)                    

ans =

   1.1364e-11

>>  norm(U'*U-eye(size(U,2)))

ans =

   2.9218e-14

>> norm(Vt*Vt'-eye(size(Vt,1)))

ans =

   4.6965e-14

>> norm(S-diag(diag(S)))

ans =

     0

>> 

For the QR wrapper:

>> [Q,R,P,I] = qr_full_rank_wrapper_intel_mkl(A);

>> norm(A(:,I) - Q*R)

ans =

   1.0397e-12

>> norm(A - Q*R*P')     

ans =

   1.0397e-12

>> norm(Q*Q' - eye(size(Q,2)))

ans =

   7.3662e-15

>> 

