Socket Code for Mex Interface
Includes example of SVD and pivoted QR from Intel MKL library
Written by Sergey Voronin. 
Last updated: August 2014.
Note: this is the initial release (v. 0.1), need more error checking
Should work on any recent Linux distribution.
Tested on opensuse 13.1 with icc 14.03 and matlab 2013b and 2014a.

========= what is this ==============================

This is a socket interface for passing data between matlab and an external program. 
This is useful for performing a large scale computation in matlab which is programmed 
with an external library (probably with parallelization). In the examples here, we create 
a code with the Intel MKL library to do the SVD and pivoted QR decompositions and are able 
to call these functions from matlab using the socket interface. This is just an example, 
the socket interface can be used for other libraries too. 
This is useful because it is often difficult to build matlab 
mex files linked with external libraries (like MKL) and running in parallel. This code 
completely separates the mex file and the compute engine (the compute engine can even 
run on a different machine, so for example you can run a matlab code on a laptop and 
have the SVD computation inside be done for you on a local server).


========= how it works ==============================

Matlab passes data (matrix) to the wrapper. Wrapper m-file calls up the socket server 
with a randomly chosen port number between 1000 and 2000 (note: if you have a firewall, 
please adjust port numbers accordingly in the wrapper). Then the wrapper 
passes the data to the mex client and the mex client communicates with the socket server 
via sockets. The server does the computation and sends the results back to the mex client 
via sockets. The mex client receives results from the compute server and passes the results 
to the matlab wrapper which returns them to the matlab session or script from which it was called.
The reads and writes between the socket server and client are performed in a synchronous 
manner using a small amount of disk i/o (i.e. server sends client data when client is ready 
and vice versa). This is done to minimize possibility of errors with transfers.

========= summary of installation and usage =========

The following steps are for the use of the included SVD and QR examples. The examples use 
the Intel MKL library. The socket code can be used with different libraries also.

1) Install Intel C Compiler (icc) and Math Kernel Library (MKL) 

Follow instructions from Intel website (basically download install package and unpack)
Afterwards, to set paths, use (in bash shell):

source /opt/intel/bin/compilervars.sh intel64


2) Compile code (svd code from mkl library, socket client/server)

The compilation looks something like below for the SVD example, similar for QR. 
You must put the proper path to your matlab installation.

icc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I "/usr/local/MATLAB/R2013b/extern/include" svd_mex_client.c socket_functions.c -L"/usr/local/MATLAB/R2013b/bin/glnxa64" -Wl,--version-script,"/usr/local/MATLAB/R2013b/extern/lib/glnxa64/mexFunction.map"  -lmex -lmx -lmat -lm -o svd_mex_client.mexa64

icc -mkl -openmp svd_with_socket_server_intel_mkl.c socket_functions.c matrix_vector_functions_intel_mkl.c  -o svd_with_socket_server_intel_mkl

You can modify the compile.sh script with proper paths and run that to compile all the examples.

3) Run the code

First make a directory called logs. This is needed for i/o communication between 
server and client.

mkdir logs/

Open Matlab in same directory as the compiled code.
make sure intel mkl paths are set properly
i.e. source /opt/intel/bin/compilervars.sh intel64

Create matrix:

>> A = rand(2000,3000);

Run svd wrapper:

>> [U,S,Vt] = svd_wrapper_intel_mkl(A);

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

Run the QR wrapper:

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

