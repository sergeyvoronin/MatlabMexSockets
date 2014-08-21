Socket Code for Mex Interface
Includes example of SVD from Intel MKL library
Written by Sergey Voronin. 
Last updated: August 2014.
Note: this is the initial release (v. 0.1), need more error checking
Should work on any recent Linux distribution.
Tested on opensuse 13.1 with icc 14.03 and matlab R2013b.

========= what is this ==============================

This is a socket interface for passing data between matlab and an external program. 
This is useful for performing a large scale computation in matlab which is programmed 
with an external library (probably with parallelization). In the example here, we create 
a code with the Intel MKL library to do the SVD, but the socket interface can be used 
for other libraries too. This is useful because it is often difficult to build matlab 
mex files linked with external libraries (like MKL) and running in parallel. This code 
completely separates the mex file and the compute engine (the compute engine can even 
run on a different machine).


========= how it works ==============================

Matlab passes data (matrix) to the wrapper. Wrapper m-file calls up the socket server 
with a randomly chosen port number between 1000 and 2000 (note: if you have a firewall, 
please adjust port numbers accordingly in the wrapper). Then the wrapper 
passes the data to the mex client and the mex client communicates with the socket server 
via sockets. The server does the computation and sends the results back to the mex client 
via sockets. The mex client receives results from the compute server and passes the results 
to the matlab wrapper which returns them to the matlab session or script from which it was called.


========= summary of installation and usage =========


1) Install Intel MKL (for this example)

Follow instructions from Intel website (basically download install package and unpack)

source /opt/intel/bin/compilervars.sh intel64


2) Compile code (svd code from mkl library, socket client/server)

Basically, the compilation looks something like below. You must put the proper path to your 
matlab installation.

icc -fpic -shared -DMATLAB_MEX_FILE -fno-omit-frame-pointer -pthread -I "/usr/local/MATLAB/R2013b/extern/include" svd_mex_client.c socket_functions.c -L"/usr/local/MATLAB/R2013b/bin/glnxa64" -Wl,--version-script,"/usr/local/MATLAB/R2013b/extern/lib/glnxa64/mexFunction.map"  -lmex -lmx -lmat -lm -o svd_mex_client.mexa64

icc -mkl -openmp svd_with_socket_server_intel_mkl.c socket_functions.c matrix_vector_functions_intel_mkl.c  -o svd_with_socket_server_intel_mkl


3) Run the code

Open Matlab in same directory
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

>> 
