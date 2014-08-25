/* max client for full rank QR ; IPC via socket server 
Sergey Voronin, 2014
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 
#include <strings.h>
#include <string.h>
#include <time.h>
#include "matrix.h"
#include "mex.h"
#include "socket_functions.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int i, j, m, n, k, T_solver_num, info, solve_type;
    int sockfd, portno, nlen, nelements, chunk_num, ind, nchunks, nleftover;
    int nrows, ncols;
    time_t start_time, end_time;
    double *ptr;
    size_t bytes_to_copy;
    mwSize dims[2];

    struct sockaddr_in serv_addr;
    struct hostent *server;
    double  *pr, *elements_local;
    FILE *fp;
    char buffer[256];
    
    int sndbuf, rcvbuf;
    int okstatus = 1;
    char *buf = (char*)malloc(10*sizeof(char));
    char *fname_server = (char*)malloc(100*sizeof(char));
    char *fname_client = (char*)malloc(100*sizeof(char));
    char *cmd_str = (char*)malloc(100*sizeof(char));
    struct timespec tspec;
 

    /* Check for proper number of arguments. */
    if (nrhs != 2) 
        mexErrMsgTxt("Two inputs required: the matrix and the port number of the server.\n");

    /* get input parameters */
    n = mxGetN(prhs[0]);
    m = mxGetM(prhs[0]);
    pr = mxGetPr(prhs[0]);
    portno = *mxGetPr(prhs[1]);

    printf("the matrix is of size %d x %d and the port number is %d\n", m, n, portno);

    // set up server and client files
    sprintf(buf,"%d",portno);
    strcpy(fname_server,"logs/server_file_");
    strcat(fname_server,buf);
    strcat(fname_server,".txt");
    
    strcpy(fname_client,"logs/client_file_");
    strcat(fname_client,buf);
    strcat(fname_client,".txt");

    printf("using log files: %s and %s\n", fname_server, fname_client);
    mark_busy(fname_client);
    


    // prepare to broadcast this stuff
    printf("prepare broadcast\n");
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if(sockfd < 0) 
        error("ERROR opening socket");
    server = gethostbyname("localhost");
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr, 
         (char *)&serv_addr.sin_addr.s_addr,
         server->h_length);
    serv_addr.sin_port = htons(portno);
    if (connect(sockfd,(struct sockaddr *)&serv_addr,sizeof(serv_addr)) < 0){
        error("ERROR connecting");
        return;
    }

    // set socket options if needed 
    //sndbuf = 50000;  /* Send buffer size */
    //rcvbuf = 50000;  /* receive buffer size */
    //nlen = setsockopt(sockfd,SOL_SOCKET,SO_SNDBUF, &sndbuf, sizeof(sndbuf));
    //nlen = setsockopt(sockfd,SOL_SOCKET,SO_RCVBUF, &rcvbuf, sizeof(rcvbuf));


    nelements = m*n;

    // broadcast
    printf("broadcast m=%d, n=%d\n", m, n);
    wait_for_ready(fname_server);
    write_int_to_socket(sockfd, &m);
    write_int_to_socket(sockfd, &n);

    
    printf("transfer data\n");
    wait_for_ready(fname_server);
    write_elements_to_socket(sockfd, nelements, pr);
    printf("done transfer. doing computation..\n");
    sleep_funct();
    wait_for_done_with_computation(fname_server);
    printf("done with computation.\n");

    
    // read dimensions
    int nrowsQk,ncolsQk,nrowsRk,ncolsRk,nrowsP,ncolsP,nrowsI;
    mark_ready(fname_client);
    read_int_from_socket(sockfd, &nrowsQk);
    read_int_from_socket(sockfd, &ncolsQk);
    read_int_from_socket(sockfd, &nrowsRk);
    read_int_from_socket(sockfd, &ncolsRk);
    read_int_from_socket(sockfd, &nrowsP);
    read_int_from_socket(sockfd, &ncolsP);
    read_int_from_socket(sockfd, &nrowsI);
    mark_busy(fname_client);


    printf("dimensions are:\n");
    printf("Qk: %d by %d\n", nrowsQk, ncolsQk);
    printf("Rk: %d by %d\n", nrowsRk, ncolsRk);
    printf("P: %d by %d\n", nrowsP, ncolsP);
    printf("I: %d by %d\n", nrowsI, 1);

    
    // read data
    double *Qkarr, *Rkarr, *Parr, *Iarr;

    // read Qk 
    sleep_funct();
    nelements = nrowsQk*ncolsQk;
    printf("read Qk which is %d by %d with nel = %d\n", nrowsQk, ncolsQk, nelements);
    Qkarr = (double*)malloc(nelements*sizeof(double));
    printf("read elements\n");
    mark_ready(fname_client);
    read_elements_from_socket(sockfd, nelements, Qkarr);
    mark_busy(fname_client);
    

    // read Rk 
    sleep_funct();
    nelements = nrowsRk*ncolsRk;
    printf("read Rk which is %d by %d with nel = %d\n", nrowsRk, ncolsRk, nelements);
    Rkarr = (double*)malloc(nelements*sizeof(double));
    printf("read elements\n");
    mark_ready(fname_client);
    read_elements_from_socket(sockfd, nelements, Rkarr); // read
    mark_busy(fname_client);
    

    // read P 
    sleep_funct();
    nelements = nrowsP*ncolsP;
    printf("read P which is %d by %d with nel = %d\n", nrowsP, ncolsP, nelements);
    Parr = (double*)malloc(nelements*sizeof(double));
    printf("read elements\n");
    mark_ready(fname_client);
    read_elements_from_socket(sockfd, nelements, Parr);
    mark_busy(fname_client);


    // read I (vector) 
    sleep_funct();
    nelements = nrowsI;
    printf("read I which is %d by %d with nel = %d\n", nrowsI, 1, nelements);
    Iarr = (double*)malloc(nelements*sizeof(double));
    printf("read elements\n");
    mark_ready(fname_client);
    read_elements_from_socket(sockfd, nelements, Iarr);
    mark_busy(fname_client);
    
    

    // create in Matlab
    printf("creating Qk: %d by %d..\n", nrowsQk, ncolsQk);
    dims[0] = nrowsQk; dims[1] = ncolsQk;
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    
    nelements = nrowsQk*ncolsQk;
    double *Qkarr_matlab = (double*)malloc(nelements*sizeof(double));
    for(i=0; i<nelements; i++){
        Qkarr_matlab[i] = Qkarr[i];
    }
    ptr = mxGetPr(plhs[0]);
    bytes_to_copy = nrowsQk*ncolsQk*sizeof(double);
    memcpy(ptr,Qkarr_matlab,bytes_to_copy);


    printf("creating Rk: %d by %d..\n", nrowsRk, ncolsRk);
    dims[0] = nrowsRk; dims[1] = ncolsRk;
    plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

    nelements = nrowsRk*ncolsRk;
    double *Rkarr_matlab = (double*)malloc(nelements*sizeof(double));
    for(i=0; i<nelements; i++){
        Rkarr_matlab[i] = Rkarr[i];
    }
    ptr = mxGetPr(plhs[1]);
    bytes_to_copy = nrowsRk*ncolsRk*sizeof(double);
    memcpy(ptr,Rkarr_matlab,bytes_to_copy);


    printf("creating P: %d by %d..\n", nrowsP, ncolsP);
    dims[0] = nrowsP; dims[1] = ncolsP;
    plhs[2] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

    nelements = nrowsP*ncolsP;
    double *Parr_matlab = (double*)malloc(nelements*sizeof(double));
    for(i=0; i<nelements; i++){
        Parr_matlab[i] = Parr[i];
    }
    ptr = mxGetPr(plhs[2]);
    bytes_to_copy = nrowsP*ncolsP*sizeof(double);
    memcpy(ptr,Parr_matlab,bytes_to_copy);


    printf("creating I: %d by %d..\n", nrowsI, 1);
    dims[0] = nrowsI; dims[1] = 1;
    plhs[3] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

    nelements = nrowsI;
    double *Iarr_matlab = (double*)malloc(nelements*sizeof(double));
    for(i=0; i<nelements; i++){
        Iarr_matlab[i] = Iarr[i];
    }
    ptr = mxGetPr(plhs[3]);
    bytes_to_copy = nrowsI*sizeof(double);
    memcpy(ptr,Iarr_matlab,bytes_to_copy);


    printf("free and exit\n"); 

    // remove log file
    strcpy(cmd_str,"rm -f ");
    strcat(cmd_str,fname_client);
    system(cmd_str);


    close(sockfd);
    free(Qkarr); free(Rkarr); free(Parr); free(Iarr);
    return;
}

