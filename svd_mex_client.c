/* max client for SVD with intel mkl library ; IPC via socket server 
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
    printf("broadcast m=%d and n=%d\n", m, n);
    wait_for_ready(fname_server);
    write_int_to_socket(sockfd, &m);
    write_int_to_socket(sockfd, &n);

    
    printf("transfer data\n");
    wait_for_ready(fname_server);
    write_elements_to_socket(sockfd, nelements, pr);
    printf("done transfer. doing computation..\n");

    

    // read dimensions
    int nrowsU,ncolsU,nrowsS,ncolsS,nrowsVt,ncolsVt;
    mark_ready(fname_client);
    read_int_from_socket(sockfd, &nrowsU);
    read_int_from_socket(sockfd, &ncolsU);
    read_int_from_socket(sockfd, &nrowsS);
    read_int_from_socket(sockfd, &ncolsS);
    read_int_from_socket(sockfd, &nrowsVt);
    read_int_from_socket(sockfd, &ncolsVt);
    mark_busy(fname_client);
    //write_int_to_socket(sockfd, &okstatus);


    printf("dimensions are:\n");
    printf("U: %d by %d\n", nrowsU, ncolsU);
    printf("S: %d by %d\n", nrowsS, ncolsS);
    printf("Vt: %d by %d\n", nrowsVt, ncolsVt);

    
    // read data
    double *Uarr, *Sarr, *Vtarr;

    // read U 
    sleep_funct();
    nelements = nrowsU*ncolsU;
    printf("read U which is %d by %d with nel = %d\n", nrowsU, ncolsU, nelements);
    Uarr = (double*)malloc(nelements*sizeof(double));
    printf("read elements\n");
    mark_ready(fname_client);
    read_elements_from_socket(sockfd, nelements, Uarr);
    mark_busy(fname_client);
    

    // read S 
    sleep_funct();
    nelements = nrowsS*ncolsS;
    printf("read S which is %d by %d with nel = %d\n", nrowsS, ncolsS, nelements);
    Sarr = (double*)malloc(nelements*sizeof(double));
    printf("read elements\n");
    mark_ready(fname_client);
    read_elements_from_socket(sockfd, nelements, Sarr); // read
    mark_busy(fname_client);
    

    // read Vt 
    sleep_funct();
    nelements = nrowsVt*ncolsVt;
    printf("read Vt which is %d by %d with nel = %d\n", nrowsVt, ncolsVt, nelements);
    Vtarr = (double*)malloc(nelements*sizeof(double));
    printf("read elements\n");
    mark_ready(fname_client);
    read_elements_from_socket(sockfd, nelements, Vtarr);
    mark_busy(fname_client);
    

    // create in Matlab
    printf("creating U: %d by %d..\n", nrowsU, ncolsU);
    dims[0] = nrowsU; dims[1] = ncolsU;
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    
    nelements = nrowsU*ncolsU;
    double *Uarr_matlab = (double*)malloc(nelements*sizeof(double));
    for(i=0; i<nelements; i++){
        Uarr_matlab[i] = Uarr[i];
    }
    ptr = mxGetPr(plhs[0]);
    bytes_to_copy = nrowsU*ncolsU*sizeof(double);
    memcpy(ptr,Uarr_matlab,bytes_to_copy);


    printf("creating S: %d by %d..\n", nrowsS, ncolsS);
    dims[0] = nrowsS; dims[1] = ncolsS;
    plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

    nelements = nrowsS*ncolsS;
    double *Sarr_matlab = (double*)malloc(nelements*sizeof(double));
    for(i=0; i<nelements; i++){
        Sarr_matlab[i] = Sarr[i];
    }
    ptr = mxGetPr(plhs[1]);
    bytes_to_copy = nrowsS*ncolsS*sizeof(double);
    memcpy(ptr,Sarr_matlab,bytes_to_copy);


    printf("creating Vt: %d by %d..\n", nrowsVt, ncolsVt);
    dims[0] = nrowsVt; dims[1] = ncolsVt;
    plhs[2] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

    nelements = nrowsVt*ncolsVt;
    double *Vtarr_matlab = (double*)malloc(nelements*sizeof(double));
    for(i=0; i<nelements; i++){
        Vtarr_matlab[i] = Vtarr[i];
    }
    ptr = mxGetPr(plhs[2]);
    bytes_to_copy = nrowsVt*ncolsVt*sizeof(double);
    memcpy(ptr,Vtarr_matlab,bytes_to_copy);


    printf("free and exit\n"); 

    // remove log file
    strcpy(cmd_str,"rm -f ");
    strcat(cmd_str,fname_client);
    system(cmd_str);

    close(sockfd);
    free(Uarr); free(Sarr); free(Vtarr);
    return;
}

