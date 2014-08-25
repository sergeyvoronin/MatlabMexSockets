/* full QR of rank k (via MKL functions) callable from matlab with Mex Socket Server
Sergey Voronin - 2014 
*/

#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include <strings.h>
#include "socket_functions.h"
#include "matrix_vector_functions_intel_mkl.h"




/* QR for mxn A with m > n */
void pivotedQR_mkl(mat *M, mat **Q, mat **R, mat **P, vec **I){
    int i,j,k,m,n,Rrows,Rcols,Qrows,Qcols;
    mat *Mwork, *Porig;
    vec *col_vec;
    m = M->nrows; n = M->ncols;
    k = min(m,n);
    
    Mwork = matrix_new(m,n);
    matrix_copy(Mwork,M);

    int *Iarr = (int*)malloc(n*sizeof(int));
    double *tau_arr = (double*)malloc(min(m,n)*sizeof(double));

    // set up dimensions
    if(m <= n){
        Qrows = k; Qcols = k;
        Rrows = k; Rcols = n;
    } else {
        Qrows = m; Qcols = k;
        Rrows = k; Rcols = k;
    }
    
    // get R
    LAPACKE_dgeqp3(CblasColMajor, Mwork->nrows, Mwork->ncols, Mwork->d, Mwork->nrows, Iarr, tau_arr);

    *R = matrix_new(Rrows,Rcols);
    for(i=0; i<Rrows; i++){
        for(j=i; j<Rcols; j++){
            matrix_set_element(*R,i,j,matrix_get_element(Mwork,i,j));
        }
    }

    // get Q
    LAPACKE_dorgqr(CblasColMajor, Mwork->nrows, Mwork->nrows, min(Mwork->nrows,Mwork->ncols), Mwork->d, Mwork->nrows, tau_arr);

    *Q = matrix_new(Qrows,Qcols);
    for(i=0; i<Qrows; i++){
        for(j=0; j<Qcols; j++){
            matrix_set_element(*Q,i,j,matrix_get_element(Mwork,i,j));
        }
    }

    // get I
    *I = vector_new(n);
    for(i=0; i<n; i++){
        vector_set_element(*I,i,Iarr[i]);
    }

    // get P
    Porig = matrix_new(n,n);
    *P = matrix_new(n,n);
    initialize_identity_matrix(Porig);
    for(i=0; i<n; i++){
        col_vec = vector_new(n);
        matrix_get_col(Porig,vector_get_element(*I,i)-1,col_vec);
        matrix_set_col(*P,i,col_vec);
        vector_delete(col_vec);
    }

    matrix_delete(Porig);
}


int main(int argc, char *argv[])
{
    int i, j, m, n, frank, sockfd, newsockfd, portno, clilen, nelements, stop = 0;
    int nlen, number, nchunks, nleftover, ind, chunk_num, nrows, ncols;
    double normM,normU,normS,normV,normP,percent_error;
    mat *M, *Qk, *Rk, *P;
    vec *I;
    time_t start_time, end_time;
    struct sockaddr_in serv_addr, cli_addr;
    double *elements_all;
    FILE *fp, *fp_server, *fp_client;

    int sndbuf, rcvbuf;
    int okstatus;
    char * buf = (char*)malloc(10*sizeof(char));
    char * fname_server = (char*)malloc(100*sizeof(char));
    char * fname_client = (char*)malloc(100*sizeof(char));
    char * cmd_str = (char*)malloc(200*sizeof(char));
    
    if (argc < 2) {
         fprintf(stderr,"ERROR, no port provided\n");
         return(1);
    }

    
    printf("setting up socket server..\n");
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0){ 
        error("ERROR opening socket");
        close(sockfd);
        return -1;
    }
    else{
        printf("socket initialized\n");
    }
    bzero((char *) &serv_addr, sizeof(serv_addr));
    portno = atoi(argv[1]);
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port = htons(portno);
    if(bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0){ 
        error("ERROR on binding");
        close(sockfd);
        return -1;
    }
    listen(sockfd,5);
    clilen = sizeof(cli_addr);
    printf("binding complete\n");

    // set up server and client files
    strcpy(fname_server,"logs/server_file_");
    strcat(fname_server,argv[1]);
    strcat(fname_server,".txt");
    
    strcpy(fname_client,"logs/client_file_");
    strcat(fname_client,argv[1]);
    strcat(fname_client,".txt");
    
    mark_busy(fname_server);


    printf("looping socket server.\n");
    while(!stop){
        newsockfd = accept(sockfd, (struct sockaddr *) &cli_addr, &clilen);
        printf("after accept\n");
        if (newsockfd < 0) 
          error("ERROR on accept");

        // if socket performance is too slow, experiment with these socket options
        //sndbuf = 50000;  /* send buffer size */
        //rcvbuf = 50000;  /* receive buffer size */
        //nlen = setsockopt(newsockfd,SOL_SOCKET,SO_SNDBUF, &sndbuf, sizeof(sndbuf));
        //nlen = setsockopt(newsockfd,SOL_SOCKET,SO_RCVBUF, &rcvbuf, sizeof(rcvbuf));

        mark_ready(fname_server);
        read_int_from_socket(newsockfd, &m);
        read_int_from_socket(newsockfd, &n);
        mark_busy(fname_server);

        nelements = m*n;
        printf("read nrows = %d and ncols = %d\n", m, n);


        // init space for all elements 
        elements_all = (double*)malloc(nelements*sizeof(double));


        // read all elements
        mark_ready(fname_server);
        read_elements_from_socket(newsockfd, nelements, elements_all);
        mark_busy(fname_server);

        
        // set up matrix using the read data
        printf("setup matrix of size %d by %d\n", m,n);
        M = matrix_new(m,n);
        M->d = elements_all;

        printf("call full QR..\n");
        pivotedQR_mkl(M, &Qk, &Rk, &P, &I);
        matrix_hard_threshold(Rk,1e-12);
        mark_done_with_computation(fname_server);
        printf("finished QR.\n");

        // record info
        /*fp = fopen("logs/qr_info.txt", "w");
        fprintf(fp,"norm(Qk,fro) = %f\n", get_matrix_frobenius_norm(Qk));
        fprintf(fp,"norm(Rk,fro) = %f\n", get_matrix_frobenius_norm(Rk));
        fclose(fp);*/

        // pass dimensions
        wait_for_ready(fname_client);
        write_int_to_socket(newsockfd, &Qk->nrows);
        write_int_to_socket(newsockfd, &Qk->ncols);
        write_int_to_socket(newsockfd, &Rk->nrows);
        write_int_to_socket(newsockfd, &Rk->ncols);
        write_int_to_socket(newsockfd, &P->nrows);
        write_int_to_socket(newsockfd, &P->ncols);
        write_int_to_socket(newsockfd, &I->nrows);

        // write Qk
        //sleep_funct();
        nrows = Qk->nrows; ncols = Qk->ncols; 
        nelements = nrows*ncols;
        wait_for_ready(fname_client);
        write_elements_to_socket(newsockfd, nelements, Qk->d);

       
        // write Rk
        //sleep_funct();
        nrows = Rk->nrows; ncols = Rk->ncols; 
        nelements = nrows*ncols;
        wait_for_ready(fname_client);
        write_elements_to_socket(newsockfd, nelements, Rk->d);
        

        // write P
        //sleep_funct();
        nrows = P->nrows; ncols = P->ncols; 
        nelements = nrows*ncols;
        wait_for_ready(fname_client);
        write_elements_to_socket(newsockfd, nelements, P->d);
        

        // write I
        //sleep_funct();
        nrows = I->nrows; ncols = 1; 
        nelements = nrows;
        wait_for_ready(fname_client);
        write_elements_to_socket(newsockfd, nelements, I->d);
        

        // close
        if (close(newsockfd) != 0)
            printf("server newsockfd closing failed!\n");
        else
            printf("server newsockfd successfully closed!\n");

        stop = 1;
    // end of socket loop
    }
   

    printf("done with socket server.");

    // remove log file
    strcpy(cmd_str,"rm -f ");
    strcat(cmd_str,fname_server);
    system(cmd_str);

    if (close(sockfd) != 0)
        printf("server sockfd closing failed!\n");
    else
        printf("server sockfd successfully closed!\n");

    matrix_delete(M);
    matrix_delete(Qk);
    matrix_delete(Rk);
    matrix_delete(P);
    vector_delete(I);
    
    return 0;
}

