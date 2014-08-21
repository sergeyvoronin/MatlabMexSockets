/* SVD with Intel MKL callable from matlab with Mex Socket Server
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



/* computes SVD: M = U*S*Vt; note Vt = V^T */
void singular_value_decomposition(mat *M, mat **U, mat **S, mat **Vt){
    int m,n,k;
    m = M->nrows; n = M->ncols;
    k = min(m,n);
    vec * work = vector_new(2*max(3*min(m, n)+max(m, n), 5*min(m,n)));
    vec * svals = vector_new(k);
    *U = matrix_new(m,k);
    *S = matrix_new(k,k);
    *Vt = matrix_new(k,n);

    LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'S', 'S', m, n, M->d, m, svals->d, (*U)->d, m, (*Vt)->d, k, work->d );

    initialize_diagonal_matrix(*S, svals);

    vector_delete(work);
    vector_delete(svals);
}


int main(int argc, char *argv[])
{
    int i, j, m, n, k, frank, T_solver_num, sockfd, newsockfd, portno, clilen, nelements, stop = 0;
    int nlen, number, nchunks, nleftover, ind, chunk_num, nrows, ncols;
    double normM,normU,normS,normV,normP,percent_error;
    mat *M, *U, *S, *Vt;
    vec *I;
    time_t start_time, end_time;
    struct sockaddr_in serv_addr, cli_addr;
    double *elements_all, *elements_local;
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

        printf("call SVD..\n");
        singular_value_decomposition(M, &U, &S, &Vt);
        printf("finished SVD.\n");

        // optionally, record info
        /*fp = fopen("svd_info.txt", "w");
        fprintf(fp,"norm(U,fro) = %f\n", get_matrix_frobenius_norm(U));
        fprintf(fp,"norm(S,fro) = %f\n", get_matrix_frobenius_norm(S));
        fprintf(fp,"norm(Vt,fro) = %f\n", get_matrix_frobenius_norm(Vt));
        fclose(fp);*/

        // pass dimensions
        wait_for_ready(fname_client);
        write_int_to_socket(newsockfd, &U->nrows);
        write_int_to_socket(newsockfd, &U->ncols);
        write_int_to_socket(newsockfd, &S->nrows);
        write_int_to_socket(newsockfd, &S->ncols);
        write_int_to_socket(newsockfd, &Vt->nrows);
        write_int_to_socket(newsockfd, &Vt->ncols);
        //read_int_from_socket(newsockfd, &okstatus);

        // write U
        //sleep_funct();
        nrows = U->nrows; ncols = U->ncols; 
        nelements = nrows*ncols;
        wait_for_ready(fname_client);
        write_elements_to_socket(newsockfd, nelements, U->d);
        //read_int_from_socket(newsockfd, &okstatus);

       
        // write S
        //sleep_funct();
        nrows = S->nrows; ncols = S->ncols; 
        nelements = nrows*ncols;
        wait_for_ready(fname_client);
        write_elements_to_socket(newsockfd, nelements, S->d);
        //read_int_from_socket(newsockfd, &okstatus);
        

        // write Vt
        //sleep_funct();
        nrows = Vt->nrows; ncols = Vt->ncols; 
        nelements = nrows*ncols;
        wait_for_ready(fname_client);
        write_elements_to_socket(newsockfd, nelements, Vt->d);
        //read_int_from_socket(newsockfd, &okstatus);
        

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
    
    return 0;
}

