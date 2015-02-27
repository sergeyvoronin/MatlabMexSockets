/* basic socket read / write functions. Sergey Voronin */

#include "socket_functions.h"

void error(char *msg){
    perror(msg);
}


int read_int_from_socket(int sockfd, int *element){
    int total = 0;        // how many bytes we've read 
    int bytesremaining = sizeof(int); // how many we have left to read
    int n;

    n = read(sockfd, element, bytesremaining);
    if(n!=bytesremaining){
        total += n;
        bytesremaining -= n;

        while(bytesremaining > 0){
            n = read(sockfd, element+total, bytesremaining);
            if (n == -1) { break; }
            total += n;
            bytesremaining -= n;
        }
    }
    return n;
}


int write_int_to_socket(int sockfd, int *element){
    int total = 0;        // how many bytes we've sent
    int bytesremaining = sizeof(int); // how many we have left to send
    int n;

    n = write(sockfd, element, bytesremaining);
    if(n!=bytesremaining){
        total += n;
        bytesremaining -= n;

        while(bytesremaining > 0) {
            n = write(sockfd, element+total, bytesremaining);
            if (n == -1) { break; }
            total += n;
            bytesremaining -= n;
        }
    }
    return n;
}




int read_double_array_from_socket(int sockfd, double *elements, int size){
    int total = 0;        // how many bytes we've read
    int bytesremaining = size; // how many we have left to read
    int n;

    n = read(sockfd, elements, size);
    if(n!=size){
        if(n>0){
            total += n;
            bytesremaining -= n;
        }

        while(bytesremaining>0) {
            n = read(sockfd, elements+total, bytesremaining);
            if (n == -1) { break; }
            total += n;
            bytesremaining -= n;
        }
    }
    return n;
}



int write_double_array_to_socket(int sockfd, double *elements, int size){
    int total = 0;        // how many bytes we've sent
    int bytesremaining = size; // how many we have left to send
    int n;

    n = write(sockfd, elements, size);
    if(n!=size){
        if(n>0){
            total += n;
            bytesremaining -= n;
        }

        while(bytesremaining>0) {
            n = write(sockfd, elements+total, bytesremaining);
            if (n == -1) { break; }
            total += n;
            bytesremaining -= n;
        }
    }
    return n;
}




void read_elements_from_socket(int sockfd, int nelements, double *elements_all){
    int i, ind, nlen, nchunks, nleftover, chunk_num;
    double * elements_local;

    // read all elements
    if(nelements <= 6000){
        nlen = read_double_array_from_socket(sockfd, elements_all, nelements*sizeof(double));
        if(nlen < 0) 
            error("ERROR writing to socket");
    }
    // read in batches of 6000
    else{
        nchunks = nelements/6000;  
        nleftover = nelements % 6000;
        printf("server: nchunks = %d and nleftover = %d\n", nchunks, nleftover);
        ind = 0; // ind into global array

        for(chunk_num=0; chunk_num < nchunks; chunk_num++){
                //printf("read chunk %d\n", chunk_num);
                elements_local = (double*)malloc(6000*sizeof(double));
                nlen = read_double_array_from_socket(sockfd, elements_local, 6000*sizeof(double));
                for(i=0; i<6000; i++){
                    elements_all[ind++] = elements_local[i];
                }
                free(elements_local);
        }
        //printf("read leftovers..\n");
        elements_local = (double*)malloc(nleftover*sizeof(double));
        nlen = read_double_array_from_socket(sockfd, elements_local, nleftover*sizeof(double));
        for(i=0; i<nleftover; i++){
            elements_all[ind++] = elements_local[i];
        }
        free(elements_local);
    }
}




void write_elements_to_socket(int sockfd, int nelements, double *elements_all){
    int i, ind, nlen, nchunks, nleftover, chunk_num;
    double * elements_local;

    // write all elements
    if(nelements <= 6000){
        printf("transfer all at once..\n");
        //nlen = write(sockfd,elements_all,nelements*sizeof(double));
        nlen = write_double_array_to_socket(sockfd, elements_all, nelements*sizeof(double));
        if(nlen < 0) 
            error("ERROR writing to socket");
    }
    // write in batches of 6000
    else{
        nchunks = nelements/6000;  
        nleftover = nelements % 6000;
        printf("server: nchunks = %d and nleftover = %d\n", nchunks, nleftover);

        ind = 0; // ind into global array
        for(chunk_num=0; chunk_num < nchunks; chunk_num++){
                //printf("transfer chunk %d\n", chunk_num);
                elements_local = (double*)malloc(6000*sizeof(double));
                for(i=0; i<6000; i++){
                    elements_local[i] = elements_all[ind++];
                }
                //nlen = write(sockfd,elements_local,6000*sizeof(double));
                nlen = write_double_array_to_socket(sockfd, elements_local, 6000*sizeof(double));
                sleep_funct2();
                free(elements_local);
                
        }
        //printf("transfer leftovers...\n");
        elements_local = (double*)malloc(nleftover*sizeof(double));
        for(i=0; i<nleftover; i++){
            elements_local[i] = elements_all[ind++];
        }
        //nlen = write(sockfd,elements_local,nleftover*sizeof(double));
        nlen = write_double_array_to_socket(sockfd, elements_local, nleftover*sizeof(double));
        free(elements_local);
    }
}



void mark_ready(char *fname){
    FILE *fp;
    fp = fopen(fname,"w");
    fprintf(fp,"1\n");
    fclose(fp); 
}


void mark_done_with_computation(char *fname){
    FILE *fp;
    fp = fopen(fname,"w");
    fprintf(fp,"2\n");
    fclose(fp); 
}


void mark_busy(char *fname){
    FILE *fp;
    fp = fopen(fname,"w");
    fprintf(fp,"0\n");
    fclose(fp); 
}


int get_mark(char *fname){
    int mark;
    FILE *fp;
    char *buf = (char*)malloc(10*sizeof(char)); 
    fp = fopen(fname,"r");
    fgets(buf,10,fp);
    mark = atoi(buf);
    fclose(fp);
    return mark; 
}


void wait_for_ready(char *fname){
    int mark = 0;
    while(mark < 1){
        mark = get_mark(fname); 
        sleep(1);
    }
}


void wait_for_done_with_computation(char *fname){
    int mark = 0;
    while(mark < 2){
        mark = get_mark(fname); 
        sleep(1);
    }
}



void sleep_funct(){
    nanosleep((struct timespec[]){{0, 360000000}}, NULL);
}

void sleep_funct2(){
    nanosleep((struct timespec[]){{0, 2600000}}, NULL);
}

