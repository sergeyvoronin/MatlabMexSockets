#include <stdio.h>
#include <time.h>



void error(char *msg);


int read_int_from_socket(int sockfd, int *element);


int write_int_to_socket(int sockfd, int *element);


int read_double_array_from_socket(int sockfd, double *elements, int size);

int write_double_array_to_socket(int sockfd, double *elements, int size);


void read_elements_from_socket(int sockfd, int nelements, double *elements_all);

void write_elements_to_socket(int sockfd, int nelements, double *elements_all);



void mark_ready(char *fname);


void mark_busy(char *fname);


int get_mark(char *fname);


void wait_for_ready(char *fname);


void mark_done_with_computation(char *fname);


void wait_for_done_with_computation(char *fname);

void sleep_funct();

void sleep_funct2();

