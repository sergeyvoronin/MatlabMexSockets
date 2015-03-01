/*
 * cula gpu code with host openmp 
 */

#include "matrix_vector_functions_nvidia_cula.h"


/* initialize new matrix and set all entries to zero */
mat * matrix_new(int nrows, int ncols)
{
    mat *M = malloc(sizeof(mat));
    //M->d = (double*)mkl_malloc(sizeof(double) * nrows*ncols, 64);
    M->d = (double*)calloc(nrows*ncols, sizeof(double));
    M->nrows = nrows;
    M->ncols = ncols;
    return M;
}


/* initialize new vector and set all entries to zero */
vec * vector_new(int nrows)
{
    vec *v = malloc(sizeof(vec));
    v->d = (double*)calloc(nrows,sizeof(double));
    v->nrows = nrows;
    return v;
}


void matrix_delete(mat *M)
{
    free(M->d);
    free(M);
}


void vector_delete(vec *v)
{
    free(v->d);
    free(v);
}


// column major format
void matrix_set_element(mat *M, int row_num, int col_num, double val){
    //M->d[row_num*(M->ncols) + col_num] = val;
    M->d[col_num*(M->nrows) + row_num] = val;
}



double matrix_get_element(mat *M, int row_num, int col_num){
    //return M->d[row_num*(M->ncols) + col_num];
    return M->d[col_num*(M->nrows) + row_num];
}


void vector_set_element(vec *v, int row_num, double val){
    v->d[row_num] = val;
}


double vector_get_element(vec *v, int row_num){
    return v->d[row_num];
}


void vector_set_data(vec *v, double *data){
    int i;
    #pragma omp parallel shared(v) private(i) 
    {
    #pragma omp for
    for(i=0; i<(v->nrows); i++){
        v->d[i] = data[i];
    }
    }
}


/* scale vector by a constant */
void vector_scale(vec *v, double scalar){
    int i;
    #pragma omp parallel shared(v,scalar) private(i) 
    {
    #pragma omp for
    for(i=0; i<(v->nrows); i++){
        v->d[i] = scalar*(v->d[i]);
    }
    }
}


/* compute euclidean norm of vector */
double vector_get2norm(vec *v){
    int i;
    double val, normval = 0;
    #pragma omp parallel shared(v,normval) private(i,val) 
    {
    #pragma omp for reduction(+:normval)
    for(i=0; i<(v->nrows); i++){
        val = v->d[i];
        normval += val*val;
    }
    }
    return sqrt(normval);
}


/* copy contents of vec s to d  */
void vector_copy(vec *d, vec *s){
    int i;
    //#pragma omp parallel for
    #pragma omp parallel shared(d,s) private(i) 
    {
    #pragma omp for 
    for(i=0; i<(s->nrows); i++){
        d->d[i] = s->d[i];
    }
    }
}




/* extract column of a matrix into a vector */
void matrix_get_col(mat *M, int j, vec *column_vec){
    int i;
    #pragma omp parallel shared(column_vec,M,j) private(i) 
    {
    #pragma omp parallel for
    for(i=0; i<M->nrows; i++){
        vector_set_element(column_vec,i,matrix_get_element(M,i,j));
    }
    }
}


/* set column of matrix to vector */
void matrix_set_col(mat *M, int j, vec *column_vec){
    int i;
    #pragma omp parallel shared(column_vec,M,j) private(i) 
    {
    #pragma omp for
    for(i=0; i<M->nrows; i++){
        matrix_set_element(M,i,j,vector_get_element(column_vec,i));
    }
    }
}


/* extract row i of a matrix into a vector */
void matrix_get_row(mat *M, int i, vec *row_vec){
    int j;
    #pragma omp parallel shared(row_vec,M,i) private(j) 
    {
    #pragma omp parallel for
    for(j=0; j<M->ncols; j++){ 
        vector_set_element(row_vec,j,matrix_get_element(M,i,j));
    }
    }
}


/* put vector row_vec as row i of a matrix */
void matrix_set_row(mat *M, int i, vec *row_vec){
    int j;
    #pragma omp parallel shared(row_vec,M,i) private(j) 
    {
    #pragma omp parallel for
    for(j=0; j<M->ncols; j++){ 
        matrix_set_element(M,i,j,vector_get_element(row_vec,j));
    }
    }
}



/* copy contents of mat S to D  */
void matrix_copy(mat *D, mat *S){
    int i;
    //#pragma omp parallel for
    #pragma omp parallel shared(D,S) private(i) 
    {
    #pragma omp for 
    for(i=0; i<((S->nrows)*(S->ncols)); i++){
        D->d[i] = S->d[i];
    }
    }
}



/* copy the first k rows of M into M_out where k = M_out->nrows (M_out pre-initialized) */
void matrix_copy_first_rows(mat *M_out, mat *M){
    int i,k;
    k = M_out->nrows;
    vec * row_vec;
    for(i=0; i<k; i++){
        row_vec = vector_new(M->ncols);
        matrix_get_row(M,i,row_vec);
        matrix_set_row(M_out,i,row_vec);
        vector_delete(row_vec);
    }
} 


/* copy the first k columns of M into M_out where k = M_out->ncols (M_out pre-initialized) */
void matrix_copy_first_columns(mat *M_out, mat *M){
    int i,k;
    k = M_out->ncols;
    vec * col_vec;
    for(i=0; i<k; i++){
        col_vec = vector_new(M->nrows);
        matrix_get_col(M,i,col_vec);
        matrix_set_col(M_out,i,col_vec);
        vector_delete(col_vec);
    }
} 




/* load matrix from file 
format:
% comment
num_rows num_columns num_nonzeros
row col nnz
.....
row col nnz
*/
mat * matrix_load_from_text_file(char *fname){
    int i, j, num_rows, num_columns, num_nonzeros, row_num, col_num;
    double nnz_val;
    char *nnz_val_str;
    char *line;
    FILE *fp;
    mat *M;
    
    line = (char*)malloc(200*sizeof(char));
    fp = fopen(fname,"r");
    fgets(line,100,fp); //read comment
    fgets(line,100,fp); //read dimensions and nnzs 
    sscanf(line, "%d %d %d", &num_rows, &num_columns, &num_nonzeros);
    M = matrix_new(num_rows,num_columns);

    // read and set elements
    nnz_val_str = (char*)malloc(50*sizeof(char));
    for(i=0; i<(num_nonzeros); i++){
        fgets(line,100,fp); 
        sscanf(line, "%d %d %s", &row_num, &col_num, nnz_val_str);
        nnz_val = atof(nnz_val_str);
        matrix_set_element(M,row_num,col_num,nnz_val);
    }
    fclose(fp);

    // clean
    free(line);
    free(nnz_val_str);

    return M;
}



/* load matrix from binary file 
 * the nonzeros are in order of double loop over rows and columns
format:
num_rows (int) 
num_columns (int)
nnz (double)
...
nnz (double)
*/
mat * matrix_load_from_binary_file(char *fname){
    int i, j, num_rows, num_columns, row_num, col_num;
    double nnz_val;
    size_t one = 1;
    FILE *fp;
    mat *M;
    
    fp = fopen(fname,"r");
    fread(&num_rows,sizeof(int),one,fp); //read m
    fread(&num_columns,sizeof(int),one,fp); //read n
    printf("initializing M of size %d by %d\n", num_rows, num_columns);
    M = matrix_new(num_rows,num_columns);
    printf("done..\n");

    // read and set elements
    for(i=0; i<num_rows; i++){
        for(j=0; j<num_columns; j++){
            fread(&nnz_val,sizeof(double),one,fp); //read nnz
            matrix_set_element(M,i,j,nnz_val);
        }
    }
    fclose(fp);

    return M;
}


/* write matrix to binary file 
 * the nonzeros are in order of double loop over rows and columns
format:
num_rows (int) 
num_columns (int)
nnz (double)
...
nnz (double)
*/
void matrix_write_to_binary_file(mat *M, char *fname){
    int i, j, num_rows, num_columns, row_num, col_num;
    double nnz_val;
    size_t one = 1;
    FILE *fp;
    num_rows = M->nrows; num_columns = M->ncols;
    
    fp = fopen(fname,"w");
    fwrite(&num_rows,sizeof(int),one,fp); //write m
    fwrite(&num_columns,sizeof(int),one,fp); //write n

    // write the elements
    for(i=0; i<num_rows; i++){
        for(j=0; j<num_columns; j++){
            nnz_val = matrix_get_element(M,i,j);
            fwrite(&nnz_val,sizeof(double),one,fp); //write nnz
        }
    }
    fclose(fp);
}



/* load vector from file 
format:
% comment
num_rows
value
.....
value
*/
vec * vector_load_from_file(char *fname){
    int i, j, num_rows;
    double nnz_val;
    char *nnz_val_str;
    char *line;
    FILE *fp;
    vec *v;

    line = (char*)malloc(200*sizeof(char));
    fp = fopen(fname,"r");
    fgets(line,100,fp); //read comment
    fgets(line,100,fp); //read dimension 
    sscanf(line, "%d", &num_rows);
    v = vector_new(num_rows);

    // read and set elements
    nnz_val_str = (char*)malloc(50*sizeof(char));
    for(i=0; i<num_rows; i++){
        fgets(line,100,fp);
        sscanf(line, "%s", nnz_val_str);
        nnz_val = atof(nnz_val_str);
        vector_set_element(v, i, nnz_val);
    }
    fclose(fp);

    // clean
    free(line);
    free(nnz_val_str);

    return v;
}



void matrix_print(mat * M){
    int i,j;
    double val;
    for(i=0; i<M->nrows; i++){
        for(j=0; j<M->ncols; j++){
            val = matrix_get_element(M, i, j);
            printf("%f  ", val);
        }
        printf("\n");
    }
}


void vector_print(vec * v){
    int i;
    double val;
    for(i=0; i<v->nrows; i++){
        val = vector_get_element(v, i);
        printf("%f\n", val);
    }
}



/* keep only upper triangular matrix part as for symmetric matrix */
void matrix_copy_symmetric(mat *S, mat *M){
    int i,j,n,m;
    m = M->nrows;
    n = M->ncols;
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            if(j>=i){
                matrix_set_element(S,i,j,matrix_get_element(M,i,j));
            }
        }
    }
}


/* initialize diagonal matrix from vector data */
void initialize_diagonal_matrix(mat *D, vec *data){
    int i;
    for(i=0; i<(D->nrows); i++){
        matrix_set_element(D,i,i,data->d[i]);
    }
}


/* invert diagonal matrix */
void invert_diagonal_matrix(mat *Dinv, mat *D){
    int i;
    for(i=0; i<(D->nrows); i++){
        matrix_set_element(Dinv,i,i,1.0/(matrix_get_element(D,i,i)));
    }
}


/* initialize a random matrix */
void initialize_random_matrix(mat *M){
    int i,m,n,N;
    double val;
    m = M->nrows;
    n = M->ncols;
    N = m*n;

    // seed 
    srand(time(NULL));

    // read and set elements
    for(i=0; i<N; i++){
        val = ((double) rand() / (RAND_MAX));
        M->d[i] = val;
    }
}


/* matrix frobenius norm */
double get_matrix_frobenius_norm(mat *M){
    int i;
    double val, normval = 0;
    #pragma omp parallel shared(M,normval) private(i,val) 
    {
    #pragma omp for reduction(+:normval)
    for(i=0; i<((M->nrows)*(M->ncols)); i++){
        val = M->d[i];
        normval += val*val;
    }
    }
    return sqrt(normval);
}



/* C = A*B ; column major */
void matrix_matrix_mult(mat *A, mat *B, mat *C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->nrows, B->ncols, A->ncols, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
    culaDgemm('N', 'N', A->nrows, B->ncols, A->ncols, alpha, A->d, A->nrows, B->d, B->nrows, beta, C->d, C->nrows);
}


/* C = A^T*B ; column major */
void matrix_transpose_matrix_mult(mat *A, mat *B, mat *C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, A->ncols, B->ncols, A->nrows, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
    culaDgemm('T', 'N', A->ncols, B->ncols, A->nrows, alpha, A->d, A->nrows, B->d, B->nrows, beta, C->d, C->nrows);
}


/* C = A*B^T ; column major */
void matrix_matrix_transpose_mult(mat *A, mat *B, mat *C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, A->nrows, B->nrows, A->ncols, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
    culaDgemm('N', 'T', A->nrows, B->nrows, A->ncols, alpha, A->d, A->nrows, B->d, B->nrows, beta, C->d, C->nrows);
}


/* y = M*x */
void matrix_vector_mult(mat *M, vec *x, vec *y){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemv (CblasRowMajor, CblasNoTrans, M->nrows, M->ncols, alpha, M->d, M->ncols, x->d, 1, beta, y->d, 1);
    culaDgemv ('N', M->nrows, M->ncols, alpha, M->d, M->nrows, x->d, 1, beta, y->d, 1);
}


/* y = M^T*x */
void matrix_transpose_vector_mult(mat *M, vec *x, vec *y){
    //gsl_blas_dgemv (CblasNoTrans, 1.0, M, x, 0.0, y);
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemv (CblasRowMajor, CblasTrans, M->nrows, M->ncols, alpha, M->d, M->ncols, x->d, 1, beta, y->d, 1);
    culaDgemv ('T', M->nrows, M->ncols, alpha, M->d, M->nrows, x->d, 1, beta, y->d, 1);
}




/* subtract b from a and save result in a  */
void vector_sub(vec *a, vec *b){
    int i;
    //#pragma omp parallel for
    #pragma omp parallel shared(a,b) private(i) 
    {
    #pragma omp for 
    for(i=0; i<(a->nrows); i++){
        a->d[i] = a->d[i] - b->d[i];
    }
    }
}




/* returns the dot product of two vectors */
double vector_dot_product(vec *u, vec *v){
    int i;
    double dotval = 0;
    for(i=0; i<u->nrows; i++){
        dotval += (u->d[i])*(v->d[i]);
    }
    return dotval;
}


/* subtract B from A and save result in A  */
void matrix_sub(mat *A, mat *B){
    int i;
    //#pragma omp parallel for
    #pragma omp parallel shared(A,B) private(i) 
    {
    #pragma omp for 
    for(i=0; i<((A->nrows)*(A->ncols)); i++){
        A->d[i] = A->d[i] - B->d[i];
    }
    }
}



/*
% project v in direction of u
function p=project_vec(v,u)
p = (dot(v,u)/norm(u)^2)*u;
*/
void project_vector(vec *v, vec *u, vec *p){
    double dot_product_val, vec_norm, scalar_val; 
    dot_product_val = vector_dot_product(v, u);
    vec_norm = vector_get2norm(u);
    scalar_val = dot_product_val/(vec_norm*vec_norm);
    vector_copy(p, u);
    vector_scale(p, scalar_val); 
}


/* build orthonormal basis matrix
Q = Y;
for j=1:k
    vj = Q(:,j);
    for i=1:(j-1)
        vi = Q(:,i);
        vj = vj - project_vec(vj,vi);
    end
    vj = vj/norm(vj);
    Q(:,j) = vj;
end
*/
void build_orthonormal_basis_from_mat(mat *A, mat *Q){
    int m,n,i,j,ind,num_ortos=2;
    double vec_norm;
    vec *vi,*vj,*p;
    m = A->nrows;
    n = A->ncols;
    vi = vector_new(m);
    vj = vector_new(m);
    p = vector_new(m);
    matrix_copy(Q, A);

    for(ind=0; ind<num_ortos; ind++){
        for(j=0; j<n; j++){
            matrix_get_col(Q, j, vj);
            for(i=0; i<j; i++){
                matrix_get_col(Q, i, vi);
                project_vector(vj, vi, p);
                vector_sub(vj, p);
            }
            vec_norm = vector_get2norm(vj);
            vector_scale(vj, 1.0/vec_norm);
            matrix_set_col(Q, j, vj);
        }
    }
}



/* Performs [Q,R] = qr(M,'0') compact QR factorization 
M is mxn ; Q is mxn ; R is min(m,n) x min(m,n) */ 
void compact_QR_factorization(mat *M, mat *Q, mat *R){
    int i,j,m,n,k;
    m = M->nrows; n = M->ncols;
    k = min(m,n);
    mat *R_full = matrix_new(m,n);
    matrix_copy(R_full,M);
    vec *tau = vector_new(m);

    // get R
    //LAPACKE_dgeqrf(CblasRowMajor, m, n, R_full->d, n, tau->d);
    culaDgeqrf(m, n, R_full->d, m, tau->d);
    
    for(i=0; i<k; i++){
        for(j=0; j<k; j++){
            if(j>=i){
                matrix_set_element(R,i,j,matrix_get_element(R_full,i,j));
            }
        }
    }

    // get Q
    matrix_copy(Q,R_full); 
    //LAPACKE_dorgqr(CblasRowMajor, m, n, n, Q->d, n, tau->d);
    culaDorgqr(m, n, n, Q->d, m, tau->d);

    // clean up
    matrix_delete(R_full);
    vector_delete(tau);
}




/* returns Q from [Q,R] = qr(M,'0') compact QR factorization 
M is mxn ; Q is mxn ; R is min(m,n) x min(m,n) */ 
void QR_factorization_getQ(mat *M, mat *Q){
    int i,j,m,n,k;
    m = M->nrows; n = M->ncols;
    k = min(m,n);
    matrix_copy(Q,M);
    vec *tau = vector_new(m);

    culaDgeqrf(m, n, Q->d, m, tau->d);
    culaDorgqr(m, n, n, Q->d, m, tau->d);

    // clean up
    vector_delete(tau);
}



/* compute evals and evecs of symmetric matrix M
   matrix S must be symmetric
*/
void compute_evals_and_evecs_of_symm_matrix(mat *S, vec *evals){
    //LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', S->nrows, S->d, S->nrows, evals->d);
    culaDsyev('V', 'U', S->nrows, S->d, S->ncols, evals->d);
}



/* computes SVD: M = U*S*V^T; note Vt is V transposed */
void singular_value_decomposition(mat *M, mat *U, mat *S, mat *Vt){
    int m,n,k;
    m = M->nrows; n = M->ncols;
    k = min(m,n);
    vec * work = vector_new(k);
    vec * svals = vector_new(k);
    //LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m, n, M->d, n, svals->d, U->d, m, Vt->d, n, work->d );
    culaDgesvd('A', 'A', k, k, M->d, k, svals->d, U->d, k, Vt->d, k);
    initialize_diagonal_matrix(S, svals);
    vector_delete(work);
    vector_delete(svals);
}



/* P = U*S*V^T */
void form_svd_product_matrix(mat *U, mat *S, mat *V, mat *P){
    int k,m,n;
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    m = P->nrows;
    n = P->ncols;
    k = S->nrows;
    mat * SVt = matrix_new(k,n);

    // form SVt = S*V^T
    matrix_matrix_transpose_mult(S,V,SVt);

    // form P = U*S*V^T
    matrix_matrix_mult(U,SVt,P);
}


/* append matrices side by side: C = [A, B] */
void append_matrices_horizontally(mat *A, mat *B, mat *C){
    int i,j;

    for(i=0; i<A->nrows; i++){
        for(j=0; j<A->ncols; j++){
            matrix_set_element(C,i,j,matrix_get_element(A,i,j));
        }
    }

    for(i=0; i<B->nrows; i++){
        for(j=0; j<B->ncols; j++){
            matrix_set_element(C,i,A->ncols + j,matrix_get_element(B,i,j));
        }
    }

    /*
    for(i=0; i<((A->nrows)*(A->ncols)); i++){
        C->d[ind] = A->d[i];
        ind++;
    }

    for(i=0; i<((B->nrows)*(B->ncols)); i++){
        C->d[ind] = B->d[i];
        ind++;
    }*/
}



/* append matrices vertically: C = [A; B] */
void append_matrices_vertically(mat *A, mat *B, mat *C){
    int i,j;

    for(i=0; i<A->nrows; i++){
        for(j=0; j<A->ncols; j++){
            matrix_set_element(C,i,j,matrix_get_element(A,i,j));
        }
    }

    for(i=0; i<B->nrows; i++){
        for(j=0; j<B->ncols; j++){
            matrix_set_element(C,A->nrows+i,j,matrix_get_element(B,i,j));
        }
    }
}





/* calculate percent error between A and B: 100*norm(A - B)/norm(A) */
double get_percent_error_between_two_mats(mat *A, mat *B){
    int m,n;
    double normA, normB, normA_minus_B;
    m = A->nrows;
    n = A->ncols;
    mat *A_minus_B = matrix_new(m,n);
    matrix_copy(A_minus_B, A);
    matrix_sub(A_minus_B, B);
    normA = get_matrix_frobenius_norm(A);
    normB = get_matrix_frobenius_norm(B);
    normA_minus_B = get_matrix_frobenius_norm(A_minus_B);
    return 100.0*normA_minus_B/normA;
}


void checkStatus(culaStatus status)
{
    char buf[256];

    if(!status)
        return;

    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);

    culaShutdown();
    exit(EXIT_FAILURE);
}


void estimate_rank_and_buildQ(mat *M, double frac_of_max_rank, double TOL, mat **Q, int *good_rank){
    int m,n,i,j,ind,maxdim;
    double vec_norm;
    mat *RN,*Y,*Qbig,*Qsmall;
    vec *vi,*vj,*p,*p1;
    m = M->nrows;
    n = M->ncols;
    maxdim = round(min(m,n)*frac_of_max_rank);

    vi = vector_new(m);
    vj = vector_new(m);
    p = vector_new(m);
    p1 = vector_new(m);

    // build random matrix
    printf("form RN..\n");
    RN = matrix_new(n, maxdim);
    initialize_random_matrix(RN);

    // multiply to get matrix of random samples Y
    printf("form Y: %d x %d..\n",m,maxdim);
    Y = matrix_new(m, maxdim);
    matrix_matrix_mult(M, RN, Y);

    // estimate rank k and build Q from Y
    printf("form Qbig..\n");
    Qbig = matrix_new(m, maxdim);

    matrix_copy(Qbig, Y);

    printf("estimate rank with TOL = %f..\n", TOL);
    *good_rank = maxdim;
    //for(ind=0; ind<num_ortos; ind++){
    int forbreak = 0;
    for(j=0; !forbreak && j<maxdim; j++){
        matrix_get_col(Qbig, j, vj);
        for(i=0; i<j; i++){
            matrix_get_col(Qbig, i, vi);
            project_vector(vj, vi, p);
            vector_sub(vj, p);
            if(vector_get2norm(p) < TOL && vector_get2norm(p1) < TOL){
                *good_rank = j;
                forbreak = 1;
                break;
            }
            vector_copy(p1,p);
        }
        vec_norm = vector_get2norm(vj);
        vector_scale(vj, 1.0/vec_norm);
        matrix_set_col(Qbig, j, vj);
    }
    //}

    printf("estimated rank = %d\n", *good_rank);
    Qsmall = matrix_new(m, *good_rank);
    *Q = matrix_new(m, *good_rank);
    matrix_copy_first_columns(Qsmall, Qbig);
    QR_factorization_getQ(Qsmall, *Q);

    matrix_delete(RN);
    matrix_delete(Y);
    matrix_delete(Qsmall);
    matrix_delete(Qbig);
}



void estimate_rank_and_buildQ2(mat *M, int kblock, double TOL, mat **Y, mat **Q, int *good_rank){
    int m,n,i,j,ind,exit_loop = 0;
    double error_norm;
    mat *RN,*Y_new,*Y_big,*QtM,*QQtM;
    vec *vi,*vj,*p,*p1;
    m = M->nrows;
    n = M->ncols;

    // build random matrix
    printf("form RN..\n");
    RN = matrix_new(n,kblock);
    initialize_random_matrix(RN);

    // multiply to get matrix of random samples Y
    printf("form Y: %d x %d..\n",m,kblock);
    *Y = matrix_new(m, kblock);
    matrix_matrix_mult(M, RN, *Y);

    ind = 0;
    while(!exit_loop){
        printf("form Q..\n");
        if(ind > 0){
            matrix_delete(*Q);
        }
        *Q = matrix_new((*Y)->nrows, (*Y)->ncols);
        QR_factorization_getQ(*Y, *Q);

        // compute QtM
        QtM = matrix_new((*Q)->ncols, M->ncols);
        matrix_transpose_matrix_mult(*Q,M,QtM);

        // compute QQtM
        QQtM = matrix_new(M->nrows, M->ncols); 
        matrix_matrix_mult(*Q,QtM,QQtM);

        error_norm = 0.01*get_percent_error_between_two_mats(QQtM, M);

        printf("Y is of size %d x %d and error_norm = %f\n", (*Y)->nrows, (*Y)->ncols, error_norm);
        *good_rank = (*Y)->ncols;
       
        // add more samples if needed
        if(error_norm > TOL){
            Y_new = matrix_new(m, kblock);
            initialize_random_matrix(RN);
            matrix_matrix_mult(M, RN, Y_new);

            Y_big = matrix_new((*Y)->nrows, (*Y)->ncols + Y_new->ncols); 
            append_matrices_horizontally(*Y, Y_new, Y_big);
            matrix_delete(*Y);
            *Y = matrix_new(Y_big->nrows,Y_big->ncols);
            matrix_copy(*Y,Y_big);
            
            matrix_delete(Y_big);
            matrix_delete(Y_new);
            matrix_delete(QtM);
            matrix_delete(QQtM);
            ind++;
        }
        else{
            matrix_delete(RN);
            exit_loop = 1;
        }    
    }
}

