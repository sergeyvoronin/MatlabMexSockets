/* high level matrix/vector functions usibg Intel MKL for blas */

#include "matrix_vector_functions_intel_mkl.h"


/* initialize new matrix and set all entries to zero */
mat * matrix_new(int nrows, int ncols)
{
    mat *M = malloc(sizeof(mat));
    //M->d = (double*)mkl_calloc(nrows*ncols, sizeof(double), 64);
    M->d = (double*)calloc(nrows*ncols, sizeof(double));
    M->nrows = nrows;
    M->ncols = ncols;
    return M;
}


/* initialize new vector and set all entries to zero */
vec * vector_new(int nrows)
{
    vec *v = malloc(sizeof(vec));
    //v->d = (double*)mkl_calloc(nrows,sizeof(double), 64);
    v->d = (double*)calloc(nrows,sizeof(double));
    v->nrows = nrows;
    return v;
}


void matrix_delete(mat *M)
{
    //mkl_free(M->d);
    free(M->d);
    free(M);
}


void vector_delete(vec *v)
{
    //mkl_free(v->d);
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


/* scale matrix by a constant */
void matrix_scale(mat *M, double scalar){
    int i;
    #pragma omp parallel shared(M,scalar) private(i) 
    {
    #pragma omp for
    for(i=0; i<((M->nrows)*(M->ncols)); i++){
        M->d[i] = scalar*(M->d[i]);
    }
    }
}


/* compute euclidean norm of vector */
double vector_get2norm(vec *v){
    int i;
    double val, normval = 0;
    //#pragma omp parallel for
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



/* hard threshold matrix entries  */
void matrix_hard_threshold(mat *M, double TOL){
    int i;
    #pragma omp parallel shared(M) private(i) 
    {
    #pragma omp for 
    for(i=0; i<((M->nrows)*(M->ncols)); i++){
        if(fabs(M->d[i]) < TOL){
            M->d[i] = 0;
        }
    }
    }
}


/* build transpose of matrix : Mt = M^T */
void matrix_build_transpose(mat *Mt, mat *M){
    int i,j;
    for(i=0; i<(M->nrows); i++){
        for(j=0; j<(M->ncols); j++){
            matrix_set_element(Mt,j,i,matrix_get_element(M,i,j)); 
        }
    }
}




/* copy contents of mat S to D */
void matrix_copy_first_columns(mat *D, mat *S, int num_columns){
    int i,j;
    for(i=0; i<(S->nrows); i++){
        for(j=0; j<num_columns; j++){
            matrix_set_element(D,i,j,matrix_get_element(S,i,j));
        }
    }
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


/* matrix max abs val */
double get_matrix_max_abs_element(mat *M){
    int i;
    double val, max = 0;
    for(i=0; i<((M->nrows)*(M->ncols)); i++){
        val = M->d[i];
        if( abs(val) > max )
            max = val;
    }
    return max;
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


/* C = A*B ; column major */
void matrix_matrix_mult(mat *A, mat *B, mat *C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A->nrows, B->ncols, A->ncols, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A->nrows, B->ncols, A->ncols, alpha, A->d, A->nrows, B->d, B->nrows, beta, C->d, C->nrows);
}


/* C = A^T*B ; column major */
void matrix_transpose_matrix_mult(mat *A, mat *B, mat *C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, A->ncols, B->ncols, A->nrows, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, A->ncols, B->ncols, A->nrows, alpha, A->d, A->nrows, B->d, B->nrows, beta, C->d, C->nrows);
}


/* C = A*B^T ; column major */
void matrix_matrix_transpose_mult(mat *A, mat *B, mat *C){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, A->nrows, B->nrows, A->ncols, alpha, A->d, A->ncols, B->d, B->ncols, beta, C->d, C->ncols);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, A->nrows, B->nrows, A->ncols, alpha, A->d, A->nrows, B->d, B->nrows, beta, C->d, C->nrows);
}


/* y = M*x ; column major */
void matrix_vector_mult(mat *M, vec *x, vec *y){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    cblas_dgemv (CblasColMajor, CblasNoTrans, M->nrows, M->ncols, alpha, M->d, M->nrows, x->d, 1, beta, y->d, 1);
}


/* y = M^T*x ; column major */
void matrix_transpose_vector_mult(mat *M, vec *x, vec *y){
    double alpha, beta;
    alpha = 1.0; beta = 0.0;
    cblas_dgemv (CblasColMajor, CblasTrans, M->nrows, M->ncols, alpha, M->d, M->nrows, x->d, 1, beta, y->d, 1);
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



/* copy only upper triangular matrix part as for symmetric matrix */
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



/* copy only upper triangular matrix part as for symmetric matrix */
void matrix_keep_only_upper_triangular(mat *M){
    int i,j,n,m;
    m = M->nrows;
    n = M->ncols;
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            if(j<i){
                matrix_set_element(M,i,j,0);
            }
        }
    }
}



/* initialize diagonal matrix from vector data */
void initialize_diagonal_matrix(mat *D, vec *data){
    int i;
    #pragma omp parallel shared(D) private(i)
    { 
    #pragma omp parallel for
    for(i=0; i<(D->nrows); i++){
        matrix_set_element(D,i,i,data->d[i]);
    }
    }
}



/* initialize identity */
void initialize_identity_matrix(mat *D){
    int i;
    #pragma omp parallel shared(D) private(i)
    { 
    #pragma omp parallel for
    for(i=0; i<(D->nrows); i++){
        matrix_set_element(D,i,i,1.0);
    }
    }
}



/* invert diagonal matrix */
void invert_diagonal_matrix(mat *Dinv, mat *D){
    int i;
    #pragma omp parallel shared(D,Dinv) private(i)
    {
    #pragma omp parallel for
    for(i=0; i<(D->nrows); i++){
        matrix_set_element(Dinv,i,i,1.0/(matrix_get_element(D,i,i)));
    }
    }
}



/* returns the dot product of two vectors */
double vector_dot_product(vec *u, vec *v){
    int i;
    double dotval = 0;
    #pragma omp parallel shared(u,v,dotval) private(i) 
    {
    #pragma omp for reduction(+:dotval)
    for(i=0; i<u->nrows; i++){
        dotval += (u->d[i])*(v->d[i]);
    }
    }
    return dotval;
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




void fill_matrix_from_first_rows(mat *M, int k, mat *M_k){
    int i;
    vec *row_vec;
    //#pragma omp parallel shared(M,M_k,k) private(i,row_vec) 
    {
    //#pragma omp for
    for(i=0; i<k; i++){
        row_vec = vector_new(M->ncols);
        matrix_get_row(M,i,row_vec);
        matrix_set_row(M_k,i,row_vec);
        vector_delete(row_vec);
    }
    }
}


void fill_matrix_from_first_columns(mat *M, int k, mat *M_k){
    int i;
    vec *col_vec;
    //#pragma omp parallel shared(M,M_k,k) private(i,col_vec) 
    {
    //#pragma omp for
    for(i=0; i<k; i++){
        col_vec = vector_new(M->nrows);
        matrix_get_col(M,i,col_vec);
        matrix_set_col(M_k,i,col_vec);
        vector_delete(col_vec);
    }
    }
}


void fill_matrix_from_last_columns(mat *M, int k, mat *M_k){
    int i,ind;
    vec *col_vec;
    ind = 0;
    for(i=k; i<M->ncols; i++){
        col_vec = vector_new(M->nrows);
        matrix_get_col(M,i,col_vec);
        matrix_set_col(M_k,ind,col_vec);
        vector_delete(col_vec);
        ind++;
    }
}


void fill_matrix_from_column_list(mat *M, vec *I, mat *M_k){
    int i,col_num;
    vec *col_vec;
    for(i=0; i<(M_k->ncols); i++){
        col_num = vector_get_element(I,i);
        col_vec = vector_new(M->nrows);
        matrix_get_col(M,col_num,col_vec);
        matrix_set_col(M_k,i,col_vec);
        vector_delete(col_vec);
    }
}



void fill_matrix_from_row_list(mat *M, vec *I, mat *M_k){
    int i,row_num;
    vec *row_vec;
    for(i=0; i<(M_k->nrows); i++){
        row_num = vector_get_element(I,i);
        row_vec = vector_new(M->ncols);
        matrix_get_row(M,row_num,row_vec);
        matrix_set_row(M_k,i,row_vec);
        vector_delete(row_vec);
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


double get_matrix_column_norm_squared(mat *M, int colnum){
    int i, m, n;
    double val,colnorm;
    m = M->nrows;
    n = M->ncols;
    colnorm = 0;
    for(i=0; i<m; i++){
        val = matrix_get_element(M,i,colnum);
        colnorm += val*val;
    }
    return colnorm;
}


double matrix_getmaxcolnorm(mat *M){
    int i,m,n;
    vec *col_vec;
    double vecnorm, maxnorm;
    m = M->nrows; n = M->ncols;
    col_vec = vector_new(m);

    maxnorm = 0;    
    #pragma omp parallel for
    for(i=0; i<n; i++){
        matrix_get_col(M,i,col_vec);
        vecnorm = vector_get2norm(col_vec);
        #pragma omp critical
        if(vecnorm > maxnorm){
            maxnorm = vecnorm;
        }
    }

    vector_delete(col_vec);
    return maxnorm;
}


void compute_matrix_column_norms(mat *M, vec *column_norms){
    int i,m,n;
    m = M->nrows; n = M->ncols;
    for(i=0; i<n; i++){
        vector_set_element(column_norms,i, get_matrix_column_norm_squared(M,i)); 
    }
}


