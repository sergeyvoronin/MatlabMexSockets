#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_vsl.h"


typedef struct {
    int nrows, ncols;
    double * d;
} mat;


typedef struct {
    int nrows;
    double * d;
} vec;



/* initialize new matrix and set all entries to zero */
mat * matrix_new(int nrows, int ncols);

/* initialize new vector and set all entries to zero */
vec * vector_new(int nrows);

void matrix_delete(mat *M);

void vector_delete(vec *v);


/* set element in column major format */
void matrix_set_element(mat *M, int row_num, int col_num, double val);

/* get element in column major format */
double matrix_get_element(mat *M, int row_num, int col_num);


/* set vector element */
void vector_set_element(vec *v, int row_num, double val);


/* get vector element */
double vector_get_element(vec *v, int row_num);


void vector_set_data(vec *v, double *data);


/* scale vector by a constant */
void vector_scale(vec *v, double scalar);


/* scale matrix by a constant */
void matrix_scale(mat *M, double scalar);


/* compute euclidean norm of vector */
double vector_get2norm(vec *v);


/* copy contents of vec s to d  */
void vector_copy(vec *d, vec *s);


/* copy contents of mat S to D  */
void matrix_copy(mat *D, mat *S);


/* hard threshold matrix entries  */
void matrix_hard_threshold(mat *M, double TOL);


/* build transpose of matrix : Mt = M^T */
void matrix_build_transpose(mat *Mt, mat *M);


/* copy contents of mat S to D */
void matrix_copy_first_columns(mat *D, mat *S, int num_columns);


/* subtract b from a and save result in a  */
void vector_sub(vec *a, vec *b);


/* subtract B from A and save result in A  */
void matrix_sub(mat *A, mat *B);


/* matrix frobenius norm */
double get_matrix_frobenius_norm(mat *M);


/* matrix max abs val */
double get_matrix_max_abs_element(mat *M);

/* print out matrix */
void matrix_print(mat * M);

/* print out vector */
void vector_print(vec * v);


/* C = A*B ; column major */
void matrix_matrix_mult(mat *A, mat *B, mat *C);


/* C = A^T*B ; column major */
void matrix_transpose_matrix_mult(mat *A, mat *B, mat *C);


/* C = A*B^T ; column major */
void matrix_matrix_transpose_mult(mat *A, mat *B, mat *C);


/* y = M*x ; column major */
void matrix_vector_mult(mat *M, vec *x, vec *y);


/* y = M^T*x ; column major */
void matrix_transpose_vector_mult(mat *M, vec *x, vec *y);


/* set column of matrix to vector */
void matrix_set_col(mat *M, int j, vec *column_vec);


/* extract column of a matrix into a vector */
void matrix_get_col(mat *M, int j, vec *column_vec);


/* extract row i of a matrix into a vector */
void matrix_get_row(mat *M, int i, vec *row_vec);


/* put vector row_vec as row i of a matrix */
void matrix_set_row(mat *M, int i, vec *row_vec);



/* keep only upper triangular matrix part of matrix in S from M */
void matrix_copy_symmetric(mat *S, mat *M);



/* initialize diagonal matrix from vector data */
void initialize_diagonal_matrix(mat *D, vec *data);



/* initialize identity */
void initialize_identity_matrix(mat *D);



/* invert diagonal matrix */
void invert_diagonal_matrix(mat *Dinv, mat *D);



/* returns the dot product of two vectors */
double vector_dot_product(vec *u, vec *v);



/*
% project v in direction of u
function p=project_vec(v,u)
p = (dot(v,u)/norm(u)^2)*u;
*/
void project_vector(vec *v, vec *u, vec *p);


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
void build_orthonormal_basis_from_mat(mat *A, mat *Q);



void fill_matrix_from_first_rows(mat *M, int k, mat *M_k);


void fill_matrix_from_first_columns(mat *M, int k, mat *M_k);


void fill_matrix_from_last_columns(mat *M, int k, mat *M_k);


void fill_matrix_from_column_list(mat *M, vec *I, mat *M_k);


void fill_matrix_from_row_list(mat *M, vec *I, mat *M_k);


/* calculate percent error between A and B: 100*norm(A - B)/norm(A) */
double get_percent_error_between_two_mats(mat *A, mat *B);


double get_matrix_column_norm_squared(mat *M, int colnum);


double matrix_getmaxcolnorm(mat *M);


void compute_matrix_column_norms(mat *M, vec *column_norms);