#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "heads.h"
#include "../COMBINATORIAL_STRATEGIES/partitioning/list.h"
#include "../COMBINATORIAL_STRATEGIES/partitioning/fib_heap.h"

#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))

/*----------------------------------------------------------------------------
 * COMMON FUNCTIONS PROTOTYPE 
 *--------------------------------------------------------------------------*/
int      daxpy                           (int n, double  a, double *x, double *y);
int      dcopy                           (int n, double *x, double *y);
double   ddot                            (int n, double *x, double *y);
int      dscal                           (int n, double  a, double *x, double *y);
int      ddiff                           (int n, double *a, double *x, double *y);
int      izero                           (int n, int *v);
int      iszero                          (double x);
char*    get_EXTENSION                   (char *fspec);
double   get_TIME                        ();
void     READ_INPUT                      (MAT* A, double** b, char* entry, int* msgs);
void 	 READ_RHS 			 (double* b, int n, FILE* f);
void     WRITE_SOLUTION                  (double* x, char* entry, char* output, int n, int* P, int* K, double* dw, int RANK, int* msgs);

/*----------------------------------------------------------------------------
 * MATRIX HEADER FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
void     MATRIX_readCSR                  (MAT* A, FILE* file);
void     MATRIX_CSRto_FULL               (MAT* A, double* a);
double   MATRIX_aij                      (MAT* A, int i, int j);
long int MATRIX_bandwidth                (MAT* A);
void     MATRIX_clean                    (MAT* A);
void     MATRIX_matvec                   (MAT* A, double* x, double* b);
void     MATRIX_forward                  (MAT* L, double* b, double* y);
void     MATRIX_backward                 (MAT* U, double* y, double* x);
void     MATRIX_writeCSR                 (MAT* A, double* f, int* s, int nP, int bandwidth, int *msgs);
void     MATRIX_ROW_permutation          (MAT* A, int* p);
void     MATRIX_COL_permutation          (MAT* A, int* p);

/*----------------------------------------------------------------------------
 * COMPARE FUNCTIONS FOR QSORT 
 *--------------------------------------------------------------------------*/
int      COMPARE_array            (const void * a, const void * b);
int      COMPARE_eig              (const void * a, const void * b);
int      COMPARE_degr_ASC         (const void * a, const void * b);
int      COMPARE_dist_degr_DES    (const void * a, const void * b);
int      COMPARE_dist_degr_ASC    (const void * a, const void * b);


/*----------------------------------------------------------------------------
 * LINKED LIST FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
LIST*    LIST_insert_IF_NOT_EXIST        (LIST* L, int x);
LIST*    LIST_insert_ORDERED             (LIST* L, int i, double elem);
LIST*    LIST_remove                     (LIST* L, int x);
LIST*    LIST_remove_first               (LIST* L);
void     LIST_print                      (LIST* L);
int      LIST_first                      (LIST* L);
void     LIST_destroy                    (LIST* L);

/*----------------------------------------------------------------------------
 * SOLVERS FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
double   GMRES_algorithm                 (MAT* A, double *x, double *b, int n, int k, int *iter, double eps);
void     GMRES                           (MAT* A, double *x, double *b, double tol, int restart, int maxiter);
double   PGMRES_algorithm                (MAT* A, MAT* L, MAT* U, double *x, double *b, int n, int k, int *iter, double eps);
void     PGMRES                          (MAT* A, MAT* L, MAT* U, double *x, double *b, double tol, int restart, int maxiter);

/*----------------------------------------------------------------------------
 * REORDERING FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
void     REORDERING_SPECTRAL             (MAT* A, double* b, int** P, int* msgs);
void     REORDERING_SPECTRAL_WGT         (MAT* A, double* b, int** P, int* msgs);
void     REORDERING_MATCHING             (MAT* A, double* b, int** Q, double** dw, int* msgs);
void     REORDERING_DE_MIN               (MAT *A, double* b, int **K, int *s, int nP, int k, int *msgs);
int      REORDERING_DE_MIN_nnzout        (MAT *A, int *s, int nP, int k);
void     PARTITIONING_MIN_MAX            (MAT* A, int nP, int** s, int** L, int* msgs); 
void     PARTITIONING_printFULL          (MAT *A, int *s, int nP, int r);
int      PARTITIONING_nnzout             (MAT *A, int *s, int nP);

/*----------------------------------------------------------------------------
 * PARALLEL FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
double   PARALLEL_ddot                   (int n, double *x, double *y);
void     PARALLEL_MATRIX_readCSR         (MAT* A_bef, MAT* A, MAT* A_aft, MAT* C, MAT* B, double** f, int bandwidth, int RANK, int NPROC, int *msgs);
void     PARALLEL_MATRIX_matvec          (MAT* A_bef, MAT* A, MAT* A_aft, double* x, double* r, int RANK, int NPROC);
void     PARALLEL_e_MATRIX_matvec        (MAT *A_bef, MAT* A, MAT* A_aft, double* x, double* r, int RANK, int NPROC);
double   PARALLEL_GMRES_algorithm        (MAT* A_bef, MAT* A, MAT* A_aft, double*  x, double*  b, int n, int k, int* iter, double eps, int RANK, int NPROC);
void     PARALLEL_GMRES                  (MAT* A_bef, MAT* A, MAT* A_aft, double*  x, double*  b, double tol, int restart, int maxiter, int RANK, int NPROC, int* msgs);
double   PARALLEL_PGMRES_algorithm       (MAT* A_bef, MAT* A, MAT* A_aft, MAT* C, MAT* B, MAT* S, void* A_pt[64], int A_iparm[64], double A_dparm[64], void* S_pt[64], int S_iparm[64], double S_dparm[64], double *x, double *b, int n, int k, int *iter, double eps, int RANK, int NPROC);
void     PARALLEL_PGMRES                 (MAT* A_bef, MAT* A, MAT* A_aft, MAT* C, MAT* B, double*  x, double*  b, double tol, int restart, int maxiter, int RANK, int NPROC, int* msgs);
void     PARALLEL_SPIKE_algorithm        (MAT *C, MAT* A, MAT* B, MAT* S, void* A_pt[64], int A_iparm[64], double A_dparm[64], void* S_pt[64], int S_iparm[64], double S_dparm[64], double* f, double* y, int RANK, int NPROC);
double*  PARALLEL_assemble_solution      (double* x, int n, int N, int RANK, int NPROC);
void     PARALLEL_REDUCED_system         (double* V, double* W, int k, MAT* S, int RANK);

/*----------------------------------------------------------------------------
 * PARDISO
 *--------------------------------------------------------------------------*/
void     pardisoinit                     (void*,int*,int*,int*,double*,int*);
void     pardiso                         (void*,int*,int*,int*,int*,int*,double*,int*,int*,int*,int*,int*,int*,double*,double*,int*,double*);
void     pardiso_chkmatrix               (int*,int*,double*,int*,int*,int*);
void     pardiso_chkvec                  (int*,int*,double*,int*);
void     pardiso_printstats              (int*,int*,double*,int*,int*,int*,double*,int*);
void     read_matrix                     (double**,int**,int**,double**,int*,int*,int*,FILE*);
void     PARDISO_bottom_tips             (MAT* A, double* B, double* V, int nrhs);
void     PARDISO_top_tips                (MAT* A, double* C, double* W, int nrhs);
void     PARDISO_numerical_factorization (MAT* A, void* pt[64], int iparm[64], double dparm[64], int nrhs);
void     PARDISO_back_substitution       (MAT* A, void* pt[64], int iparm[64], double dparm[64], double* f, double* g, int nrhs);
void     PARDISO_release_memory          (MAT* A, void* pt[64], int iparm[64], double dparm[64], int nrhs);

