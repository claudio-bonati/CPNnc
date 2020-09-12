#ifndef FLAVOUR_MATRIX_C
#define FLAVOUR_MATRIX_C

#include<complex.h>

#include"./include/flavour_matrix.h"
#include"./include/macro.h"
#include"./include/vec.h"

// initialize A_{ij}=conj(v_i) v_j - (conj(v_a)v_a/NFLAVOUR) delta_{ij}
void init_FMatrix(FMatrix * restrict A, Vec const * const v);


// A=0
void zero_FMatrix(FMatrix * restrict A);


// A=B
void equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A+=B
void plus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A-=B
void minus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A-=B^{dag}
void minus_equal_dag_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A*=r
void times_equal_real_FMatrix(FMatrix * restrict A, double r);

// A*=c
void times_equal_complex_FMatrix(FMatrix * restrict A, double complex c);

// A*=B
void times_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// real part of the trace
double retr_FMatrix(FMatrix const * const restrict A);


// l2 norm of the matrix
double norm_FMatrix(FMatrix const * const restrict A);





#endif
