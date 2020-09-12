#ifndef FLAVOUR_MATRIX_H
#define FLAVOUR_MATRIX_H

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"../include/vec.h"


typedef struct FMatrix {
   double complex comp[NFLAVOUR*NFLAVOUR] __attribute__((aligned(DOUBLE_ALIGN)));
} FMatrix;
//
//  the element [i][j] can be obtained by using m(i,j) defined in macro.h
//


// initialize A_{ij}=conj(v_i) v_j - (conj(v_a)v_a/NFLAVOUR) delta_{ij}
inline void init_FMatrix(FMatrix * restrict A, Vec const * const v)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v->comp), DOUBLE_ALIGN);
  #endif

  int i, j;
  double aux;

  for(i=0; i<NFLAVOUR; i++)
     {
     for(j=0; j<NFLAVOUR; j++)
        {
        A->comp[m(i,j)] = conj(v->comp[i])*v->comp[j];
        }
     }

  aux = norm_Vec(v);
  aux = aux*aux/((double) NFLAVOUR);

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[m(i,i)] -= aux;
     }
  }


// A=0
inline void zero_FMatrix(FMatrix * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR*NFLAVOUR; i++)
     {
     A->comp[i]=0.0;
     }
  }


// A=B
inline void equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR*NFLAVOUR; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A+=B
inline void plus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR*NFLAVOUR; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A-=B
inline void minus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR*NFLAVOUR; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A-=B^{dag}
inline void minus_equal_dag_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NFLAVOUR; i++)
     {
     for(j=0; j<NFLAVOUR; j++)
        {
        A->comp[m(i,j)]-=conj(B->comp[m(j,i)]);
        }
     }
  }


// A*=r
inline void times_equal_real_FMatrix(FMatrix * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR*NFLAVOUR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=c
inline void times_equal_complex_FMatrix(FMatrix * restrict A, double complex c)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR*NFLAVOUR; i++)
     {
     A->comp[i]*=c;
     }
  }


// A*=B
inline void times_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex aux[NFLAVOUR] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NFLAVOUR; i++)
     {
     for(j=0; j<NFLAVOUR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NFLAVOUR; j++)
        {
        sum=0.0;
        for(k=0; k<NFLAVOUR; k++)
           {
           sum+=aux[k]*(B->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// real part of the trace
inline double retr_FMatrix(FMatrix const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double tr;

  tr=0.0;
  for(i=0; i<NFLAVOUR; i++)
     {
     tr+=creal(A->comp[m(i,i)]);
     }

  return tr;
  }


// l2 norm of the matrix
inline double norm_FMatrix(FMatrix const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris=0.0;

  for(i=0; i<NFLAVOUR*NFLAVOUR; i++)
     {
     ris+=(cabs(A->comp[i])*cabs(A->comp[i]));
     }

  return sqrt(ris);
  }



#endif
