#ifndef VEC_C
#define VEC_C

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<string.h>

#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/vec.h"


// A=1
void one_Vec(Vec * restrict A);


// A=0
void zero_Vec(Vec * restrict A);


// A=B
void equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A=conj(B)
void equal_cc_Vec(Vec * restrict A, Vec const * const restrict B);


// A=Re(B)
void repart_Vec(Vec * restrict A, Vec const * const restrict B);


// A=Im(B)
void impart_Vec(Vec * restrict A, Vec const * const restrict B);


// A+=B
void plus_equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A-=B
void minus_equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A*=r
void times_equal_real_Vec(Vec * restrict A, double r);


// A*=c
void times_equal_complex_Vec(Vec * restrict A, double complex c);


// l2 norm
double norm_Vec(Vec const * const restrict A);


// unitarize
void unitarize_Vec(Vec * restrict A);


// random vector (normalized)
void rand_vec_Vec(Vec * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]=2.0*casuale()-1.0 + (2.0*casuale()-1.0)*I;
     }

  unitarize_Vec(A);
  }


// scalar product v_1^{\dag}v_2
double complex scal_prod_Vec(Vec const * const A, Vec const * const B);


// random rotation close to identity
void rand_rot_Vec(Vec * restrict A, Vec const * const restrict B, double epsilon)
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

  int i, j, counter;
  double c0, c1, c2, c3, norm;
  double complex tmp1, tmp2;

  equal_Vec(A, B);

  for(counter=0; counter<NFLAVOUR; counter++)
     {
     i = (int)(casuale()*((double)NFLAVOUR - MIN_VALUE));
     j = i + 1 + (int)(casuale()*((double)(NFLAVOUR-1) - MIN_VALUE));
     j = j % NFLAVOUR;

     tmp1=A->comp[i];
     tmp2=A->comp[j];

     c0=1.0+(2.0*casuale()-1.0)*epsilon;
     c1=(2.0*casuale()-1.0)*epsilon;
     c2=(2.0*casuale()-1.0)*epsilon;
     c3=(2.0*casuale()-1.0)*epsilon;
     norm=sqrt(c0*c0+c1*c1+c2*c2+c3*c3);
     c0/=norm;
     c1/=norm;
     c2/=norm;
     c3/=norm;

     A->comp[i] = (c0+I*c1)*tmp1 -(c2-I*c3)*tmp2;
     A->comp[j] = (c2+I*c3)*tmp1 +(c0-I*c1)*tmp2;
     }

  for(counter=0; counter<NFLAVOUR; counter++)
     {
     A->comp[counter] *= cexp(I*(2.0*casuale()-1.0)*epsilon);
     }
  }



// random rotation close to identity for a single couple of indices
void rand_rot_single_Vec(Vec * restrict A, Vec const * const restrict B, double epsilon)
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
  double c0, c1, c2, c3, norm;
  double complex tmp1, tmp2;

  equal_Vec(A, B);

  i = (int)(casuale()*((double)NFLAVOUR - MIN_VALUE));
  j = i + 1 + (int)(casuale()*((double)(NFLAVOUR-1) - MIN_VALUE));
  j = j % NFLAVOUR;

  tmp1=A->comp[i];
  tmp2=A->comp[j];

  c0=1.0+(2.0*casuale()-1.0)*epsilon;
  c1=(2.0*casuale()-1.0)*epsilon;
  c2=(2.0*casuale()-1.0)*epsilon;
  c3=(2.0*casuale()-1.0)*epsilon;
  norm=sqrt(c0*c0+c1*c1+c2*c2+c3*c3);
  c0/=norm;
  c1/=norm;
  c2/=norm;
  c3/=norm;

  A->comp[i] = (c0+I*c1)*tmp1 -(c2-I*c3)*tmp2;
  A->comp[j] = (c2+I*c3)*tmp1 +(c0-I*c1)*tmp2;

  A->comp[i] *= cexp(I*(2.0*casuale()-1.0)*epsilon);
  A->comp[j] *= cexp(I*(2.0*casuale()-1.0)*epsilon);
  }


// print on screen
int print_on_screen_Vec(Vec const * const A)
  {
  int i, err;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=printf("%.16f %.16f ", creal(A->comp[i]), cimag(A->comp[i]));
     if(err<0)
       {
       fprintf(stderr, "Problem in writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }
  fprintf(stdout, "\n");

  return 0;
  }


// print on file
int print_on_file_Vec(FILE *fp, Vec const * const A)
  {
  int i, err;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fprintf(fp, "%.16f %.16f ", creal(A->comp[i]), cimag(A->comp[i]));
     if(err<0)
       {
       fprintf(stderr, "Problem in writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }
  fprintf(fp, "\n");

  return 0;
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_Vec(FILE *fp, Vec const * const A)
  {
  int i;
  size_t err;
  double re, im;

  for(i=0; i<NFLAVOUR; i++)
     {
     re=creal(A->comp[i]);
     im=cimag(A->comp[i]);

     err=fwrite(&re, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     err=fwrite(&im, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_Vec(FILE *fp, Vec const * const A)
  {
  int i;
  size_t err;
  double re, im;

  for(i=0; i<NFLAVOUR; i++)
     {
     re=creal(A->comp[i]);
     im=cimag(A->comp[i]);

     SwapBytesDouble(&re);
     SwapBytesDouble(&im);

     err=fwrite(&re, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     err=fwrite(&im, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file in bigendian
int print_on_binary_file_bigen_Vec(FILE *fp, Vec const * const A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_Vec(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_Vec(fp, A);
    }

  return err;
  }


// read from file
int read_from_file_Vec(FILE *fp, Vec *A)
  {
  int i, err;
  double re, im;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fscanf(fp, "%lg", &re);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     err=fscanf(fp, "%lg", &im);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     A->comp[i]=re+im*I;
     }

  return 0;
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_Vec(FILE *fp, Vec *A)
  {
  size_t err;
  int i;
  double aux[2], re, im;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fread(&re, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     err=fread(&im, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     aux[0]=re;
     aux[1]=im;

     memcpy((void *)&(A->comp[i]), (void*)aux, sizeof(aux));
     //equivalent to A->comp[i]=re+im*I;
     }

  return 0;
  }


// read from binary file changing endianness
int read_from_binary_file_swap_Vec(FILE *fp, Vec *A)
  {
  int i;
  size_t err;
  double aux[2], re, im;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fread(&re, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     err=fread(&im, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     SwapBytesDouble(&re);
     SwapBytesDouble(&im);
     aux[0]=re;
     aux[1]=im;

     memcpy((void *)&(A->comp[i]), (void*)aux, sizeof(aux));
     //equivalent to A->comp[i]=re+im*I;
     }

  return 0;
  }


// read from binary file written in bigendian
int read_from_binary_file_bigen_Vec(FILE *fp, Vec *A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_Vec(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_Vec(fp, A);
    }

  return err;
  }

#endif

