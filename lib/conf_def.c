#ifndef CONF_DEF_C
#define CONF_DEF_C

#include"../include/macro.h"

#include<openssl/md5.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/endianness.h"
#include"../include/conf.h"
#include"../include/endianness.h"
#include"../include/flavour_matrix.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"


void init_conf(Conf *GC,
               GParam const * const param)
  {
  long r, j;
  int err;

  // allocate the lattice
  err=posix_memalign((void**) &(GC->theta), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC->theta[r]), (size_t) DOUBLE_ALIGN, (size_t )STDIM * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  err=posix_memalign((void**) &(GC->phi), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(Vec));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  err=posix_memalign((void**) &(GC->Qh), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(FMatrix));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize lattice
  if(param->d_start==0) // ordered start
    {
    Vec v1, v2;

    GC->update_index=0;

    one_Vec(&v1);

    for(r=0; r<(param->d_volume); r++)
       {
       rand_rot_Vec(&v2, &v1, 0.05);
       equal_Vec(&(GC->phi[r]), &v2);

       for(j=0; j<STDIM; j++)
          {
          GC->theta[r][j]=0.0+0.05*(2.0*casuale()-1.0);
          }
       }
    }

  if(param->d_start==1)  // random start
    {
    Vec v1;

    GC->update_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       rand_vec_Vec(&v1);
       equal_Vec(&(GC->phi[r]), &v1);

       for(j=0; j<STDIM; j++)
          {
          GC->theta[r][j]=2.0*casuale()-1.0;
          }
       }
    }

  if(param->d_start==2) // initialize from stored conf
    {
    read_conf(GC, param);
    }

  #ifdef LINKS_FIXED_TO_ONE
    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<STDIM; j++)
          {
          GC->theta[r][j]=0.0;
          }
       }
  #endif

  #ifdef HARD_TEMPORAL_GAUGE
    for(r=0; r<(param->d_volume); r++)
       {
       GC->theta[r][0]=0.0;
       }
  #endif
  }


void init_conf_from_conf(Conf *GCnew,
                         Conf const * const GCold,
                         GParam const * const param)
  {
  long r, j;
  int err;

  // allocate the lattice
  err=posix_memalign((void**) &(GCnew->theta), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GCnew->theta[r]), (size_t) DOUBLE_ALIGN, (size_t )STDIM * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  err=posix_memalign((void**) &(GCnew->phi), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(Vec));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  err=posix_memalign((void**) &(GCnew->Qh), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(FMatrix));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize lattice
  GCnew->update_index=GCold->update_index;

  for(r=0; r<(param->d_volume); r++)
     {
     equal_Vec(&(GCnew->phi[r]), &(GCold->phi[r]));

     for(j=0; j<STDIM; j++)
        {
        GCnew->theta[r][j]=GCold->theta[r][j];
        }
     }
  }


void read_conf(Conf *GC, GParam const * const param)
  {
  FILE *fp;
  int i, dimension, tmp_i;
  int err, mu;
  long r;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_conf_file, "r"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else // read the txt header of the configuration
    {
    err=fscanf(fp, "%d", &dimension);
    if(err!=1)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(dimension != STDIM)
      {
      fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n",
              dimension, STDIM);
      exit(EXIT_FAILURE);
      }

    for(i=0; i<STDIM; i++)
       {
       err=fscanf(fp, "%d", &tmp_i);
       if(err!=1)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       if(tmp_i != param->d_size[i])
         {
         fprintf(stderr, "The size of the configuration lattice does not coincide with the one of the global parameter\n");
         exit(EXIT_FAILURE);
         }
       }

    err=fscanf(fp, "%ld %s\n", &(GC->update_index), md5sum_old);
    if(err!=2)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    fclose(fp);
    }

  fp=fopen(param->d_conf_file, "rb"); // open the configuration file in binary
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header
    err=0;
    while(err!='\n')
         {
         err=fgetc(fp);
         }

    for(r=0; r<param->d_volume; r++)
       {
       err=read_from_binary_file_bigen_Vec(fp, &(GC->phi[r]));
       if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }

       for(mu=0; mu<STDIM; mu++)
          {
          err=read_from_binary_file_bigen_double(fp, &(GC->theta[r][mu]));
          if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    fclose(fp);

    // compute the new md5sum and check for consistency
    compute_md5sum_conf(md5sum_new, GC, param);
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
      {
      fprintf(stderr, "The computed md5sum %s does not match the stored %s for the file %s (%s, %d)\n", md5sum_new, md5sum_old, param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  }


void free_conf(Conf *GC, GParam const * const param)
  {
  long i;

  for(i=0; i<(param->d_volume); i++)
     {
     free(GC->theta[i]);
     }
  free(GC->theta);

  free(GC->phi);

  free(GC->Qh);
  }


// save a configuration in ILDG-like format
void write_conf_on_file_with_name(Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile)
  {
  long r;
  int i, mu, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_conf(md5sum, GC, param);

  fp=fopen(namefile, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(fp, "%d ", param->d_size[i]);
       }
    fprintf(fp, "%ld %s\n", GC->update_index, md5sum);
    }
  fclose(fp);

  fp=fopen(namefile, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(r=0; r<param->d_volume; r++)
       {
       err=print_on_binary_file_bigen_Vec(fp, &(GC->phi[r]) );
       if(err!=0)
         {
         fprintf(stderr, "Error in writing the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }

       for(mu=0; mu<STDIM; mu++)
          {
          err=print_on_binary_file_bigen_double(fp, GC->theta[r][mu] );
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    fclose(fp);
    }
  }


void write_conf_on_file(Conf const * const GC, GParam const * const param)
  {
  write_conf_on_file_with_name(GC, param, param->d_conf_file);
  }


void write_conf_on_file_back(Conf const * const GC, GParam const * const param)
  {
  char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
  static int counter=0;

  strcpy(name, param->d_conf_file);
  if(counter==0)
    {
    sprintf(aux, "_back0");
    }
  else
    {
    sprintf(aux, "_back1");
    }
  strcat(name, aux);

  write_conf_on_file_with_name(GC, param, name);

  counter=1-counter;
  }


// compute the md5sum of the configuration and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_conf(char *res, Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long r;
  double a, b;
  int mu, k;

  MD5_Init(&mdContext);
  for(r=0; r<param->d_volume; r++)
     {
     for(k=0; k<NFLAVOUR; k++)
        {
        a=creal((GC->phi[r]).comp[k]);
        b=cimag((GC->phi[r]).comp[k]);

        if(endian()==0)
          {
          SwapBytesDouble(&a);
          SwapBytesDouble(&b);
          }
        MD5_Update(&mdContext, &a, sizeof(double));
        MD5_Update(&mdContext, &b, sizeof(double));
        }

     for(mu=0; mu<STDIM; mu++)
        {
        a=GC->theta[r][mu];
        if(endian()==0)
          {
          SwapBytesDouble(&a);
          }
        MD5_Update(&mdContext, &a, sizeof(double));
        }
     }
  MD5_Final(c, &mdContext);

  for(k = 0; k < MD5_DIGEST_LENGTH; k++)
     {
     sprintf(&(res[2*k]), "%02x", c[k]);
     }
  }



#endif
