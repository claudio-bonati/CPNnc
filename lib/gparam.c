#ifndef GPARAM_C
#define GPARAM_C

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/endianness.h"
#include"../include/gparam.h"
#include"../include/macro.h"


// remove from input file white/empty lines and comments
// comments start with the carachter #
void remove_white_line_and_comments(FILE *input)
  {
  int temp_i;

  temp_i=getc(input);
  if(temp_i=='\n' || temp_i==' ' || temp_i=='\043') // scan for white lines and comments
    {
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\n' || temp_i==' ') // white line
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i=='\n' || temp_i==' ');
      }
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\043')  // comment, \043 = ascii oct for #
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i!='\n');
      }
    else
      {
      ungetc(temp_i, input);
      }

    remove_white_line_and_comments(input);
    }
  else
    {
    ungetc(temp_i, input);
    }
  }


void readinput(char *in_file, GParam *param)
    {
    FILE *input;
    char str[STD_STRING_LENGTH], temp_str[STD_STRING_LENGTH];
    double temp_d;
    int temp_i, i;
    int err, end=1;
    unsigned int temp_ui;

    param->d_K=0.0;
    param->d_phmass=0.0;

    input=fopen(in_file, "r");  // open the input file
    if(input==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    else
      {
      while(end==1)   // slide the file
           {
           remove_white_line_and_comments(input);

           err=fscanf(input, "%s", str);
           if(err!=1)
             {
             fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
             printf("err=%d\n", err);
             exit(EXIT_FAILURE);
             }

           if(strncmp(str, "size", 4)==0)
             {
             for(i=0; i<STDIM; i++)
                {
                err=fscanf(input, "%d", &temp_i);
                if(err!=1)
                  {
                  fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                  exit(EXIT_FAILURE);
                  }
                param->d_size[i]=temp_i;
                }
             }

           else if(strncmp(str, "J", 1)==0)
                  { 
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_J=temp_d;
                  }
           else if(strncmp(str, "K", 1)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_K=temp_d;
                  }
           else if(strncmp(str, "phmass", 6)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_phmass=temp_d;
                  }
            else if(strncmp(str, "gaugefixpar", 11)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_gaugefixpar=temp_d;
                  }

           else if(strncmp(str, "sample", 6)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_sample=temp_i;
                  }
           else if(strncmp(str, "thermal", 7)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_thermal=temp_i;
                  }
           else if(strncmp(str, "overrelax", 9)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_overrelax=temp_i;
                  }
           else if(strncmp(str, "measevery", 9)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_measevery=temp_i;
                  }
          
           else if(strncmp(str, "start", 5)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_start=temp_i;
                  }
           else if(strncmp(str, "saveconf_back_every", 19)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_saveconf_back_every=temp_i;
                  }

           else if(strncmp(str, "epsilon_metro_link", 18)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_epsilon_metro_link=temp_d;
                  }
           else if(strncmp(str, "epsilon_metro_site", 18)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_epsilon_metro_site=temp_d;
                  }

           else if(strncmp(str, "conf_file", 9)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_conf_file, temp_str);
                  }
           else if(strncmp(str, "data_file", 9)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_data_file, temp_str);
                  }
           else if(strncmp(str, "log_file", 8)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_log_file, temp_str);
                  }

           else if(strncmp(str, "randseed", 8)==0)
                  { 
                  err=fscanf(input, "%u", &temp_ui);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_randseed=temp_ui;
                  }

           else
             {
             fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, in_file, __FILE__, __LINE__);
             exit(EXIT_FAILURE);
             }

           remove_white_line_and_comments(input);

           // check if the read line is the last one
           temp_i=getc(input);
           if(temp_i==EOF)
             {
             end=0;
             }
           else
             {
             ungetc(temp_i, input);
             }
           }

      fclose(input);

      err=0;
      for(i=0; i<STDIM; i++)
         {
         if(param->d_size[i]==1)
           {
           err=1;
           }
         }
      if(err==1)
        {
        fprintf(stderr, "Error: all sizes has to be larger than 1: the totally reduced case is not implemented! (%s, %d)\n", __FILE__, __LINE__);
        }

      init_derived_constants(param);
      }
    }


void init_derived_constants(GParam *param)
  {
  int i;

  // derived constants
  param->d_volume=1;
  for(i=0; i<STDIM; i++)
     {
     (param->d_volume)*=(param->d_size[i]);
     }

  param->d_inv_vol=1.0/((double) param->d_volume);
  }


// initialize data file
void init_data_file(FILE **dataf, GParam const * const param)
  {
  int i;

  if(param->d_start==2)
    {
    *dataf=fopen(param->d_data_file, "r");
    if(*dataf!=NULL) // file exists
      {
      fclose(*dataf);
      *dataf=fopen(param->d_data_file, "a");
      }
    else
      {
      *dataf=fopen(param->d_data_file, "w");
      fprintf(*dataf, "%d ", STDIM);
      for(i=0; i<STDIM; i++)
         {
         fprintf(*dataf, "%d ", param->d_size[i]);
         }
      fprintf(*dataf, "\n");
      }
    }
  else
    {
    *dataf=fopen(param->d_data_file, "w");
    fprintf(*dataf, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(*dataf, "%d ", param->d_size[i]);
       }
    fprintf(*dataf, "\n");
    }
  fflush(*dataf);
  }


// print simulation parameters
void print_parameters(GParam const * const param, time_t time_start, time_t time_end, double acc_site, double acc_link)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-------------------------------+\n");
    fprintf(fp, "| Simulation details for cpn_nc |\n");
    fprintf(fp, "+-------------------------------+\n\n");

    #ifdef DEBUG
     fprintf(fp, "DEBUG mode\n");
    #endif

    #ifdef CSTAR_BC
     fprintf(fp, "C^* BOUNDARY CONDITIONS\n");
    #else
     fprintf(fp, "PERIODIC BOUNDARY CONDITIONS\n");
    #endif

    #ifdef LINKS_FIXED_TO_ONE
     fprintf(fp, "LINKS FIXED TO 1\n");
    #endif

    #ifdef HARD_TEMPORAL_GAUGE
     fprintf(fp, "HARD_TEMPORAL GAUGE\n");
    #endif

    #ifdef SOFT_TEMPORAL_GAUGE
     fprintf(fp, "SOFT_TEMPORAL GAUGE\n");
    #endif

    #ifdef SOFT_LORENZ_GAUGE
     fprintf(fp, "SOFT_LORENZ GAUGE\n");
    #endif

    fprintf(fp,"\n");

    fprintf(fp, "number of flavours: %d\n", NFLAVOUR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "J: %.10lf\n", param->d_J);
    fprintf(fp, "K: %.10lf\n", param->d_K);
    fprintf(fp, "phmass: %.10lf\n", param->d_phmass);
    #ifdef GAUGE_FIX
      fprintf(fp, "gaugefixpar: %.10lf\n", param->d_gaugefixpar);
    #endif

    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "epsilon_metro_site: %.10lf\n", param->d_epsilon_metro_site);
    fprintf(fp, "metropolis acceptance site: %.10lf\n", acc_site);

    fprintf(fp, "epsilon_metro_link: %.10lf\n", param->d_epsilon_metro_link);
    fprintf(fp, "metropolis acceptance link: %.10lf\n", acc_link);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


#endif

