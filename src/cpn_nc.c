#ifndef CPN_NC_C
#define CPN_NC_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char *in_file)
    {
    Conf GC;
    Geometry geo;
    GParam param;

    long count;
    FILE *datafilep;

    time_t time1, time2;
    double acc_link, acc_site;
    double acc_link_local, acc_site_local;

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize geometry
    init_indexing_lex();
    init_geometry(&geo, &param);

    // initialize configuration
    init_conf(&GC, &param);

    // acceptance
    acc_link=0.0;
    acc_site=0.0;

    // montecarlo
    time(&time1);
    // count starts from 1 to avoid problems using %
    for(count=1; count < param.d_sample + 1; count++)
       {
       update(&GC, &geo, &param, &acc_site_local, &acc_link_local);

       if(count>param.d_thermal)
         {
         acc_site+=acc_site_local;
         acc_link+=acc_link_local;
         }

       if(count<param.d_thermal)
         {
         if(acc_site_local>0.33)
           {
           param.d_epsilon_metro_site*=1.1;

           if(param.d_epsilon_metro_site>3.0)
             {
             param.d_epsilon_metro_site=3.0;
             }
           }
         else
           {
           param.d_epsilon_metro_site*=0.9;
           }

         if(acc_link_local>0.33)
           {
           param.d_epsilon_metro_link*=1.1;
           }
         else
           {
           param.d_epsilon_metro_link*=0.9;
           }
         }

       if(count % param.d_measevery ==0 && count >= param.d_thermal)
         {
         perform_measures(&GC, &param, &geo, datafilep);
         }

       // save configuration for backup
       if(param.d_saveconf_back_every!=0)
         {
         if(count % param.d_saveconf_back_every == 0 )
           {
           // simple
           write_conf_on_file(&GC, &param);

           // backup copy
           write_conf_on_file_back(&GC, &param);
           }
         }
       }
    time(&time2);
    // montecarlo end

    acc_site/=(double)(param.d_sample-param.d_thermal);
    acc_link/=(double)(param.d_sample-param.d_thermal);

    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      write_conf_on_file(&GC, &param);
      }

    // print simulation details
    print_parameters(&param, time1, time2, acc_site, acc_link);

    // free configuration
    free_conf(&GC, &param);

    // free geometry
    free_geometry(&geo, &param);
    }


void print_template_input(void)
  {
  FILE *fp;

  fp=fopen("template_input.in", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "size 4 4 4\n");
    fprintf(fp,"\n");
    fprintf(fp, "J 5.705\n");
    fprintf(fp, "K 2.0\n");
    fprintf(fp, "phmass 0.0\n");
    fprintf(fp,"\n");
    fprintf(fp, "sample    10\n");
    fprintf(fp, "thermal   0\n");
    fprintf(fp, "overrelax 5\n");
    fprintf(fp, "measevery 1\n");
    fprintf(fp,"\n");
    fprintf(fp, "epsilon_metro_site 0.8\n");
    fprintf(fp, "epsilon_metro_link 0.8\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                   0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every     5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp, "\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    #(0=time)\n");
    fclose(fp);
    }
  }


int main (int argc, char **argv)
    {
    char in_file[50];

    if(argc != 2)
      {
      printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tN_f (number of flavours): %d\n", NFLAVOUR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef CSTAR_BC 
        printf("\n\tC^* BOUNDARY CONDITIONS\n");
      #else
        printf("\n\tPERIODIC BOUNDARY CONDITIONS\n");
      #endif

      #ifdef LINKS_FIXED_TO_ONE
         printf("\n\tLINKS FIXED TO 1\n");
      #endif

      #ifdef TEMPORAL_GAUGE
         printf("\n\tTEMPORAL GAUGE\n");
      #endif

      #ifdef LORENZ_GAUGE
         printf("\n\tLORENZ GAUGE\n");
      #endif

      #ifdef COMPACT_MASS
         printf("\n\tCOMPACT MASS\n");
      #endif

      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      print_template_input();


      return EXIT_SUCCESS;
      }
    else
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
        }
      else
        {
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    }

#endif
