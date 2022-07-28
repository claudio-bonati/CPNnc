#ifndef CONF_H
#define CONF_H

#include"macro.h"

#include<complex.h>
#include<openssl/md5.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"gparam.h"
#include"geometry.h"
#include"vec.h"


typedef struct Conf {
  long update_index;

  double **theta;    // [volume][STDIM]
  Vec *phi;          // [volume]

  FMatrix *Qh;          // [volume]
  } Conf;


// in conf_def.c
void init_conf(Conf *GC,
               GParam const * const param);
void init_conf_from_conf(Conf *GCnew,
                         Conf const * const GCold,
                         GParam const * const param);
void read_conf(Conf *GC,
               GParam const * const param);
void free_conf(Conf *GC,
               GParam const * const param);
void write_conf_on_file_with_name(Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile);
void write_conf_on_file(Conf const * const GC,
                        GParam const * const param);
void write_conf_on_file_back(Conf const * const GC,
                             GParam const * const param);
void compute_md5sum_conf(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                         Conf const * const GC,
                         GParam const * const param);


// in conf_update.c
void calcstaples_for_phi(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         Vec *staple);
int metropolis_for_phi(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r);
void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            long r);
double plaqstaples_for_link(Conf *GC,
                            Geometry const * const geo,
                            long r,
                            int i);
double lorenzstaples_for_link(Conf *GC,
                              Geometry const * const geo,
                              long r,
                              int i);
int metropolis_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i);
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc_site,
            double *acc_link);

// in conf_meas.c
double plaquette_single(Conf const * const GC,
                        Geometry const * const geo,
                        long r,
                        int i,
                        int j);
double plaquettesq_single(Conf const * const GC,
                          Geometry const * const geo,
                          long r,
                          int i,
                          int j);
double plaquettesq(Conf const * const GC,
                   Geometry const * const geo,
                   GParam const * const param);
double higgs_interaction(Conf const * const GC,
                         Geometry const * const geo,
                         GParam const * const param);
double local_gauge_div(Conf const * const GC,
                       Geometry const * const geo,
                       long int r);
double lorenz_gauge_violation(Conf const * const GC,
                              Geometry const * const geo,
                              GParam const * const param);
void compute_flavour_observables(Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp);
void compute_gauge_correlators(Conf const * const GC,
                               GParam const * const param,
                               double *tildeG1_p0,
                               double *tildeG1_pmin,
                               double *tildeG2_p0,
                               double *tildeG2_pmin,
                               double *tildeG3_p0,
                               double *tildeG3_pmin,
                               double *disc_p0,
                               double *disc_pmin);
void local_fix_lorenz_gauge(Conf *GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            long int r);
void fix_lorenz_gauge(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo);
void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep);

#endif
