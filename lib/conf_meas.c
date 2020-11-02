#ifndef CONF_MEAS_C
#define CONF_MEAS_C

#include"../include/macro.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/flavour_matrix.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/conf.h"

// computation of the plaquette in position r and positive directions i,j
double plaquette_single(Conf const * const GC,
                        Geometry const * const geo,
                        long r,
                        int i,
                        int j)
   {

//
//       ^ i
//       |  (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> j
//       r  (1)
//

   double ris;

   #ifdef CSTAR_BC
     ris = GC->theta[r][j];  // (1)
     ris += bcsitep(geo, r, j)*(GC->theta[nnp(geo, r, j)][i]); // (2)
     ris -= bcsitep(geo, r, i)*(GC->theta[nnp(geo, r, i)][j]); // (3)
     ris -= GC->theta[r][i];
   #else
     ris = GC->theta[r][j];  // (1)
     ris += GC->theta[nnp(geo, r, j)][i]; // (2)
     ris -= GC->theta[nnp(geo, r, i)][j]; // (3)
     ris -= GC->theta[r][i];
   #endif

   return ris;
   }


// computation of the (squared) plaquette in position r and positive directions i,j
double plaquettesq_single(Conf const * const GC,
                          Geometry const * const geo,
                          long r,
                          int i,
                          int j)
   {

//
//       ^ i
//       |  (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> j
//       r  (1)
//

   double ris;

   #ifdef CSTAR_BC
     ris = GC->theta[r][j];  // (1)
     ris += bcsitep(geo, r, j)*(GC->theta[nnp(geo, r, j)][i]); // (2)
     ris -= bcsitep(geo, r, i)*(GC->theta[nnp(geo, r, i)][j]); // (3)
     ris -= GC->theta[r][i];
   #else
     ris = GC->theta[r][j];  // (1)
     ris += GC->theta[nnp(geo, r, j)][i]; // (2)
     ris -= GC->theta[nnp(geo, r, i)][j]; // (3)
     ris -= GC->theta[r][i];
   #endif

   return ris*ris;
   }


// compute the average plaquette (squared)
double plaquettesq(Conf const * const GC,
                   Geometry const * const geo,
                   GParam const * const param)
   {
   long r;
   double ris=0.0;

   for(r=0; r<(param->d_volume); r++)
      {
      double tmp;
      int i, j;

      i=0;
      tmp=0.0;
     
      for(i=0; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            tmp+=plaquettesq_single(GC, geo, r, i, j);
            }
         }

      ris+=tmp;
      }

   ris*=param->d_inv_vol;
   ris/=((double) STDIM*((double) STDIM-1.0)/2.0);

   return ris;
   }


// compute the average value of Re[ phi_x^{dag} lambda_{x,mu} phi_{x+mu} ] with lambda_{x,mu}=e^{i theta_{x,mu}}
double higgs_interaction(Conf const * const GC,
                         Geometry const * const geo,
                         GParam const * const param)
  {
  int i;
  long r;
  double aux, ris=0.0;
  Vec v1;

  for(r=0; r<(param->d_volume); r++)
     {
     aux=0.0;

     for(i=0; i<STDIM; i++)
        {
        #ifdef CSTAR_BC
          if(bcsitep(geo, r, i)==1)
            {
            equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]));
            }
          else
            {
            equal_cc_Vec(&v1, &(GC->phi[nnp(geo, r, i)]));
            }
        #else
          equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]));
        #endif
        times_equal_complex_Vec(&v1, cexp(I*GC->theta[r][i]) );

        aux+= creal(scal_prod_Vec(&(GC->phi[r]), &v1) );
        }

     ris+=aux;
     }

  ris/=(double) STDIM;
  ris*=param->d_inv_vol;

  return ris;
  }


// compute flavour related observables
//
// GC->Qh needs to be initialized before calling this function
//
// tildeG0=Tr[(\sum_x Q_x)(\sum_y Q_y)]/volume
// tildeGminp=ReTr[(\sum_x Q_xe^{ipx})(\sum_y Q_ye^{-ipy)]/volume
//
// tildeG0 is the susceptibility, tildeGminp is used to compute the 2nd momentum correlation function
//
void compute_flavour_observables(Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)param->d_size[1];
  FMatrix Q, Qp, Qmp, tmp1, tmp2;

  // Q =sum_x Q_x
  // Qp=sum_x e^{ipx}Q_x
  // Qmp=sum_x e^{-ipx}Q_x

  zero_FMatrix(&Q);
  zero_FMatrix(&Qp);
  zero_FMatrix(&Qmp);
  for(r=0; r<(param->d_volume); r++)
     {
     equal_FMatrix(&tmp1, &(GC->Qh[r]));
     equal_FMatrix(&tmp2, &tmp1);

     plus_equal_FMatrix(&Q, &tmp1);

     si_to_cart(coord, r, param);

     times_equal_complex_FMatrix(&tmp1, cexp(I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qp, &tmp1);

     times_equal_complex_FMatrix(&tmp2, cexp(-I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qmp, &tmp2);
     }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0=retr_FMatrix(&tmp1)*param->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp=retr_FMatrix(&tmp1)*param->d_inv_vol;
  }


// compute field strenght correlations
//
// GC->F needs to be initialized before calling this function
//
// tildeF0, tildeFminp_long, tildeFminp_per are used to compute the 2nd momentum correlation function
// for the longitudinal and perpendicular components
//
void compute_plaq_corrlengths(Conf const * const GC,
                              GParam const * const param,
                              double *tildeF0,
                              double *tildeFminp_long,
                              double *tildeFminp_perp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)param->d_size[1];
  double tF;
  double complex tFl, tFp;

  #if STDIM>2
    if(param->d_size[1]!=param->d_size[2])
      {
      fprintf(stderr, "The function compute_plaq_corrlengths requires size[1]==size[2] (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
  #endif

  tF=0.0;
  tFl=0.0+0.0*I;
  tFp=0.0+0.0*I;

  // compute wall sums
  for(r=0; r<param->d_volume; r++)
     {
     si_to_cart(coord, r, param);

     tFl+=((double complex) GC->F[r]) * cexp(I * (double)coord[1] * p);
     #if STDIM>2
       tFp+=((double complex) GC->F[r]) * cexp(I * (double)coord[2] * p);
     #endif

     tF+=GC->F[r];
     }

  *tildeF0=tF*tF*param->d_inv_vol;
  *tildeFminp_long=creal(tFl*conj(tFl))*param->d_inv_vol;
  *tildeFminp_perp=creal(tFp*conj(tFp))*param->d_inv_vol;
  }


void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep)
   {
   long r;

   double tildeG0, tildeGminp;
   double scalar_coupling, plaqsq;
   double tildeF0, tildeFminp_long, tildeFminp_perp;


   for(r=0; r<(param->d_volume); r++)
      {
      init_FMatrix(&(GC->Qh[r]), &(GC->phi[r]));
      }

   compute_flavour_observables(GC,
                               param,
                               &tildeG0,
                               &tildeGminp);

   scalar_coupling=higgs_interaction(GC, geo, param);
   plaqsq=plaquettesq(GC, geo, param);

   for(r=0; r<param->d_volume; r++)
      {
      GC->F[r] = plaquette_single(GC, geo, r, 0, 1);
      }

   compute_plaq_corrlengths(GC,
                            param,
                            &tildeF0,
                            &tildeFminp_long,
                            &tildeFminp_perp);

   fprintf(datafilep, "%.12g %.12g %.12g %.12g ", tildeG0, tildeGminp, scalar_coupling, plaqsq);
   fprintf(datafilep, "%.12g %.12g %.12g", tildeF0, tildeFminp_long, tildeFminp_perp);
   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


#endif












