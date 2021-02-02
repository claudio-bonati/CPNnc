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

// compute some susceptibilities related to gauge fields
//
// chiA = (1/volume) | \sum_x \vec{A}_x e^(ip*x) |^2
// with p = (pi/L0, pi/L1, pi/L_2...) to take into account C* bc
//
// chiA2 = (1/volume) (sum_x \vec{A}_x^2)^2
void compute_gauge_susc(Conf const * const GC,
                        GParam const * const param,
                        double *chiA,
                        double *chiA2)
  {
  double complex vc[STDIM];
  double suma2;
  int i, r, coord[STDIM];

  suma2=0.0;
  for(i=0; i<STDIM; i++)
     {
     vc[i]=0.0+0.0*I;
     }

  for(r=0; r<(param->d_volume); r++)
     {
     si_to_cart(coord, r, param);

     for(i=0; i<STDIM; i++)
        {
        suma2+=(GC->theta[r][i])*(GC->theta[r][i]);

        vc[i]+=(GC->theta[r][i])*cexp(I * (double)coord[i] * PI / (double)param->d_size[i]);
        }
     }

  *chiA2=suma2*suma2/(double)param->d_volume;

  *chiA=0;
  for(i=0; i<STDIM; i++)
     {
     *chiA += cabs(vc[i])*cabs(vc[i])/(double) param->d_volume;
     }
  }

void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep)
   {
   long r;

   double tildeG0, tildeGminp;
   double scalar_coupling, plaqsq;

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

   fprintf(datafilep, "%.12g %.12g %.12g %.12g ", tildeG0, tildeGminp, scalar_coupling, plaqsq);

   #ifdef TEMPORAL_GAUGE
   double chiA, chiA2;
   compute_gauge_susc(GC,
                      param,
                      &chiA,
                      &chiA2);
   fprintf(datafilep, "%.12g %.12g ", chiA, chiA2);
   #endif

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


#endif












