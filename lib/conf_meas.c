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


// compute the violation of the Lorentz condition: \sum_{x} (\sum_{mu}\partial_{\mu} theta_{\mu})^2
double lorenz_gauge_violation(Conf const * const GC,
                              Geometry const * const geo,
                              GParam const * const param)
   {
   long r;
   int i;
   double ris, tmp;

   ris=0.0;
   for(r=0; r<param->d_volume; r++)
      {
      tmp=0.0;
      for(i=0; i<STDIM; i++)
         {
         #ifdef CSTAR_BC
           tmp += bcsitep(geo, r, i)*(GC->theta[nnp(geo, r, i)][i]);
           tmp -= GC->theta[r][i];
         #else
           tmp += GC->theta[nnp(geo, r, i)][i];
           tmp -= GC->theta[r][i];
         #endif
         }
      ris+=tmp*tmp;
      }

   return ris;
   }


// compute flavour related observables
//
// GC->Qh needs to be initialized before calling this function
//
// tildeG0=ReTr[(\sum_x Q_x)(\sum_y Q_y)]/volume
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
   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


// compute gauge dependent correlators
//
// tildeG1_p1=Re[(\sum_x A_{x,0}e^{ip_1x})(\sum_y A_{y,0}e^{-ip_1y)]/volume
// with p_1=(pi/L_0, pi/L_1, pi/L_2...)
// and analogously with p_2=(3pi/L_0, pi/L_1, pi/L_2...)
//
// tildeG2_p1=Re[(\sum_x A_{x,1}e^{ip_1x})(\sum_y A_{y,1}e^{-ip_1y)]/volume
// with p_1=(pi/L_0, pi/L_1, pi/L_2...)
// and analogously with p_2=(3pi/L_0, pi/L_1, pi/L_2...)
//
// B_x=sum_{mu} A_{x,mu}^2
// tildeG3_p0=Re[(\sum_x B_x)(\sum_y B_y)]/volume
// tildeG3_pmin=Re[(\sum_x B_xe^{i p_min*x})(\sum_y B_ye^{-ip_min*y)]/volume
// p_min=(2pi/L_0,0,0,...)
//
// disc_p0=[\sum_x B_x]/volume
// disc_pmin=Re[\sum_x B_x e^{i*p_min*x}]/volume
//
void compute_gauge_correlators(Conf const * const GC,
                               GParam const * const param,
                               double *tildeG1_p1,
                               double *tildeG1_p2,
                               double *tildeG2_p1,
                               double *tildeG2_p2,
                               double *tildeG3_p0,
                               double *tildeG3_pmin,
                               double *disc_p0,
                               double *disc_pmin)
  {
  int i, coord[STDIM];
  long r;
  double p1[STDIM], sc, B, forG3_p0;
  double complex forG1_p1, forG1_p2, forG2_p1, forG2_p2, forG3_pmin;

  for(i=0; i<STDIM; i++)
     {
     p1[i]=PI/(double)param->d_size[i];
     }

  forG1_p1=0.0+I*0.0;
  forG1_p2=0.0+I*0.0;
  forG2_p1=0.0+I*0.0;
  forG2_p2=0.0+I*0.0;
  forG3_p0=0.0+I*0.0;
  forG3_pmin=0.0+I*0.0;

  for(r=0; r<(param->d_volume); r++)
     {
     si_to_cart(coord, r, param);

     sc=0.0;
     for(i=0; i<STDIM; i++)
        {
        sc+=coord[i]*p1[i];
        }
     forG1_p1+=GC->theta[r][0]*cexp(I*sc);
     forG2_p1+=GC->theta[r][1]*cexp(I*sc);

     sc+=2.0*PI/param->d_size[0]*coord[0];
     forG1_p2+=GC->theta[r][0]*cexp(I*sc);
     forG2_p2+=GC->theta[r][1]*cexp(I*sc);

     B=0;
     for(i=0; i<STDIM; i++)
        {
        B+=(GC->theta[r][i])*(GC->theta[r][i]);
        }
     forG3_p0+=B;
     forG3_pmin+=B*cexp(I*((double) coord[0])*2.0*PI/(double)param->d_size[0]);
     }

  *tildeG1_p1=creal(forG1_p1*conj(forG1_p1))*param->d_inv_vol;
  *tildeG1_p2=creal(forG1_p2*conj(forG1_p2))*param->d_inv_vol;

  *tildeG2_p1=creal(forG2_p1*conj(forG2_p1))*param->d_inv_vol;
  *tildeG2_p2=creal(forG2_p2*conj(forG2_p2))*param->d_inv_vol;

  *tildeG3_p0=forG3_p0*forG3_p0*param->d_inv_vol;
  *tildeG3_pmin=creal(forG3_pmin*conj(forG3_pmin))*param->d_inv_vol;

  *disc_p0=forG3_p0*param->d_inv_vol;
  *disc_pmin=creal(forG3_pmin)*param->d_inv_vol;
  }


void perform_gaugedep_measures(Conf *GC,
                               GParam const * const param,
                               FILE *datafilep)
   {
   int i;
   double meas[8];

   compute_gauge_correlators(GC,
                             param,
                             &meas[0],
                             &meas[1],
                             &meas[2],
                             &meas[3],
                             &meas[4],
                             &meas[5],
                             &meas[6],
                             &meas[7]);
   for(i=0; i<8; i++)
      {
      fprintf(datafilep, "%.12g ", meas[i]);
      }
   fprintf(datafilep, "\n");
   fflush(datafilep);
   }



#endif
