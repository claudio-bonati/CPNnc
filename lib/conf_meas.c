#ifndef CONF_MEAS_C
#define CONF_MEAS_C

#include"../include/macro.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/conf.h"
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

  zero_Vec(&v1);

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

// compute \sum_{mu}(theta_{r,\mu}-theta_{r-mu,mu})
double local_gauge_div(Conf const * const GC,
                       Geometry const * const geo,
                       long int r)
   {
   int i;
   double ris;

   ris=0.0;
   for(i=0; i<STDIM; i++)
      {
      #ifdef CSTAR_BC
        ris += GC->theta[r][i];
        ris -= bcsitem(geo, r, i)*(GC->theta[nnm(geo, r, i)][i]);
      #else
        ris += GC->theta[r][i];
        ris -= GC->theta[nnm(geo, r, i)][i];
      #endif
      }

   return ris;
   }


// compute the violation of the Lorentz condition: \sum_{x} (\sum_{mu}\partial_{\mu} theta_{\mu})^2
double lorenz_gauge_violation(Conf const * const GC,
                              Geometry const * const geo,
                              GParam const * const param)
   {
   long r;
   double ris, tmp;

   ris=0.0;
   for(r=0; r<param->d_volume; r++)
      {
      tmp=local_gauge_div(GC, geo, r);
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
  const double p = 2.0*PI/(double)param->d_size[0];
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

     times_equal_complex_FMatrix(&tmp1, cexp(I*((double)coord[0])*p));
     plus_equal_FMatrix(&Qp, &tmp1);

     times_equal_complex_FMatrix(&tmp2, cexp(-I*((double)coord[0])*p));
     plus_equal_FMatrix(&Qmp, &tmp2);
     }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0=retr_FMatrix(&tmp1)*param->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp=retr_FMatrix(&tmp1)*param->d_inv_vol;
  }


// for C*
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
// C(r)=Re(phi[r])+e^{ip1*x}Im(phi[r]) is periodic
// tildeG4_p0=rescalprod[(\sum_x C_x),(\sum_y C_y)]/volume
// tildeG4_pmin=rescalprod[(\sum_x C_xe^{i p_min*x}),(\sum_y C_ye^{-ip_min*y)]/volume
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
                               double *disc_pmin,
                               double *tildeG4_p0,
                               double *tildeG4_pmin)
  {
  #ifndef CSTAR_BC
    fprintf(stderr, "This function can be used only with C* b.c.! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  #endif

  int i, coord[STDIM];
  long r;
  double p1[STDIM], sc, sc2, B, forG3_p0;
  double complex forG1_p1, forG1_p2, forG2_p1, forG2_p2, forG3_pmin;
  Vec forG4_p0, forG4_pmin, revec, imvec;

  // min momentum for antiperiodic boundary conditions
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
  zero_Vec(&forG4_p0);
  zero_Vec(&forG4_pmin);

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

     // momentum p1+(2pi/L,0....0)
     sc2=sc+2.0*PI/param->d_size[0]*coord[0];

     forG1_p2+=GC->theta[r][0]*cexp(I*sc2);
     forG2_p2+=GC->theta[r][1]*cexp(I*sc2);

     B=0;
     for(i=0; i<STDIM; i++)
        {
        B+=(GC->theta[r][i])*(GC->theta[r][i]);
        }
     forG3_p0+=B;
     forG3_pmin+=B*cexp(I*((double) coord[0])*2.0*PI/(double)param->d_size[0]);

     repart_Vec(&revec, &(GC->phi[r]));
     impart_Vec(&imvec, &(GC->phi[r]));
     times_equal_complex_Vec(&imvec, I*cexp(I*sc));
     plus_equal_Vec(&revec, &imvec); // now revec=Re[phi]+e^{i*p1*r}I*Im[phi]
     plus_equal_Vec(&forG4_p0, &revec);
     times_equal_complex_Vec(&revec, cexp(I*((double) coord[0])*2.0*PI/(double)param->d_size[0]) );
     plus_equal_Vec(&forG4_pmin, &revec);
     }

  *tildeG1_p1=creal(forG1_p1*conj(forG1_p1))*param->d_inv_vol;
  *tildeG1_p2=creal(forG1_p2*conj(forG1_p2))*param->d_inv_vol;

  *tildeG2_p1=creal(forG2_p1*conj(forG2_p1))*param->d_inv_vol;
  *tildeG2_p2=creal(forG2_p2*conj(forG2_p2))*param->d_inv_vol;

  *tildeG3_p0=forG3_p0*forG3_p0*param->d_inv_vol;
  *tildeG3_pmin=creal(forG3_pmin*conj(forG3_pmin))*param->d_inv_vol;

  *disc_p0=forG3_p0*param->d_inv_vol;
  *disc_pmin=creal(forG3_pmin)*param->d_inv_vol;

  *tildeG4_p0=creal(scal_prod_Vec(&forG4_p0, &forG4_p0))*param->d_inv_vol;
  *tildeG4_pmin=creal(scal_prod_Vec(&forG4_pmin, &forG4_pmin))*param->d_inv_vol;
  }


// out=A*in, where ''A'' is the matrix needed to fix the Lorenz gauge
void matrix_apply(double *out,
                  double const * const in,
                  GParam const * const param,
                  Geometry const * const geo)
  {
  int i, j;
  long r, r1;
  double tmp;

  for(r=0; r<param->d_volume; r++)
     {
     out[r]=4*STDIM*STDIM*in[r];

     tmp=0.0;
     for(i=0; i<STDIM; i++)
        {
        #ifdef CSTAR_BC
          tmp+=bcsitep(geo, r, i)*in[nnp(geo, r, i)];
        #else
          tmp+=in[nnp(geo, r, i)];
        #endif
        }
     out[r]-=4.0*STDIM*tmp;

     tmp=0.0;
     for(i=0; i<STDIM; i++)
        {
        #ifdef CSTAR_BC
          tmp+=bcsitem(geo, r, i)*in[nnm(geo, r, i)];
        #else
          tmp+=in[nnm(geo, r, i)];
        #endif
        }
     out[r]-=4.0*STDIM*tmp;

     // +i -j
     tmp=0.0;
     for(i=0; i<STDIM; i++)
        {
        r1=nnp(geo, r, i);
        for(j=0; j<STDIM; j++)
           {
           #ifdef CSTAR_BC
             tmp+=bcsitep(geo, r, i)*bcsitem(geo, r1, j)*in[nnm(geo, r1, j)];
           #else
             tmp+=in[nnm(geo, r1, j)];
           #endif
           }
        }
     out[r]+=tmp;

     // +i+j
     tmp=0.0;
     for(i=0; i<STDIM; i++)
        {
        r1=nnp(geo, r, i);
        for(j=0; j<STDIM; j++)
           {
           #ifdef CSTAR_BC
             tmp+=bcsitep(geo, r, i)*bcsitep(geo, r1, j)*in[nnp(geo, r1, j)];
           #else
             tmp+=in[nnp(geo, r1, j)];
           #endif
           }
        }
     out[r]+=tmp;

     //-i-j
     tmp=0.0;
     for(i=0; i<STDIM; i++)
        {
        r1=nnm(geo, r, i);
        for(j=0; j<STDIM; j++)
           {
           #ifdef CSTAR_BC
             tmp+=bcsitem(geo, r, i)*bcsitem(geo, r1, j)*in[nnm(geo, r1, j)];
           #else
             tmp+=in[nnm(geo, r1, j)];
           #endif
           }
        }
     out[r]+=tmp;

     // -i+j
     tmp=0.0;
     for(i=0; i<STDIM; i++)
        {
        r1=nnm(geo, r, i);
        for(j=0; j<STDIM; j++)
           {
           #ifdef CSTAR_BC
             tmp+=bcsitem(geo, r, i)*bcsitep(geo, r1, j)*in[nnp(geo, r1, j)];
           #else
             tmp+=in[nnp(geo, r1, j)];
           #endif
           }
        }
     out[r]+=tmp;
     }
  }


// solve Ax=b using conjugate gradient with ''A'' the matrix needed to fix the Lorenz gauge
// *x has to be initialized before calling this function
void CG_solver(double *x,
               double const * const b,
               GParam const * const param,
               Geometry const * const geo,
               double soglia)
  {
  int k, err;
  long r;
  double alphak, betak, *pk, *rk, *Apk, norm1, norm1new, norm2;
  const int maxiteration=10000;

  err=posix_memalign((void**)&pk, (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in CG_solver! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  err=posix_memalign((void**)&rk, (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in CG_solver! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  err=posix_memalign((void**)&Apk, (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in CG_solver! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  matrix_apply(pk, x, param, geo);  // pk = A*x
  for(r=0; r<param->d_volume; r++)  // pk=rk=b-A*x
     {
     pk[r]=b[r]-pk[r];
     rk[r]=pk[r];
     }

  norm1=0.0;     // norm_1=(r_k, r_k)
  for(r=0; r<param->d_volume; r++)
     {
     norm1+=rk[r]*rk[r];
     }

  k=0; // iteration index
  while(sqrt(norm1)>soglia && k<maxiteration)
       {
       matrix_apply(Apk, pk, param, geo);  // Apk = A*p_k

       norm2=0.0;  // norm2=(p_k, A p_k)
       for(r=0; r<param->d_volume; r++)
          {
          norm2+=pk[r]*Apk[r];
          }

       alphak=norm1/norm2; // alpha_k=(r_k,r_k)/(p_k, A*p_k)

       for(r=0; r<param->d_volume; r++)
          {
          x[r]+=alphak*pk[r];    // x_k->x_{k+1}
          rk[r]-=alphak*Apk[r];  // r_k->r_{k+1}
          }

       norm1new=0.0;  // norm1new = (r_{k+1}, r_{k+1})
       for(r=0; r<param->d_volume; r++)
          {
          norm1new+=rk[r]*rk[r];
          }

       betak=norm1new/norm1;  // beta_k = (r_{k+1}, r_{k+1})/(r_k, r_k)
       norm1=norm1new;        // norm1=(r_{k+1},r_{k+1})

       for(r=0; r<param->d_volume; r++)
          {
          pk[r]=rk[r]+betak*pk[r];  // p_k->p_{k+1}
          }

       #ifdef DEBUG
         printf("CG: iteration %d (max %d), normalized residual %g\n", k, maxiteration, sqrt(norm1)/soglia);
       #endif

       k+=1;
       }

  if(k>=maxiteration)
    {
    fprintf(stderr, "Max number of iterations reached in CG_solver! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef DEBUG
    matrix_apply(pk, x, param, geo);  // pk =A*x

    for(r=0; r<param->d_volume; r++)  // rk = b-A*x
       {
       rk[r]=b[r]-pk[r];
       }

    norm1=0.0;
    for(r=0; r<param->d_volume; r++)
       {
       norm1+=rk[r]*rk[r];
       }
    printf("CG final normalized residual %g\n\n", sqrt(norm1)/soglia);

    if(norm1 >= 2*soglia)
      {
      fprintf(stderr, "CG solver final precision test FAILED! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
  #endif

  free(pk);
  free(rk);
  free(Apk);
  }


// fix lorentz gauge via conjugate gradient
void fix_lorenz_gauge_conjgrad(Conf *GC,
                               GParam const * const param,
                               Geometry const * const geo)
  {
  int i, err;
  long int r;
  double *aux, *b, *sol;
  const double soglia=MIN_VALUE*sqrt((double)param->d_volume);
  double sogliaCG, test;

  #ifdef DEBUG
    double test1, test2, test1new;

    test1=plaquettesq(GC, geo, param);
  #endif

  err=posix_memalign((void**)&aux, (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in lorenz_gauge_conjgrad! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  err=posix_memalign((void**)&b, (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in lorenz_gauge_conjgrad! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  err=posix_memalign((void**)&sol, (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in lorenz_gauge_conjgrad! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // aux is the vector of the divergences
  for(r=0; r<param->d_volume; r++)
     {
     aux[r]=local_gauge_div(GC, geo, r);
     }

  // b is the r.h.s. of the linear equation to be solved
  for(r=0; r<param->d_volume; r++)
     {
     b[r]=-2*STDIM*aux[r];
     for(i=0; i<STDIM; i++)
        {
        #ifdef CSTAR_BC
          b[r]+=bcsitep(geo,r,i)*aux[nnp(geo, r, i)];
          b[r]+=bcsitem(geo,r,i)*aux[nnm(geo, r, i)];
        #else
          b[r]+=aux[nnp(geo, r, i)];
          b[r]+=aux[nnm(geo, r, i)];
        #endif
        }
     }

  // free the auxilliary vector
  free(aux);

  // initialize to zero vector
  for(r=0; r<param->d_volume; r++)
     {
     sol[r]=0.0;
     }

  sogliaCG=soglia;
  test=lorenz_gauge_violation(GC, geo, param);
  while(test/soglia>1)
       {
       // find the gauge transformation
       CG_solver(sol, b, param, geo, sqrt(sogliaCG));

       // apply the gauge transformation
       for(r=0; r<param->d_volume; r++)
          {
          for(i=0; i<STDIM; i++)
             {
             #ifdef CSTAR_BC
               GC->theta[r][i]+=sol[r];
               GC->theta[r][i]-=bcsitep(geo, r, i)*sol[nnp(geo, r, i)];
             #else
               GC->theta[r][i]+=sol[r];
               GC->theta[r][i]-=sol[nnp(geo, r, i)];
             #endif
             }
          times_equal_complex_Vec(&(GC->phi[r]), cexp(I*sol[r]));
          }

       test=lorenz_gauge_violation(GC, geo, param);

       sogliaCG/=25.0;
       }

  #ifdef DEBUG
    test1new=plaquettesq(GC, geo, param);
    test2=lorenz_gauge_violation(GC, geo, param);

    if(fabs(test1-test1new)>MIN_VALUE)
      {
      fprintf(stderr, "Problems in fix_lorenz_gauge_conjgrad! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    printf("Normalized Lorenz gauge violation after Lorenz gauge fixing: %g\n\n", test2/soglia);
  #endif


  free(b);
  free(sol);
  }


// perform the measures
void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep)
   {
   long r;
   int i;

   double tildeG0, tildeGminp;
   double scalar_coupling, plaqsq;
   double meas[10];

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

   #ifdef HARD_LORENZ_GAUGE  // fix lorenz gauge in measures
     Conf GCbis;

     init_conf_from_conf(&GCbis, GC, param);
     fix_lorenz_gauge_conjgrad(&GCbis, param, geo);


     // gauge dependent measures
     compute_gauge_correlators(&GCbis,
                               param,
                               &meas[0],
                               &meas[1],
                               &meas[2],
                               &meas[3],
                               &meas[4],
                               &meas[5],
                               &meas[6],
                               &meas[7],
                               &meas[8],
                               &meas[9]);


     free_conf(&GCbis, param);
   #else
     // gauge dependent measures
     compute_gauge_correlators(GC,
                               param,
                               &meas[0],
                               &meas[1],
                               &meas[2],
                               &meas[3],
                               &meas[4],
                               &meas[5],
                               &meas[6],
                               &meas[7],
                               &meas[8],
                               &meas[9]);
   #endif

   for(i=0; i<10; i++)
      {
      fprintf(datafilep, "%.12g ", meas[i]);
      }

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

#endif
