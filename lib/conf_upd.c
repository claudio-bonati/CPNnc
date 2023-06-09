#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdlib.h>

#include"../include/conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

// compute the staple for the phi field, such that
// \sum_{mu>0 and <0} lambda_{x,mu} conj(\phi_x).\phi_{x+\mu} = conj(\phi_x).staple
// with lambda_{x,mu}=e^{i theta[x][mu]}
void calcstaples_for_phi(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         Vec *staple)
  {
  int i;
  double complex lambda;
  Vec v1, v2;

  zero_Vec(staple);

  for(i=0; i<STDIM; i++)
     {
     // forward
     #ifdef CSTAR_BC
       if(bcsitep(geo, r, i)==1)
         {
         equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]) );
         }
       else
         {
         equal_cc_Vec(&v1, &(GC->phi[nnp(geo, r, i)]) );
         }
     #else
       equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]) );
     #endif

     times_equal_complex_Vec(&v1, cexp(I*GC->theta[r][i]));
     plus_equal_Vec(staple, &v1);

     // backward
     #ifdef CSTAR_BC
       if(bcsitem(geo, r, i)==1)
         {
         equal_Vec(&v2, &(GC->phi[nnm(geo, r, i)]) );
         lambda=cexp(-I*GC->theta[nnm(geo, r, i)][i]);
         }
       else
         {
         equal_cc_Vec(&v2, &(GC->phi[nnm(geo, r, i)]) );
         lambda=cexp(I*GC->theta[nnm(geo, r, i)][i]);
         }
     #else
       equal_Vec(&v2, &(GC->phi[nnm(geo, r, i)]) );
       lambda=cexp(-I*GC->theta[nnm(geo, r, i)][i]);
     #endif
     times_equal_complex_Vec(&v2, lambda);

     plus_equal_Vec(staple, &v2);
     }
  }



// perform an update of the phi field with metropolis
// retrn 0 if the trial state is rejected and the number of accepted trial
int metropolis_for_phi(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r)
  {
  int acc=0;
  double old_energy, new_energy;
  Vec staple, new_vector;

  calcstaples_for_phi(GC, geo, r, &staple);

  old_energy=-2.0 * (double) NFLAVOUR * param->d_J * creal(scal_prod_Vec(&(GC->phi[r]), &staple));

  rand_rot_single_Vec(&new_vector, &(GC->phi[r]), param->d_epsilon_metro_site);

  new_energy=-2.0 * (double) NFLAVOUR * param->d_J * creal(scal_prod_Vec(&new_vector, &staple));

  #ifdef DEBUG
  double old_energy_aux, new_energy_aux;
  Vec old_vector;
  old_energy_aux= -2.0 * (double) NFLAVOUR * param->d_J * higgs_interaction(GC, geo, param) * (double) STDIM * (double) param->d_volume;
  equal_Vec(&old_vector, &(GC->phi[r]));
  equal_Vec(&(GC->phi[r]), &new_vector);
  new_energy_aux= -2.0 * (double) NFLAVOUR * param->d_J * higgs_interaction(GC, geo, param) * (double) STDIM * (double) param->d_volume;
  equal_Vec(&(GC->phi[r]), &old_vector);
  //printf("%g %g\n", old_energy-new_energy, old_energy-new_energy -(old_energy_aux-new_energy_aux));
  if(fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux))>1.0e-10 )
    {
    fprintf(stderr, "Problem in energy in metropolis for phi (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  if(old_energy>new_energy)
    {
    equal_Vec(&(GC->phi[r]), &new_vector);
    acc+=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         equal_Vec(&(GC->phi[r]), &new_vector);
         acc+=1;
         }

  #ifdef DEBUG
  if(fabs(norm_Vec(&(GC->phi[r]))-1)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in metropolis (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  return acc;
  }


// perform an update of the phi field with overrelaxation
void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            long r)
  {
  double norm;
  double complex aux1;
  Vec staple, newlink;

  calcstaples_for_phi(GC, geo, r, &staple);
  norm=norm_Vec(&staple);

  #ifdef DEBUG
  double prod_before=creal(scal_prod_Vec(&(GC->phi[r]), &staple));
  #endif

  if(norm>MIN_VALUE)
    {
    equal_Vec(&newlink, &staple);
    times_equal_real_Vec(&newlink, 1./norm);
    aux1=scal_prod_Vec(&newlink, &(GC->phi[r]) );
    times_equal_complex_Vec(&newlink, 2.0*aux1);

    minus_equal_Vec(&newlink, &(GC->phi[r]));

    equal_Vec(&(GC->phi[r]), &newlink);
    }
  else
    {
    rand_vec_Vec(&(GC->phi[r]) );
    }

  #ifdef DEBUG
  double prod_after=creal(scal_prod_Vec(&(GC->phi[r]), &staple));

  if(fabs(prod_before-prod_after)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in overrelaxation (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  if(fabs(norm_Vec(&(GC->phi[r]))-1)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in overrelaxation (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif
  }


// staples for the plaquette component of the action
// sum (plaq)^2 = 2*(STDIM-1)*theta^2 + 2*theta*plaqstaple + independent of theta
double plaqstaples_for_link(Conf *GC,
                            Geometry const * const geo,
                            long r,
                            int i)
  {
  int j, k;
  double ris;
  long r1;

  ris=0.0;

  for(j=i+1; j<i+STDIM; j++)
     {
     k=j%STDIM;

//              ^ i
//         (6)  |  (3)
//      +---<---+--->---+
//      |       |       |
//   (5)V       ^       V (2)
//      |       |       |
//      +--->---+---<---+---> k
//     r1  (4)  r  (1)

       #ifdef CSTAR_BC
         //         (1)                          (2)                                           (3)
         ris+= -GC->theta[r][k] -bcsitep(geo, r, k)*(GC->theta[nnp(geo, r, k)][i]) +bcsitep(geo, r, i)*(GC->theta[nnp(geo,r,i)][k]);
         r1=nnm(geo, r, k);
         //          (4)                                           (5)                                   (6)
         ris += bcsitem(geo,r,k)*(GC->theta[r1][k]) -bcsitem(geo,r,k)*(GC->theta[r1][i]) -bcsitep(geo, r1, i)*bcsitem(geo,r,k)*(GC->theta[nnp(geo, r1, i)][k]);
       #else
         ris+= -GC->theta[r][k] -GC->theta[nnp(geo, r, k)][i] +GC->theta[nnp(geo,r,i)][k];
         r1=nnm(geo, r, k);
         ris += GC->theta[r1][k] -GC->theta[r1][i] -GC->theta[nnp(geo, r1, i)][k];
       #endif
       }

    return ris;
    }


// staples for the lorenz gauge term of the action
// sum_x (\sum_j[\theta_{x,j}-theta_{x-j,j}])^2 = 2*theta_{r,i}^2 + 2*theta_{r,i}*lorenzstap + indip(theta_{r,i})
double lorenzstaples_for_link(Conf *GC,
                              Geometry const * const geo,
                              long r,
                              int i)
  {
  double ris;
  long r1;
  int k, j;

  #ifdef CSTAR_BC
    int sign;

    ris=0.0;
    // from r
    for(k=1; k<STDIM; k++)
       {
       j=(i+k)%STDIM;
       ris+= GC->theta[r][j] - bcsitem(geo,r,j)*GC->theta[nnm(geo,r,j)][j];
       }
    ris-= bcsitem(geo, r, i)*GC->theta[nnm(geo,r,i)][i];

    r1=nnp(geo,r,i);
    sign=bcsitep(geo,r,i);
    // from r+i
    for(k=1; k<STDIM; k++)
       {
       j=(i+k)%STDIM;
       ris-= sign*GC->theta[r1][j] - sign*bcsitem(geo,r1,j)*GC->theta[nnm(geo,r1,j)][j];
       }
    ris-=sign*GC->theta[r1][i];
  #else
    ris=0.0;
    // from r
    for(k=1; k<STDIM; k++)
       {
       j=(i+k)%STDIM;
       ris+= GC->theta[r][j] - GC->theta[nnm(geo,r,j)][j];
       }
    ris-=GC->theta[nnm(geo,r,i)][i];

    r1=nnp(geo,r,i);
    // from r+i
    for(k=1; k<STDIM; k++)
       {
       j=(i+k)%STDIM;
       ris-= GC->theta[r1][j] - GC->theta[nnm(geo,r1,j)][j];
       }
    ris-=GC->theta[r1][i];
  #endif

  return ris;
  }


// perform an update with metropolis of the link variables
// return 0 if the trial state is rejected and 1 otherwise
int metropolis_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i)
  {
  double old_energy, new_energy;
  double old_theta, new_theta;
  double complex sc;
  double pstaple;
  int acc=0;

  #ifdef SOFT_LORENZ_GAUGE
  double lstaple=lorenzstaples_for_link(GC, geo, r, i);
  #endif

  Vec v1;

  #ifdef CSTAR_BC
    if(bcsitep(geo, r, i)==1)
      {
      equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
    else
      {
      equal_cc_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
  #else
    equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
  #endif

  sc=scal_prod_Vec(&(GC->phi[r]), &v1);

  old_theta=GC->theta[r][i];

  if(fabs(param->d_K)>MIN_VALUE)
    {
    pstaple=plaqstaples_for_link(GC, geo, r, i);
    }
  else
    {
    pstaple=0.0;
    }

  old_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*cexp(I*old_theta) );
  old_energy+=0.5*param->d_K*(2.0*((double)STDIM-1.0)*old_theta*old_theta + 2.0*old_theta*pstaple);
  // we used sum (plaq)^2 = 2*(STDIM-1)*theta^2 + 2*theta*plaqstaple + independent of theta
  old_energy+=0.5 * param->d_phmass * param->d_phmass * old_theta * old_theta;
  #ifdef SOFT_LORENZ_GAUGE
    // sum_x (sum_i[\theta_{x,i}-theta_{x-i,i}])^2 = 2 theta^2 + 2*theta*lorenzstap + indip. theta
    old_energy+= param->d_gaugefixpar * old_theta*old_theta + param->d_gaugefixpar*old_theta*lstaple;
  #endif
  #ifdef SOFT_TEMPORAL_GAUGE
    if(i==0)
      {
      old_energy+= (param->d_gaugefixpar/2.0)*old_theta*old_theta;
      }
  #endif

  new_theta = old_theta + param->d_epsilon_metro_link*(2.0*casuale()-1);

  new_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*cexp(I*new_theta) );
  new_energy+=0.5*param->d_K*(2.0*((double)STDIM-1.0)*new_theta*new_theta + 2.0*new_theta*pstaple);
  new_energy+=0.5 * param->d_phmass * param->d_phmass * new_theta * new_theta;
  #ifdef SOFT_LORENZ_GAUGE
    new_energy+= param->d_gaugefixpar*new_theta*new_theta + param->d_gaugefixpar*new_theta*lstaple;
  #endif
  #ifdef SOFT_TEMPORAL_GAUGE
    if(i==0)
      {
      new_energy+= (param->d_gaugefixpar/2.0)*new_theta*new_theta;
      }
  #endif

  #ifdef DEBUG
    double old_energy_aux, new_energy_aux;
    old_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
    old_energy_aux +=(0.5*param->d_K)*plaquettesq(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
    old_energy_aux += 0.5 * param->d_phmass * param->d_phmass * old_theta * old_theta;
    #ifdef SOFT_LORENZ_GAUGE
      old_energy_aux += 0.5 * param->d_gaugefixpar * lorenz_gauge_violation(GC, geo, param);
    #endif
    #ifdef SOFT_TEMPORAL_GAUGE
      if(i==0)
        {
        old_energy_aux += (param->d_gaugefixpar/2.0)*old_theta*old_theta;
        }
    #endif

    GC->theta[r][i] = new_theta;
    new_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
    new_energy_aux += (0.5*param->d_K)*plaquettesq(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
    new_energy_aux += 0.5 * param->d_phmass * param->d_phmass * new_theta * new_theta;
    #ifdef SOFT_LORENZ_GAUGE
      new_energy_aux += 0.5 * param->d_gaugefixpar * lorenz_gauge_violation(GC, geo, param);
    #endif
    #ifdef SOFT_TEMPORAL_GAUGE
      if(i==0)
        {
        new_energy_aux += (param->d_gaugefixpar/2.0)*new_theta*new_theta;
        }
    #endif
    GC->theta[r][i] = old_theta;

    //printf("%g %g\n", old_energy-new_energy, old_energy-new_energy -(old_energy_aux-new_energy_aux));
    if(fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux))>1.0e-10 )
      {
      fprintf(stderr, "Problem in energy in metropolis for link %g (%s, %d)\n", fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux)), __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
  #endif

  if(old_energy>new_energy)
    {
    GC->theta[r][i] = new_theta;
    acc=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         GC->theta[r][i] = new_theta;
         acc=1;
         }

  return acc;
  }


// perform an update with metropolis of the link variables of 2pi step
// return 0 if the trial state is rejected and 1 otherwise
int metropolis_for_link_twopi(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i)
  {
  double old_energy, new_energy;
  double old_theta, new_theta;
  double complex sc;
  double pstaple;
  double tmp;
  int shift, maxshift=2;
  int acc=0;

  #ifdef SOFT_LORENZ_GAUGE
  double lstaple=lorenzstaples_for_link(GC, geo, r, i);
  #endif

  Vec v1;

  #ifdef CSTAR_BC
    if(bcsitep(geo, r, i)==1)
      {
      equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
    else
      {
      equal_cc_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
      }
  #else
    equal_Vec(&v1, &(GC->phi[nnp(geo,r,i)]));
  #endif

  sc=scal_prod_Vec(&(GC->phi[r]), &v1);

  old_theta=GC->theta[r][i];

  if(fabs(param->d_K)>MIN_VALUE)
    {
    pstaple=plaqstaples_for_link(GC, geo, r, i);
    }
  else
    {
    pstaple=0.0;
    }

  old_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*cexp(I*old_theta) );
  old_energy+=0.5*param->d_K*(2.0*((double)STDIM-1.0)*old_theta*old_theta + 2.0*old_theta*pstaple);
  // we used sum (plaq)^2 = 2*(STDIM-1)*theta^2 + 2*theta*plaqstaple + independent of theta
  old_energy+=0.5 * param->d_phmass * param->d_phmass * old_theta * old_theta;
  #ifdef SOFT_LORENZ_GAUGE
    // sum_x (sum_i[\theta_{x,i}-theta_{x-i,i}])^2 = 2 theta^2 + 2*theta*lorenzstap + indip. theta
    old_energy+= param->d_gaugefixpar * old_theta*old_theta + param->d_gaugefixpar*old_theta*lstaple;
  #endif
  #ifdef SOFT_TEMPORAL_GAUGE
    if(i==0)
      {
      old_energy+= (param->d_gaugefixpar/2.0)*old_theta*old_theta;
      }
  #endif

  if(casuale()>0.5)
    {
    tmp=1;
    }
  else
    {
    tmp=-1;
    }
  shift=1+(int)((double)(maxshift+0.1)*casuale());
  new_theta = old_theta + tmp*PI2*(double)shift;

  new_energy=-2.0*(double)NFLAVOUR*(param->d_J)*creal(sc*cexp(I*new_theta) );
  new_energy+=0.5*param->d_K*(2.0*((double)STDIM-1.0)*new_theta*new_theta + 2.0*new_theta*pstaple);
  new_energy+=0.5 * param->d_phmass * param->d_phmass * new_theta * new_theta;
  #ifdef SOFT_LORENZ_GAUGE
    new_energy+= param->d_gaugefixpar*new_theta*new_theta + param->d_gaugefixpar*new_theta*lstaple;
  #endif
  #ifdef SOFT_TEMPORAL_GAUGE
    if(i==0)
      {
      new_energy+= (param->d_gaugefixpar/2.0)*new_theta*new_theta;
      }
  #endif

  #ifdef DEBUG
    double old_energy_aux, new_energy_aux;
    old_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
    old_energy_aux +=(0.5*param->d_K)*plaquettesq(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
    old_energy_aux += 0.5 * param->d_phmass * param->d_phmass * old_theta * old_theta;
    #ifdef SOFT_LORENZ_GAUGE
      old_energy_aux += 0.5 * param->d_gaugefixpar * lorenz_gauge_violation(GC, geo, param);
    #endif
    #ifdef SOFT_TEMPORAL_GAUGE
      if(i==0)
        {
        old_energy_aux += (param->d_gaugefixpar/2.0)*old_theta*old_theta;
        }
    #endif

    GC->theta[r][i] = new_theta;
    new_energy_aux = -2.0 * (double)NFLAVOUR *(param->d_J)*higgs_interaction(GC, geo, param)*(double)STDIM * (double)param->d_volume;
    new_energy_aux += (0.5*param->d_K)*plaquettesq(GC, geo, param)*(double)STDIM*((double)STDIM-1.0)/2.0 *(double) param->d_volume;
    new_energy_aux += 0.5 * param->d_phmass * param->d_phmass * new_theta * new_theta;
    #ifdef SOFT_LORENZ_GAUGE
      new_energy_aux += 0.5 * param->d_gaugefixpar * lorenz_gauge_violation(GC, geo, param);
    #endif
    #ifdef SOFT_TEMPORAL_GAUGE
      if(i==0)
        {
        new_energy_aux += (param->d_gaugefixpar/2.0)*new_theta*new_theta;
        }
    #endif
    GC->theta[r][i] = old_theta;

    //printf("%g %g\n", old_energy-new_energy, old_energy-new_energy -(old_energy_aux-new_energy_aux));
    if(fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux))>1.0e-10 )
      {
      fprintf(stderr, "Problem in energy in metropolis for link %g (%s, %d)\n", fabs(old_energy-new_energy -(old_energy_aux-new_energy_aux)), __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
  #endif

  if(old_energy>new_energy)
    {
    GC->theta[r][i] = new_theta;
    acc=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         GC->theta[r][i] = new_theta;
         acc=1;
         }

  return acc;
  }


// perform a complete update
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc_site,
            double *acc_link,
            double *acc_link_twopi)
   {
   long r, asum_site, asum_link, asum_link_twopi;
   int j, dir;

   // metropolis on links
   asum_link=0;
   asum_link_twopi=0;

   #ifndef LINKS_FIXED_TO_ONE
     #ifndef HARD_TEMPORAL_GAUGE
     for(r=0; r<param->d_volume; r++)
        {
        for(dir=0; dir<STDIM; dir++)
           {
           asum_link+=metropolis_for_link(GC, geo, param, r, dir);
           asum_link_twopi+=metropolis_for_link_twopi(GC, geo, param, r, dir);
           }
        }
     #else
      for(r=0; r<param->d_volume; r++)
        {
        for(dir=1; dir<STDIM; dir++)
           {
           asum_link+=metropolis_for_link(GC, geo, param, r, dir);
           asum_link_twopi+=metropolis_for_link_twopi(GC, geo, param, r, dir);
           }
        }
     #endif
   #endif
   *acc_link=((double)asum_link)*param->d_inv_vol;
   *acc_link_twopi=((double)asum_link_twopi)*param->d_inv_vol;

   #ifndef HARD_TEMPORAL_GAUGE
   *acc_link/=(double)STDIM;
   *acc_link_twopi/=(double)STDIM;
   #else
   *acc_link/=(double)(STDIM-1);
   *acc_link_twopi/=(double)(STDIM-1);
   #endif

   // metropolis on phi
   asum_site=0;
   for(r=0; r<param->d_volume; r++)
      {
      asum_site+=metropolis_for_phi(GC, geo, param, r);
      }
   *acc_site=((double)asum_site)*param->d_inv_vol;

   // overrelax on phi
   for(j=0; j<param->d_overrelax; j++)
      {
      for(r=0; r<(param->d_volume); r++)
         {
         overrelaxation_for_phi(GC, geo, r);
         }
      }

   // final unitarization
   for(r=0; r<(param->d_volume); r++)
      {
      unitarize_Vec(&(GC->phi[r]));
      }

   // this prevents theta to overflow if d_K=0 (no kinetic term for theta) and d_phmass=0
   #if ( !defined(SOFT_LORENZ_GAUGE) && !defined(HARD_LORENZ_GAUGE) )
     #if !defined(SOFT_TEMPORAL_GAUGE)               // for HARD_TEMPORAL_GAUGE or no gauge fixing at all
       if(fabs(param->d_K)<MIN_VALUE && fabs(param->d_phmass) < MIN_VALUE)
         {
         for(r=0; r<param->d_volume; r++)
            {
            for(dir=0; dir<STDIM; dir++)
               {
               while(GC->theta[r][dir]>PI2)
                    {
                    GC->theta[r][dir]-=PI2;
                    }
               while(GC->theta[r][dir]<-PI2)
                    {
                    GC->theta[r][dir]+=PI2;
                    }
               }
            }
         }
     #else // if SOFT_TEMPORAL_GAUGE
       if(fabs(param->d_K)<MIN_VALUE && fabs(param->d_phmass) < MIN_VALUE)
         {
         for(r=0; r<param->d_volume; r++)
            {
            for(dir=1; dir<STDIM; dir++)
               {
               while(GC->theta[r][dir]>PI2)
                    {
                    GC->theta[r][dir]-=PI2;
                    }
               while(GC->theta[r][dir]<-PI2)
                    {
                    GC->theta[r][dir]+=PI2;
                    }
               }
            }
         }
     #endif
   #endif

   GC->update_index++;
   }


#endif
