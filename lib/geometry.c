#ifndef GEOMETRY_C
#define GEOMETRY_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>

#include"../include/geometry.h"
#include"../include/gparam.h"

// single index = lexicographic index
long (*cart_to_si)(int const * const cartcoord, GParam const * const param)=&cart_to_lex; // cartesian coordinates -> single index
void (*si_to_cart)(int *cartcoord, long si, GParam const * const param)=&lex_to_cart;     // single index -> cartesian coordinates

// initialize geometry
void init_geometry(Geometry *geo, GParam const * const param)
  {
  int i, value, valuep, valuem, err;
  long r, rm, rp;
  int cartcoord[STDIM];

  // allocate memory
  err=posix_memalign((void**)&(geo->d_nnp), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_nnm), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_bcsitep), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(int *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_bcsitem), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(int *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(geo->d_nnp[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     err=posix_memalign((void**)&(geo->d_nnm[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     err=posix_memalign((void**)&(geo->d_bcsitep[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(int));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     err=posix_memalign((void**)&(geo->d_bcsitem[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(int));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // INITIALIZE
  for(r=0; r<param->d_volume; r++)
     {
     si_to_cart(cartcoord, r, param);

     for(i=0; i<STDIM; i++)
        {
        value=cartcoord[i];

        geo->d_bcsitep[r][i]=1;
        geo->d_bcsitem[r][i]=1;

        valuep=value+1;
        if(valuep >= param->d_size[i])
          {
          valuep-=param->d_size[i];
          geo->d_bcsitep[r][i]=-1;
          }
        cartcoord[i]=valuep;
        rp=cart_to_si(cartcoord, param);
        geo->d_nnp[r][i]=rp;

        valuem=value-1;
        if(valuem<0)
          {
          valuem+=param->d_size[i];
          geo->d_bcsitem[r][i]=-1;
          }
        cartcoord[i]=valuem;
        rm=cart_to_si(cartcoord, param);
        geo->d_nnm[r][i]=rm;

        cartcoord[i]=value;
        }

     } // end of loop on r

  #ifdef DEBUG
    test_geometry(geo, param);
  #endif
  }  


// free memory
void free_geometry(Geometry *geo, GParam const * const param)
  {
  long r;

  for(r=0; r<param->d_volume; r++)
     {
     free(geo->d_nnp[r]);
     free(geo->d_nnm[r]);
     free(geo->d_bcsitep[r]);
     free(geo->d_bcsitem[r]);
     }
  free(geo->d_nnp);
  free(geo->d_nnm);
  free(geo->d_bcsitep);
  free(geo->d_bcsitem);
  }


// next neighbour in + direction
long nnp(Geometry const * const geo, long r, int i);


// next neighbour in - direction
long nnm(Geometry const * const geo, long r, int i);


// return -1 if nnp[r][i] is beyond the boundary, else +1
int bcsitep(Geometry const * const geo, long r, int i);


// return -1 if nnm[r][i] is beyond the boundary, else +1
int bcsitem(Geometry const * const geo, long r, int i);


void test_geometry(Geometry const * const geo, GParam const * const param)
  {
  long si, ris_test, si_bis;
  int dir, cart[STDIM];

  // test of lex_to_cart <-> cart_to_lex
  for(si=0; si < param->d_volume; si++)
     {
     lex_to_cart(cart, si, param);
     ris_test=cart_to_lex(cart, param);

     if(si != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of nnp <-> nnm
  for(si=0; si < param->d_volume; si++)
     {
     for(dir=0; dir<STDIM; dir++)
        {
        si_bis=nnp(geo, si, dir);
        ris_test=nnm(geo, si_bis, dir);

        if(si != ris_test)
          {
          fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }
  }



//------------ these are not to be used outside geometry.c ----------------

// cartesian coordinates -> lexicographic index
long cart_to_lex(int const * const cartcoord, GParam const * const param)
  {
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=0; i<STDIM; i++)
     {
     ris+=cartcoord[i]*aux;
     aux*=param->d_size[i];
     }

  // ris = cartcoord[0]
  //      +cartcoord[1]*size[0]
  //      +cartcoord[2]*size[0]*size[1]
  //      +...
  //      +cartcoord[STDIM-1]*size[0]*size[1]*...*size[STDIM-2]

  return ris;
  }


// lexicographic index -> cartesian coordinates
void lex_to_cart(int *cartcoord, long lex, GParam const * const param)
  {
  int i;
  long aux[STDIM];

  aux[0]=1;
  for(i=1; i<STDIM; i++)
     {
     aux[i]=aux[i-1]*param->d_size[i-1];
     }
  // aux[0]=1
  // aux[1]=size[0]
  // aux[2]=size[0]*size[1]
  // ...
  // aux[STDIM-1]=size[0]*size[1]*...*size[STDIM-2]

  for(i=STDIM-1; i>=0; i--)
     {
     cartcoord[i]=(int) (lex/aux[i]);
     lex-=aux[i]*cartcoord[i];
     }
  }


#endif
